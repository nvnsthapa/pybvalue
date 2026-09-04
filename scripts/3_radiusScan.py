"""
Step 3 - choose the sampling radius.

This is step 1 of the Schorlemmer et al. (2004) workflow, and it is not
boilerplate: the paper is explicit that the optimal radius depends on the local
seismotectonic fabric and data availability, so it has to be re-derived for
every new region rather than copied from Parkfield.

The rule is "take the largest radius that still resolves the heterogeneity".
Larger radii cover more nodes but smooth anomalies away, so the scan reports
both coverage and contrast at each radius and the choice is made on the two
together.

Outputs into output/3_radiusScan/:
    <name>_radiusScan.png      one panel per radius, shared colour scale
    <name>_radiusTradeoff.png  coverage against contrast
    <name>_radiusScan.json     the numbers behind both figures

Example
-------
python scripts/3_radiusScan.py --prepared parkfield --mc 1.3
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib                                                  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                    # noqa: E402

from src import bmap, bplot, fmd, geometry                         # noqa: E402

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUT_DIR = os.path.join(PROJECT_DIR, "output", "2_prepare")
OUTPUT_DIR = os.path.join(PROJECT_DIR, "output", "3_radiusScan")

# The radii of the paper's Figure 4.
DEFAULT_RADII = (2, 3, 4, 5, 6, 7, 8, 10, 20)

# Landmarks the paper names, for orientation on the Parkfield section.
PARKFIELD_FEATURES = {"Middle Mountain asperity": (70.0, 10.0),
                      "creeping section": (63.5, 3.0)}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Step 3: map b over a range of radii and pick one.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--prepared", required=True,
                        help="basename in output/2_prepare (e.g. parkfield)")
    parser.add_argument("--radii", nargs="+", type=float, default=list(DEFAULT_RADII))
    parser.add_argument("--mc", type=float,
                        help="completeness (default: whatever step 2 used)")
    parser.add_argument("--nmin", type=int, default=50,
                        help="minimum events above Mc per node (default 50)")
    parser.add_argument("--spacing", type=float, default=0.5,
                        help="node spacing, km (default 0.5)")
    parser.add_argument("--depth", nargs=2, type=float, default=(0.0, 16.0),
                        metavar=("MIN", "MAX"))
    parser.add_argument("--features", action="store_true",
                        help="annotate the Parkfield landmarks")
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    name = args.prepared

    catalog_path = os.path.join(INPUT_DIR, f"{name}_section.csv")
    meta_path = os.path.join(INPUT_DIR, f"{name}_prepare.json")
    if not os.path.exists(catalog_path):
        raise SystemExit(f"error: {catalog_path} not found. Run step 2 first.")

    events = pd.read_csv(catalog_path, parse_dates=["time"])
    with open(meta_path) as handle:
        meta = json.load(handle)
    if meta.get("geometry") != "section":
        raise SystemExit("error: this step maps cross sections; "
                         f"{name} was prepared as '{meta.get('geometry')}'")

    section = geometry.CrossSection(meta["p1"], meta["p2"], meta["width_km"])
    mc = args.mc if args.mc is not None else meta["mc_used"]
    binsize = meta.get("magnitude_binsize") or None
    regional = meta.get("regional_fit", {})
    centre = regional.get("b") or 1.0

    print(f"{name}: {len(events):,} events on a {section.length_km:.1f} km section")
    print(f"  Mc={mc}  Nmin={args.nmin}  spacing={args.spacing} km  "
          f"depth {args.depth[0]}..{args.depth[1]} km")
    print(f"  regional b = {centre:.3f} (colour scale centres here)\n")

    results = bmap.radius_scan(events, section, args.radii, mc, nmin=args.nmin,
                               binsize=binsize, spacing_km=args.spacing,
                               depth_range=tuple(args.depth))

    # ---- pick a radius ----------------------------------------------------
    usable = [r for r in results if r["n_resolved"] > 0]
    if not usable:
        raise SystemExit("error: no radius resolved any node; lower --nmin "
                         "or check Mc")
    similarity = pattern_similarity(results)
    for result, value in zip(results, similarity):
        result["similarity_to_finest"] = value

    best = choose_radius(usable)
    print(f"\n{'r (km)':>7} {'coverage':>9} {'contrast':>9} {'corr vs finest':>15}")
    for result, value in zip(results, similarity):
        mark = "  <- knee" if result is best else ""
        print(f"{result['radius_km']:>7.0f} {result['coverage']*100:>8.1f}% "
              f"{result['b_contrast']:>9.3f} {value:>15.3f}{mark}")
    print(f"\nsuggested radius: {best['radius_km']:.0f} km "
          f"(coverage {best['coverage']*100:.0f}%, contrast {best['b_contrast']:.3f})")
    print("  the paper chose 5 km for Parkfield, calling 4-5 km optimal.")
    print("  Similarity decays smoothly with no sharp break, so no automatic")
    print("  rule recovers that exactly - read the panels before committing.")

    # ---- figures ----------------------------------------------------------
    scan_png = plot_scan(name, results, centre, section, args)
    trade_png = plot_tradeoff(name, usable, best)

    payload = {
        "prepared": name, "mc": mc, "nmin": args.nmin,
        "spacing_km": args.spacing, "depth_range": list(args.depth),
        "regional_b": centre, "suggested_radius_km": best["radius_km"],
        "radii": [{k: v for k, v in r.items()
                   if not isinstance(v, np.ndarray)} for r in results],
    }
    with open(os.path.join(OUTPUT_DIR, f"{name}_radiusScan.json"), "w") as handle:
        json.dump(payload, handle, indent=2, default=float)

    print(f"\nwrote {os.path.relpath(scan_png, PROJECT_DIR)}")
    print(f"      {os.path.relpath(trade_png, PROJECT_DIR)}")
    print(f"      {os.path.relpath(os.path.join(OUTPUT_DIR, name + '_radiusScan.json'), PROJECT_DIR)}")
    return 0


def pattern_similarity(results):
    """
    Correlate every map against the finest one, over nodes resolved in both.

    This is the quantity behind the paper's argument that radii of 2-5 km give
    "a nearly identical pattern": if a larger radius still correlates with the
    finest map, it is not destroying structure, only filling in coverage.
    """
    reference = results[0]["b"]
    similarity = []
    for result in results:
        both = np.isfinite(reference) & np.isfinite(result["b"])
        similarity.append(float(np.corrcoef(reference[both], result["b"][both])[0, 1])
                          if both.sum() > 10 else np.nan)
    return similarity


def choose_radius(results):
    """
    The knee of the contrast-against-radius curve.

    Contrast falls steeply while the radius is still small enough to merge
    genuinely different volumes, then flattens once the map is simply smooth.
    The knee - the point furthest from the chord joining the first and last
    radius - is where that transition happens, and is the largest radius still
    buying resolution.

    This is deliberately *not* tuned to reproduce the paper's 5 km. The
    correlation between maps decays smoothly with no sharp break, so no
    automatic rule recovers their choice exactly; they picked it by looking at
    the panels and said so. Treat this as where to start looking, not an answer.
    """
    if len(results) < 3:
        return results[-1]
    radii = np.array([r["radius_km"] for r in results], dtype=float)
    contrast = np.array([r["b_contrast"] for r in results], dtype=float)

    x = (radii - radii[0]) / (radii[-1] - radii[0])
    span = contrast[0] - contrast[-1]
    if span <= 0:
        return results[-1]
    y = (contrast - contrast[-1]) / span
    # distance from the straight line joining the first and last point
    return results[int(np.argmax(np.abs(x + y - 1)))]


def plot_scan(name, results, centre, section, args):
    """Small multiples: one panel per radius, one shared colour scale."""
    columns = 3
    rows = int(np.ceil(len(results) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(5.4 * columns, 2.5 * rows),
                             sharex=True, sharey=True, constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()
    norm = bplot.bvalue_norm(centre)
    cmap = bplot.bvalue_cmap()
    mesh = None

    for index, (ax, result) in enumerate(zip(axes, results)):
        mesh = bplot.draw_section_map(
            ax, result, norm, cmap,
            show_xlabel=index >= len(results) - columns,
            show_ylabel=index % columns == 0)
        bplot.annotate_panel(ax, f"({chr(97+index)})  r = {result['radius_km']:.0f} km",
                             coverage=result["coverage"])
        if args.features:
            bplot.mark_features(ax, PARKFIELD_FEATURES)

    for ax in axes[len(results):]:
        ax.set_visible(False)

    bar = fig.colorbar(mesh, ax=axes.tolist(), fraction=.02, pad=.01,
                       extend="both")
    # The regional value is marked by the rule on the bar; naming it in the
    # label keeps it off the tick column, where it collided with the ticks.
    bar.set_label(f"b-value   (rule = regional {centre:.2f})", fontsize=10)
    bar.ax.axhline(centre, color="#1a1a1a", lw=1.4)

    fig.suptitle(f"{name}: b-value against sampling radius   "
                 f"(Mc = {results[0]['mc']}, Nmin = {results[0]['nmin']}, "
                 f"nodes {args.spacing} km)   grey = fewer than Nmin events",
                 fontsize=12, y=1.02)
    path = os.path.join(OUTPUT_DIR, f"{name}_radiusScan.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_tradeoff(name, results, best):
    """
    What the choice of radius is actually made on.

    Coverage is a percentage and contrast is in b units, so they get separate
    panels rather than being crushed onto one axis. Contrast and similarity do
    share a panel: both are dimensionless 0-1 measures of how much structure
    survives, so they are directly comparable.
    """
    radii = [r["radius_km"] for r in results]
    coverage = [r["coverage"] * 100 for r in results]
    contrast = [r["b_contrast"] for r in results]
    similarity = [r.get("similarity_to_finest", np.nan) for r in results]

    fig, (one, two, three) = plt.subplots(1, 3, figsize=(15, 4.3))

    def mark_best(ax):
        ax.axvline(best["radius_km"], color="#4a4a4a", ls="--", lw=1.3, zorder=1)

    one.plot(radii, coverage, "-o", color="#3d7ea6", lw=2, ms=7)
    mark_best(one)
    one.set_xlabel("sampling radius (km)")
    one.set_ylabel("nodes resolved (%)")
    one.set_title("coverage rises with radius", fontsize=11)
    one.set_ylim(0, 104)
    one.grid(alpha=.3)
    one.annotate(f"knee at r = {best['radius_km']:.0f} km",
                 (best["radius_km"], 8), fontsize=9, color="#4a4a4a",
                 textcoords="offset points", xytext=(7, 0))

    two.plot(radii, contrast, "-s", color="#8c2d16", lw=2, ms=7,
             label="contrast (p95 - p5 of b)")
    two.plot(radii, similarity, "-^", color="#3f6b34", lw=2, ms=7,
             label="correlation with finest map")
    mark_best(two)
    two.set_xlabel("sampling radius (km)")
    two.set_ylabel("structure retained (dimensionless)")
    two.set_title("structure is smoothed away", fontsize=11)
    two.set_ylim(0, 1.05)
    two.legend(fontsize=9, framealpha=.9, loc="upper right")
    two.grid(alpha=.3)

    three.plot(coverage, contrast, "-", color="#b9b6b1", lw=1.4, zorder=1)
    three.scatter(coverage, contrast, c=radii, cmap="Blues", s=110, zorder=2,
                  edgecolor="#4a4a4a", linewidth=.9, vmin=-6)
    for radius, cover, contra in zip(radii, coverage, contrast):
        three.annotate(f"{radius:.0f} km", (cover, contra), fontsize=8.5,
                       textcoords="offset points", xytext=(7, 5), color="#4a4a4a")
    three.set_xlabel("nodes resolved (%)")
    three.set_ylabel("contrast (p95 - p5 of b)")
    three.set_title("the trade-off itself", fontsize=11)
    three.grid(alpha=.3)

    fig.suptitle(f"{name}: choosing the sampling radius   "
                 f"(the paper chose 5 km at Parkfield)", fontsize=12)
    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, f"{name}_radiusTradeoff.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


if __name__ == "__main__":
    sys.exit(main())
