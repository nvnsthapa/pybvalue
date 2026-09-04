"""
Step 5 - is the b-value pattern stationary?

Steps 2 and 3 of the Schorlemmer et al. (2004) workflow. Split the catalog into
two abutting periods, map b in each on the same grid, difference them, and ask
at every node whether the change is larger than chance using the Utsu (1992)
test.

The answer at Parkfield is "almost everywhere, yes": the spatial pattern of b
is stationary over decades, which is what licenses using it as a forecast. The
number to reproduce is 34 significant nodes out of 2950, under 1.2%.

Significance follows the paper:
    dAIC < 2  not significant       dAIC > 2  significant (Pb <= 0.05)
                                    dAIC > 5  highly significant (Pb <= 0.01)

Outputs into output/5_stationarity/:
    <name>_stationarity_<split>.png   six panels: b1, b2, db, s1, s2, log Pb
    <name>_stationarity_<split>.json  counts and the largest changes
    <name>_divisionScan.png           significance against where the split falls

Examples
--------
python scripts/5_stationarity.py --prepared parkfield --split 1992-01-01
python scripts/5_stationarity.py --prepared parkfield --split 1996-01-01
python scripts/5_stationarity.py --prepared parkfield --scan
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
from matplotlib.colors import Normalize, TwoSlopeNorm              # noqa: E402

from src import bmap, bplot, fmd, geometry                              # noqa: E402

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUT_DIR = os.path.join(PROJECT_DIR, "output", "2_prepare")
OUTPUT_DIR = os.path.join(PROJECT_DIR, "output", "5_stationarity")

PARKFIELD_FEATURES = {"Middle Mountain asperity": (70.0, 10.0)}

SIGNIFICANT, HIGHLY_SIGNIFICANT = 2.0, 5.0     # thresholds on dAIC


def parse_args():
    parser = argparse.ArgumentParser(
        description="Step 5: test whether the b-value pattern is stationary.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--prepared", required=True)
    parser.add_argument("--split", default="1992-01-01",
                        help="date dividing the two periods (default 1992-01-01)")
    parser.add_argument("--radius", type=float, default=5.0)
    parser.add_argument("--mc", type=float)
    parser.add_argument("--nmin", type=int, default=50)
    parser.add_argument("--spacing", type=float, default=0.5)
    parser.add_argument("--depth", nargs=2, type=float, default=(0.0, 16.0),
                        metavar=("MIN", "MAX"))
    parser.add_argument("--scan", action="store_true",
                        help="also sweep the division date year by year")
    parser.add_argument("--features", action="store_true")
    return parser.parse_args()


def load(name):
    events = pd.read_csv(os.path.join(INPUT_DIR, f"{name}_section.csv"),
                         parse_dates=["time"])
    with open(os.path.join(INPUT_DIR, f"{name}_prepare.json")) as handle:
        meta = json.load(handle)
    if meta.get("geometry") != "section":
        raise SystemExit(f"error: {name} was prepared as '{meta.get('geometry')}'")
    return events, meta


def split_and_compare(events, section, split, args, mc, binsize):
    """Map b either side of `split` and run the Utsu test node by node."""
    moment = pd.Timestamp(split)
    early = events[events["time"] < moment]
    late = events[events["time"] >= moment]
    if min(len(early), len(late)) < args.nmin:
        return None

    common = dict(nmin=args.nmin, binsize=binsize, spacing_km=args.spacing,
                  depth_range=tuple(args.depth))
    first = bmap.map_section(early, section, args.radius, mc, **common)
    second = bmap.map_section(late, section, args.radius, mc, **common)
    comparison = bmap.compare_maps(first, second)
    return early, late, first, second, comparison


def main():
    args = parse_args()
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    name = args.prepared

    events, meta = load(name)
    section = geometry.CrossSection(meta["p1"], meta["p2"], meta["width_km"])
    mc = args.mc if args.mc is not None else meta["mc_used"]
    binsize = meta.get("magnitude_binsize") or None
    regional_b = meta["regional_fit"]["b"]

    print(f"{name}: {len(events):,} events, r = {args.radius} km, Mc = {mc}, "
          f"Nmin = {args.nmin}")

    outcome = split_and_compare(events, section, args.split, args, mc, binsize)
    if outcome is None:
        raise SystemExit(f"error: too few events either side of {args.split}")
    early, late, first, second, comparison = outcome

    label = pd.Timestamp(args.split).strftime("%Y")
    print(f"\nsplit at {args.split}")
    print(f"  period 1  {early['time'].min().date()} .. {early['time'].max().date()}"
          f"  {len(early):,} events   {first['n_resolved']:,} nodes resolved")
    print(f"  period 2  {late['time'].min().date()} .. {late['time'].max().date()}"
          f"  {len(late):,} events   {second['n_resolved']:,} nodes resolved")

    for tag, result in (("period 1", first), ("period 2", second)):
        sigma = result["sigma"][np.isfinite(result["sigma"])]
        if sigma.size:
            print(f"  sigma(b) {tag}: mean {sigma.mean():.3f}, "
                  f"90% below {np.percentile(sigma, 90):.3f}")

    drift = completeness_check(early, late, mc, binsize, comparison)

    compared = comparison["n_compared"]
    significant = comparison["n_significant"]
    highly = comparison["n_highly_significant"]
    share = 100 * significant / compared if compared else float("nan")
    print(f"\n  nodes compared        {compared:,}")
    print(f"  significant           {significant:,}  ({share:.2f}%)   dAIC > 2")
    print(f"  highly significant    {highly:,}   dAIC > 5")
    print("  paper (split 1992):   34 of 2,950  (<1.2%)")

    delta = comparison["delta_b"]
    finite = np.isfinite(delta)
    if finite.any():
        print(f"  delta b range         {np.nanmin(delta):+.3f} .. "
              f"{np.nanmax(delta):+.3f}  (median {np.nanmedian(delta):+.3f})")

    hotspots = largest_changes(first, comparison)
    if hotspots:
        print("\n  largest significant changes:")
        for spot in hotspots:
            print(f"    {spot['along_km']:6.1f} km along, {spot['depth_km']:4.1f} km deep"
                  f"   delta b = {spot['delta_b']:+.2f}   dAIC = {spot['d_aic']:.1f}")

    payload = {
        "prepared": name, "split": args.split, "radius_km": args.radius,
        "mc": mc, "nmin": args.nmin, "regional_b": regional_b,
        "period1": {"start": str(early["time"].min()), "end": str(early["time"].max()),
                    "n_events": len(early), "n_nodes": first["n_resolved"]},
        "period2": {"start": str(late["time"].min()), "end": str(late["time"].max()),
                    "n_events": len(late), "n_nodes": second["n_resolved"]},
        "n_compared": compared, "n_significant": significant,
        "n_highly_significant": highly, "percent_significant": share,
        "completeness_check": drift,
        "largest_changes": hotspots,
    }
    with open(os.path.join(OUTPUT_DIR, f"{name}_stationarity_{label}.json"),
              "w") as handle:
        json.dump(payload, handle, indent=2, default=float)

    figure = plot_panels(name, label, first, second, comparison, regional_b, args)
    print(f"\nwrote {os.path.relpath(figure, PROJECT_DIR)}")

    if args.scan:
        scan_png = run_scan(name, events, section, args, mc, binsize)
        print(f"      {os.path.relpath(scan_png, PROJECT_DIR)}")
    return 0


def completeness_check(early, late, mc, binsize, comparison):
    """
    Warn when the two periods differ in b catalog-wide, not just locally.

    The Utsu test asks whether one node changed. It cannot tell a local
    tectonic change from a network-wide shift in what the catalog records: if
    completeness improves between the periods, *every* node's b moves the same
    way and the test flags the lot. That shows up as a regional b difference of
    the same size and sign as the typical node difference, which is exactly
    what this compares.

    Args:
        early, late (pd.DataFrame): the two period catalogs
        mc (float), binsize (float): completeness and magnitude bin width
        comparison (dict): output of bmap.compare_maps
    Returns:
        dict: per-period regional fits, the shift, and whether it dominates
    """
    first = fmd.fit_gr(early["magnitude"], mc, binsize=binsize)
    second = fmd.fit_gr(late["magnitude"], mc, binsize=binsize)
    shift = second["b"] - first["b"]
    typical = float(np.nanmedian(comparison["delta_b"]))
    spread = float(np.nanmedian(np.abs(comparison["delta_b"])))
    dominated = bool(abs(shift) > 0.5 * spread) if spread > 0 else False

    print(f"\n  regional b period 1   {first['b']:.3f} +- {first['sigma']:.3f} "
          f"(N = {first['n']:,})")
    print(f"  regional b period 2   {second['b']:.3f} +- {second['sigma']:.3f} "
          f"(N = {second['n']:,})")
    print(f"  catalog-wide shift    {shift:+.3f}    "
          f"median node change {typical:+.3f}")

    if dominated:
        print("\n  ! The whole catalog shifts by about as much as a typical node "
              "does.\n"
              "    The Utsu test cannot separate that from local change, so these\n"
              "    counts are inflated by a completeness difference between the\n"
              "    periods rather than by tectonics. Re-run with a higher --mc\n"
              "    and see whether the shift, and the counts, go away.")
    return {"b_period1": first["b"], "b_period2": second["b"],
            "regional_shift": shift, "median_node_change": typical,
            "shift_dominates": dominated}


def largest_changes(first, comparison, limit=4, separation_km=8.0):
    """
    The strongest significant changes, thinned so one broad patch reports once.

    Without the separation filter the top hits are all neighbouring nodes of a
    single anomaly, which reads as four findings instead of one.
    """
    d_aic = comparison["d_aic"]
    candidates = np.argwhere(np.isfinite(d_aic) & (d_aic > SIGNIFICANT))
    if candidates.size == 0:
        return []
    order = np.argsort([-d_aic[tuple(index)] for index in candidates])

    chosen = []
    for position in order:
        row, column = candidates[position]
        along = float(first["grid_along"][row, column])
        depth = float(first["grid_depth"][row, column])
        if any(abs(along - spot["along_km"]) < separation_km for spot in chosen):
            continue
        chosen.append({"along_km": along, "depth_km": depth,
                       "delta_b": float(comparison["delta_b"][row, column]),
                       "d_aic": float(d_aic[row, column]),
                       "p_b": float(comparison["p_b"][row, column])})
        if len(chosen) == limit:
            break
    return chosen


def plot_panels(name, label, first, second, comparison, regional_b, args):
    """The paper's Figure 8 layout: b, b, delta-b, sigma, sigma, significance."""
    fig, axes = plt.subplots(3, 2, figsize=(16, 9.4), sharex=True, sharey=True,
                             constrained_layout=True)

    norm_b = bplot.bvalue_norm(regional_b)
    cmap_b = bplot.bvalue_cmap()
    for ax, result, title in ((axes[0, 0], first, "(a)  b, period 1"),
                              (axes[0, 1], second, "(b)  b, period 2")):
        mesh = bplot.draw_section_map(ax, result, norm_b, cmap_b, key="b",
                                      show_xlabel=False)
        ax.set_title(f"{title}   ({result['n_resolved']:,} nodes)",
                     fontsize=10.5, loc="left")
    bar = fig.colorbar(mesh, ax=axes[0, :].tolist(), fraction=.02, pad=.006,
                       extend="both")
    bar.set_label(f"b-value  (rule = regional {regional_b:.2f})", fontsize=9)
    bar.ax.axhline(regional_b, color="#1a1a1a", lw=1.3)

    # difference: diverging about zero, in a different hue pair from b itself
    delta = np.ma.masked_invalid(comparison["delta_b"])
    reach = float(np.nanpercentile(np.abs(comparison["delta_b"]), 98)) or 0.1
    mesh = axes[1, 0].pcolormesh(first["along"], first["depth"], delta,
                                 cmap=bplot.delta_cmap(),
                                 norm=TwoSlopeNorm(vcenter=0, vmin=-reach, vmax=reach),
                                 shading="nearest", rasterized=True)
    axes[1, 0].set_title("(c)  change in b  (period 2 - period 1)",
                         fontsize=10.5, loc="left")
    bar = fig.colorbar(mesh, ax=axes[1, 0], fraction=.04, pad=.006, extend="both")
    bar.set_label("$\\Delta b$", fontsize=9)

    # significance, with the paper's two thresholds drawn on
    log_p = np.ma.masked_invalid(comparison["log_p"])
    mesh = axes[1, 1].pcolormesh(first["along"], first["depth"], log_p,
                                 cmap=bplot.sequential_cmap("probability").reversed(),
                                 norm=Normalize(-4, 0), shading="nearest",
                                 rasterized=True)
    axes[1, 1].set_title(f"(d)  Utsu test    {comparison['n_significant']:,} of "
                         f"{comparison['n_compared']:,} nodes significant "
                         f"({100*comparison['n_significant']/max(comparison['n_compared'],1):.1f}%)",
                         fontsize=10.5, loc="left")
    filled = np.where(np.isfinite(comparison["d_aic"]), comparison["d_aic"], -9)
    for threshold, style in ((SIGNIFICANT, "-"), (HIGHLY_SIGNIFICANT, "--")):
        axes[1, 1].contour(first["grid_along"], first["grid_depth"], filled,
                           levels=[threshold], colors="#17445f",
                           linewidths=1.3, linestyles=style, zorder=6)
    bar = fig.colorbar(mesh, ax=axes[1, 1], fraction=.04, pad=.006, extend="min")
    bar.set_label("$\\log_{10} P_b$", fontsize=9)
    bar.ax.axhline(-1.3, color="#17445f", lw=1.3)
    bar.ax.axhline(-1.9, color="#17445f", lw=1.3, ls="--")

    sigma_all = np.concatenate([first["sigma"][np.isfinite(first["sigma"])],
                                second["sigma"][np.isfinite(second["sigma"])]])
    top = float(np.percentile(sigma_all, 98)) if sigma_all.size else 0.3
    cmap_sigma = bplot.sequential_cmap("sigma")
    for ax, result, title in ((axes[2, 0], first, "(e)  $\\sigma(b)$, period 1"),
                              (axes[2, 1], second, "(f)  $\\sigma(b)$, period 2")):
        mesh = ax.pcolormesh(result["along"], result["depth"],
                             np.ma.masked_invalid(result["sigma"]), cmap=cmap_sigma,
                             norm=Normalize(0, top), shading="nearest",
                             rasterized=True)
        ax.set_title(title, fontsize=10.5, loc="left")
        ax.set_xlabel("distance along section from P1 (km)")
    bar = fig.colorbar(mesh, ax=axes[2, :].tolist(), fraction=.02, pad=.006,
                       extend="max")
    bar.set_label("$\\sigma(b)$", fontsize=9)

    for ax in axes.ravel():
        ax.set_ylim(first["depth"].max(), first["depth"].min())
        if args.features:
            bplot.mark_features(ax, PARKFIELD_FEATURES)
    for ax in axes[:, 0]:
        ax.set_ylabel("depth (km)")

    fig.suptitle(f"{name}: stationarity of b across a {label} split   "
                 f"(r = {args.radius:.0f} km, Nmin = {args.nmin}; "
                 f"solid contour $\\Delta$AIC = 2, dashed = 5)", fontsize=12.5)
    path = os.path.join(OUTPUT_DIR, f"{name}_stationarity_{label}.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def run_scan(name, events, section, args, mc, binsize):
    """
    Sweep the division date - the paper's "numerous possible catalog divisions".

    A pattern that is genuinely stationary stays stationary wherever the split
    is placed; a spike at one particular date is the interesting case.
    """
    start, end = events["time"].min(), events["time"].max()
    years = range(start.year + 3, end.year - 2)
    print(f"\nscanning division dates {start.year + 3}..{end.year - 3}")

    rows = []
    for year in years:
        outcome = split_and_compare(events, section, f"{year}-01-01", args, mc, binsize)
        if outcome is None:
            continue
        _, _, _, _, comparison = outcome
        compared = comparison["n_compared"]
        if not compared:
            continue
        rows.append({"year": year, "compared": compared,
                     "significant": comparison["n_significant"],
                     "highly": comparison["n_highly_significant"],
                     "percent": 100 * comparison["n_significant"] / compared})
        print(f"  {year}  {rows[-1]['significant']:5,} / {compared:5,} "
              f"= {rows[-1]['percent']:5.2f}%")

    if not rows:
        raise SystemExit("error: no division produced a comparison")
    table = pd.DataFrame(rows)

    fig, ax = plt.subplots(figsize=(10, 4.4))
    ax.plot(table.year, table.percent, "-o", color="#8c2d16", lw=2, ms=6,
            label="significant ($\\Delta$AIC > 2)")
    ax.plot(table.year, 100 * table.highly / table.compared, "-s",
            color="#17445f", lw=2, ms=5,
            label="highly significant ($\\Delta$AIC > 5)")
    ax.axhline(1.2, color="#4a4a4a", ls=":", lw=1.4)
    ax.annotate("paper: <1.2% at the 1992 split", (table.year.iloc[0], 1.2),
                fontsize=9, color="#4a4a4a", textcoords="offset points",
                xytext=(4, 5))
    ax.set_xlabel("year the catalog is divided")
    ax.set_ylabel("nodes with a significant change (%)")
    ax.set_title(f"{name}: stationarity against where the split falls "
                 f"(r = {args.radius:.0f} km)", fontsize=11)
    ax.legend(fontsize=9, framealpha=.9)
    ax.grid(alpha=.3)
    fig.tight_layout()

    path = os.path.join(OUTPUT_DIR, f"{name}_divisionScan.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    table.to_csv(os.path.join(OUTPUT_DIR, f"{name}_divisionScan.csv"), index=False)
    return path


if __name__ == "__main__":
    sys.exit(main())
