"""
Step 4 - the b-value map at the chosen sampling radius.

Maps b, its uncertainty and the sample size on the section plane, then turns
the a and b values into the quantity the paper is ultimately after: the annual
probability of a damaging earthquake at every node.

That last figure is the paper's argument in one picture. Using a single
regional b value, the highest M >= 6 probability sits in the *creeping* section
- which contradicts the observations. Letting b vary in space moves the peak
onto the Middle Mountain asperity, where the 1966 Parkfield event nucleated and
where the large events in the catalog actually are.

Outputs into output/4_bvalueMap/:
    <name>_bvalueMap.png     b, sigma(b) and N on the section
    <name>_probability.png   annual P(M >= M0): constant b against varying b
    <name>_bvalueMap.json    the numbers, including where each peak falls
    <name>_bvalueMap.npz     the raw grids, for step 5

Example
-------
python scripts/4_bvalueMap.py --prepared parkfield --radius 5
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
from matplotlib.colors import LogNorm, Normalize                   # noqa: E402

from src import bmap, bplot, geometry                              # noqa: E402

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUT_DIR = os.path.join(PROJECT_DIR, "output", "2_prepare")
OUTPUT_DIR = os.path.join(PROJECT_DIR, "output", "4_bvalueMap")

PARKFIELD_FEATURES = {"Middle Mountain asperity": (70.0, 10.0),
                      "creeping section": (63.5, 3.0)}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Step 4: map b at one radius, and the hazard that follows.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--prepared", required=True)
    parser.add_argument("--radius", type=float, default=5.0,
                        help="sampling radius in km (default 5, the paper's choice)")
    parser.add_argument("--mc", type=float)
    parser.add_argument("--nmin", type=int, default=50)
    parser.add_argument("--spacing", type=float, default=0.5)
    parser.add_argument("--depth", nargs=2, type=float, default=(0.0, 16.0),
                        metavar=("MIN", "MAX"))
    parser.add_argument("--m0", type=float, default=6.0,
                        help="target magnitude for the probability map")
    parser.add_argument("--big", type=float, default=4.5,
                        help="overlay events at or above this magnitude")
    parser.add_argument("--features", action="store_true",
                        help="annotate the Parkfield landmarks")
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    name = args.prepared

    events = pd.read_csv(os.path.join(INPUT_DIR, f"{name}_section.csv"),
                         parse_dates=["time"])
    with open(os.path.join(INPUT_DIR, f"{name}_prepare.json")) as handle:
        meta = json.load(handle)
    if meta.get("geometry") != "section":
        raise SystemExit(f"error: {name} was prepared as '{meta.get('geometry')}'")

    section = geometry.CrossSection(meta["p1"], meta["p2"], meta["width_km"])
    mc = args.mc if args.mc is not None else meta["mc_used"]
    binsize = meta.get("magnitude_binsize") or None
    regional_b = meta["regional_fit"]["b"]

    span = events["time"].max() - events["time"].min()
    duration = span.days / 365.25

    print(f"{name}: {len(events):,} events, {duration:.1f} yr "
          f"({events['time'].min().date()} .. {events['time'].max().date()})")
    print(f"  r = {args.radius} km, Mc = {mc}, Nmin = {args.nmin}, "
          f"nodes {args.spacing} km")
    print(f"  regional b = {regional_b:.3f}\n")

    result = bmap.map_section(events, section, args.radius, mc, nmin=args.nmin,
                              binsize=binsize, spacing_km=args.spacing,
                              depth_range=tuple(args.depth))
    print(f"resolved {result['n_resolved']:,} of {result['n_nodes']:,} nodes "
          f"({result['coverage']*100:.1f}%)")
    print(f"  b        {result['b_min']:.2f} .. {result['b_max']:.2f}  "
          f"(median {result['b_median']:.2f})")
    finite_sigma = result["sigma"][np.isfinite(result["sigma"])]
    print(f"  sigma(b) mean {finite_sigma.mean():.3f}, "
          f"90% below {np.percentile(finite_sigma, 90):.3f}")

    # ---- hazard: does letting b vary move the peak? -----------------------
    a_constant = bmap.constant_b_avalue(result, regional_b)
    _, prob_varying = bmap.recurrence(result["a"], result["b"], args.m0, duration)
    _, prob_constant = bmap.recurrence(a_constant, regional_b, args.m0, duration)
    prob_varying = np.where(np.isfinite(result["b"]), prob_varying, np.nan)
    prob_constant = np.where(np.isfinite(result["b"]), prob_constant, np.nan)

    peaks = {}
    for label, field in (("constant b", prob_constant), ("varying b", prob_varying)):
        flat = np.nanargmax(field)
        row, column = np.unravel_index(flat, field.shape)
        peaks[label] = {
            "along_km": float(result["grid_along"][row, column]),
            "depth_km": float(result["grid_depth"][row, column]),
            "annual_probability": float(field[row, column]),
        }
        print(f"\npeak annual P(M>={args.m0}) with {label}: "
              f"{field[row, column]:.4f} at {peaks[label]['along_km']:.1f} km "
              f"along, {peaks[label]['depth_km']:.1f} km deep")

    big = events[events["magnitude"] >= args.big]
    print(f"\n{len(big)} events with M >= {args.big} on the section")
    if len(big):
        print(f"  along {big['along'].min():.1f} .. {big['along'].max():.1f} km, "
              f"Mmax {big['magnitude'].max():.1f}")

    # ---- save -------------------------------------------------------------
    np.savez_compressed(
        os.path.join(OUTPUT_DIR, f"{name}_bvalueMap.npz"),
        along=result["along"], depth=result["depth"], b=result["b"],
        sigma=result["sigma"], a=result["a"], n=result["n"],
        prob_varying=prob_varying, prob_constant=prob_constant)

    payload = {k: v for k, v in result.items() if not isinstance(v, np.ndarray)}
    payload.update({"prepared": name, "duration_years": round(duration, 2),
                    "regional_b": regional_b, "m0": args.m0,
                    "probability_peaks": peaks,
                    "n_big_events": int(len(big)), "big_threshold": args.big})
    with open(os.path.join(OUTPUT_DIR, f"{name}_bvalueMap.json"), "w") as handle:
        json.dump(payload, handle, indent=2, default=float)

    maps_png = plot_maps(name, result, events, big, regional_b, args)
    prob_png = plot_probability(name, result, prob_constant, prob_varying,
                                big, peaks, regional_b, args)
    print(f"\nwrote {os.path.relpath(maps_png, PROJECT_DIR)}")
    print(f"      {os.path.relpath(prob_png, PROJECT_DIR)}")
    return 0


def plot_maps(name, result, events, big, regional_b, args):
    """b, sigma(b) and N, stacked on a shared section axis."""
    fig, axes = plt.subplots(3, 1, figsize=(13, 9.6), sharex=True,
                             constrained_layout=True)

    norm_b = bplot.bvalue_norm(regional_b)
    mesh = bplot.draw_section_map(axes[0], result, norm_b, bplot.bvalue_cmap(),
                                  key="b", show_xlabel=False)
    bar = fig.colorbar(mesh, ax=axes[0], fraction=.028, pad=.008, extend="both")
    bar.set_label(f"b-value   (rule = regional {regional_b:.2f})", fontsize=9)
    bar.ax.axhline(regional_b, color="#1a1a1a", lw=1.4)
    axes[0].set_title(f"b-value    r = {args.radius:.0f} km, Mc = {result['mc']}, "
                      f"Nmin = {result['nmin']}    "
                      f"{result['coverage']*100:.0f}% of nodes resolved",
                      fontsize=11, loc="left")
    if len(big):
        axes[0].scatter(big["along"], big["depth"], s=90, marker="*",
                        facecolor="#fdfdfd", edgecolor="#1a1a1a", linewidth=1.1,
                        zorder=7, label=f"M $\\geq$ {args.big} ({len(big)})")
        axes[0].legend(loc="lower left", fontsize=8.5, framealpha=.92,
                       markerscale=1.1)

    sigma = np.ma.masked_invalid(result["sigma"])
    top = float(np.nanpercentile(result["sigma"], 98)) if sigma.count() else 1.0
    mesh = axes[1].pcolormesh(result["along"], result["depth"], sigma,
                              cmap=bplot.sequential_cmap("sigma"),
                              norm=Normalize(0, top), shading="nearest",
                              rasterized=True)
    bar = fig.colorbar(mesh, ax=axes[1], fraction=.028, pad=.008, extend="max")
    bar.set_label("$\\sigma(b)$  (Shi & Bolt)", fontsize=9)
    axes[1].set_title("uncertainty in b", fontsize=11, loc="left")

    counts = np.ma.masked_less(result["n"], result["nmin"])
    mesh = axes[2].pcolormesh(result["along"], result["depth"], counts,
                              cmap=bplot.sequential_cmap("count"),
                              shading="nearest", rasterized=True)
    bar = fig.colorbar(mesh, ax=axes[2], fraction=.028, pad=.008)
    bar.set_label(f"events with M $\\geq$ {result['mc']}", fontsize=9)
    axes[2].set_title(f"sample size per node (grey = below Nmin = {result['nmin']})",
                      fontsize=11, loc="left")

    for ax in axes:
        ax.set_ylim(result["depth"].max(), result["depth"].min())
        ax.set_ylabel("depth (km)")
        if args.features:
            bplot.mark_features(ax, PARKFIELD_FEATURES)
    axes[-1].set_xlabel("distance along section from P1 (km)")
    fig.suptitle(f"{name}: b-value mapping on the section", fontsize=13)

    path = os.path.join(OUTPUT_DIR, f"{name}_bvalueMap.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_probability(name, result, prob_constant, prob_varying, big, peaks,
                     regional_b, args):
    """
    The paper's Figure 7: annual P(M >= M0) with a constant b, then with b
    varying in space. Both panels share one log scale - the comparison is the
    whole point, so they cannot have separate scales.
    """
    both = np.concatenate([prob_constant[np.isfinite(prob_constant)],
                           prob_varying[np.isfinite(prob_varying)]])
    both = both[both > 0]
    if both.size == 0:
        raise SystemExit("error: no finite probabilities to plot")
    norm = LogNorm(vmin=max(both.min(), both.max() * 1e-4), vmax=both.max())
    cmap = bplot.sequential_cmap("probability")

    fig, axes = plt.subplots(2, 1, figsize=(13, 7), sharex=True,
                             constrained_layout=True)
    panels = [
        (prob_constant, f"(a)  constant b = {regional_b:.2f}, a varying", "constant b"),
        (prob_varying, "(b)  a and b both varying in space", "varying b"),
    ]
    mesh = None
    for ax, (field, title, key) in zip(axes, panels):
        mesh = ax.pcolormesh(result["along"], result["depth"],
                             np.ma.masked_invalid(field), cmap=cmap, norm=norm,
                             shading="nearest", rasterized=True)
        ax.set_ylim(result["depth"].max(), result["depth"].min())
        ax.set_ylabel("depth (km)")
        ax.set_title(title, fontsize=11, loc="left")
        peak = peaks[key]
        ax.plot(peak["along_km"], peak["depth_km"], marker="X", ms=13,
                mfc="#17445f", mec="white", mew=1.6, zorder=8)
        ax.annotate(f"peak {peak['annual_probability']:.3f}/yr",
                    (peak["along_km"], peak["depth_km"]), fontsize=9,
                    textcoords="offset points", xytext=(12, 2), zorder=8,
                    color="#17445f", fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.22", fc="white",
                              ec="#17445f", alpha=.9))
        if len(big):
            ax.scatter(big["along"], big["depth"], s=90, marker="*",
                       facecolor="#fdfdfd", edgecolor="#1a1a1a", linewidth=1.1,
                       zorder=7, label=f"M $\\geq$ {args.big}")
        if args.features:
            bplot.mark_features(ax, PARKFIELD_FEATURES)
    if len(big):
        axes[0].legend(loc="lower left", fontsize=8.5, framealpha=.92)

    bar = fig.colorbar(mesh, ax=axes.tolist(), fraction=.028, pad=.008)
    bar.set_label(f"annual probability of one or more M $\\geq$ {args.m0}",
                  fontsize=9)
    axes[-1].set_xlabel("distance along section from P1 (km)")
    fig.suptitle(f"{name}: where the hazard peaks depends on whether b is "
                 f"allowed to vary", fontsize=13)

    path = os.path.join(OUTPUT_DIR, f"{name}_probability.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


if __name__ == "__main__":
    sys.exit(main())
