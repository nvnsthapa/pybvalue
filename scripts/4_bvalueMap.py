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

Outputs into output/<study>/4_bvalueMap/ (CONFIG['study'] must match the
study step 2 wrote the prepared catalog under):
    <name>_bvalueMap.png     section view: b, sigma(b) and N on the section
    <name>_probability.png   annual P(M >= M0): constant b against varying b
    <name>_bvalueMap.json    the numbers, including where each peak falls
    <name>_bvalueMap.npz     the raw grids, for step 5

Map view instead writes <name>_bvalueMap.png as b, b-positive, sigma(b) and
per-node Mc (the four panels needed to read and trust the map), plus - unless
CONFIG['diagnostics_figure'] is False - a second, smaller
<name>_bvalueMap_diagnostics.png with b-positive - b and N per node, the two
supporting panels that are not needed every time but matter when the
difference map is what actually shows a completeness artefact.

No command-line flags: edit CONFIG below, then

    python scripts/4_bvalueMap.py

Example - what CONFIG should look like:
    CONFIG["prepared"] = "parkfield"
    CONFIG["radius"] = 5
"""

import json
import os
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib                                                  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                    # noqa: E402
from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm                   # noqa: E402

from src import bmap, bplot, geometry                              # noqa: E402

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# INPUT_DIR / OUTPUT_DIR are computed in main() from CONFIG['study']:
# output/<study>/2_prepare, output/<study>/4_bvalueMap
INPUT_DIR = OUTPUT_DIR = None

PARKFIELD_FEATURES = {"Middle Mountain asperity": (70.0, 10.0),
                      "creeping section": (63.5, 3.0)}


# Edit this, then just run the script - no command-line flags.
CONFIG = {
    # output/<study>/... - must match CONFIG['study'] used in step 2 for
    # this prepared catalog.
    "study": "elsalvador",

    "prepared": 'elsalvador',          # required: basename in output/<study>/2_prepare, e.g. "parkfield"
    "radius": 12.0,              # sampling radius, km (5 = the paper's choice)
    "mc": 1.42,                 # completeness; None = whatever step 2 used
    "nmin": 50,
    "spacing": 0.5,
    "depth": (0.0, 16.0),        # (min, max) km - section geometry only
    "m0": 6.0,                   # target magnitude for the probability map (section only)
    "big": 4.5,                  # overlay events at or above this magnitude (section only)
    # map view only: 'fixed' (one Mc for the whole map), 'maxcurv' (per-node
    # max curvature) or 'ks' (per-node min-KS-distance completeness)
    "mc_method": "fixed",
    "features": False,          # annotate the Parkfield landmarks (section view only)
    "wells": False,              # overlay data/wells/well_data.csv (map view only;
                                # a no-op anywhere outside Salton regardless)
    "diagnostics_figure": True,  # also write the b+-b / N-events figure (map view only)
}


def main():
    global INPUT_DIR, OUTPUT_DIR
    args = SimpleNamespace(**CONFIG)
    if not args.prepared:
        raise SystemExit("error: set CONFIG['prepared'] at the top of this script.")
    INPUT_DIR = os.path.join(PROJECT_DIR, "output", args.study, "2_prepare")
    OUTPUT_DIR = os.path.join(PROJECT_DIR, "output", args.study, "4_bvalueMap")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    name = args.prepared

    events = pd.read_csv(os.path.join(INPUT_DIR, f"{name}_section.csv"),
                         parse_dates=["time"])
    with open(os.path.join(INPUT_DIR, f"{name}_prepare.json")) as handle:
        meta = json.load(handle)
    if meta.get("geometry") != "section":
        return run_map_view(name, events, meta, args)

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


def run_map_view(name, events, meta, args):
    """
    Map-view b-value mapping, with the completeness diagnostic alongside.

    Produces b and b-positive on the same scale (figure 1), and, unless
    disabled, their difference alongside N per node (figure 2). The
    difference matters even when its figure is off: b-positive is largely
    insensitive to a moving detection threshold, so where it disagrees with
    classic b the sample is incomplete rather than tectonically distinct -
    the region-wide median difference is always printed below for that
    reason, whether or not the map of it gets drawn.
    """
    region = geometry.MapRegion(meta["bbox"])
    binsize = meta.get("magnitude_binsize") or None
    mc = args.mc if args.mc is not None else meta["mc_used"]
    # mc_method changes the map (see PLAN.md, Phase B: maxcurv and KS give
    # materially different per-node Mc and coverage), so it goes in the output
    # name - otherwise a second run under a different method silently
    # overwrites the first, and re-running last year's command for the wrong
    # method looks like it worked while quietly discarding a result.
    out_name = f"{name}_{args.mc_method}"

    print(f"{name}: {len(events):,} events, map view {region.bbox}")
    print(f"  r = {args.radius} km, Mc = {mc} ({args.mc_method}), "
          f"Nmin = {args.nmin}, nodes {args.spacing} km")

    result = bmap.map_region(events, region, args.radius, mc=mc, nmin=args.nmin,
                             binsize=binsize, spacing_km=args.spacing,
                             mc_method=args.mc_method, positive=True)

    print(f"\nresolved {result['n_resolved']:,} of {result['n_nodes']:,} nodes "
          f"({result['coverage']*100:.1f}%)")
    print(f"  b          {result['b_min']:.2f} .. {result['b_max']:.2f}  "
          f"(median {result['b_median']:.3f})")
    plus = result["b_positive"][np.isfinite(result["b_positive"])]
    if plus.size:
        print(f"  b-positive {plus.min():.2f} .. {plus.max():.2f}  "
              f"(median {np.median(plus):.3f})")
    print(f"  median (b+ - b)  {result['median_b_difference']:+.3f}")
    node_mc = result["mc_node"][np.isfinite(result["mc_node"])]
    if node_mc.size:
        print(f"  per-node Mc  {node_mc.min():.2f} .. {node_mc.max():.2f}  "
              f"(median {np.median(node_mc):.2f})")

    gap = abs(result["median_b_difference"])
    if gap > 0.05:
        print(f"\n  ! b and b-positive differ by {gap:.3f} region-wide.\n"
              "    That is a completeness signal, not a tectonic one: the\n"
              "    classic estimator is biased low where the catalog is\n"
              "    incomplete. Read the difference panel, and treat b-positive\n"
              "    as the more trustworthy of the two maps.")

    wells = load_wells() if args.wells else None
    figures = plot_region_panels(out_name, result, region, wells, args)

    np.savez_compressed(
        os.path.join(OUTPUT_DIR, f"{out_name}_bvalueMap.npz"),
        x=result["x"], y=result["y"], b=result["b"], sigma=result["sigma"],
        a=result["a"], n=result["n"], mc_node=result["mc_node"],
        b_positive=result["b_positive"], b_difference=result["b_difference"])
    payload = {k: v for k, v in result.items() if not isinstance(v, np.ndarray)}
    payload.update({"prepared": name, "mc_method": args.mc_method,
                    "geometry": "map", "bbox": list(region.bbox)})
    with open(os.path.join(OUTPUT_DIR, f"{out_name}_bvalueMap.json"), "w") as handle:
        json.dump(payload, handle, indent=2, default=float)

    print("\nwrote " + "\n      ".join(os.path.relpath(f, PROJECT_DIR) for f in figures))
    return 0


def load_wells():
    """Geothermal well locations, if present, for overlay on the maps."""
    path = os.path.join(PROJECT_DIR, "data", "wells", "well_data.csv")
    if not os.path.exists(path):
        return None
    wells = pd.read_csv(path).dropna(subset=["Latitude", "Longitude"])
    return wells if len(wells) else None


def plot_region_panels(name, result, region, wells, args):
    """
    Figure 1 (always): b, b-positive, sigma(b), per-node Mc - the map, a
    robustness check on it, and the two things that say how much to trust it
    at any given node.

    Figure 2 (args.diagnostics_figure): b-positive - b difference, and N per
    node - supporting diagnostics, not needed to read the map every time, but
    kept available since the difference panel is often where a completeness
    artefact (as opposed to a real tectonic signal) actually shows up.

    Returns the list of figure paths written (one or two).
    """
    written = [_plot_main_panels(name, result, region, wells, args)]
    if args.diagnostics_figure:
        written.append(_plot_diagnostic_panels(name, result, region, wells, args))
    return written


def _plot_main_panels(name, result, region, wells, args):
    fig, axes = plt.subplots(2, 2, figsize=(15, 12.5), sharex=True, sharey=True,
                             constrained_layout=True)

    # Centre the diverging scale on the median of the *mapped* values, not on
    # the regional fit: where the catalog is incomplete the regional b is
    # biased, and centring on it would paint almost the whole map "above
    # average" for a reason that has nothing to do with the ground. The
    # range itself is the 5th-95th percentile of b and b-positive combined
    # (both share this scale) - not a fixed +-0.45 guess, which either
    # saturates (real spread wider) or wastes most of the ramp (narrower).
    centre = float(result["b_median"])
    combined = np.concatenate([result["b"][np.isfinite(result["b"])],
                               result["b_positive"][np.isfinite(result["b_positive"])]])
    norm = bplot.bvalue_norm(centre, values=combined)
    cmap = bplot.bvalue_cmap()

    mesh = bplot.draw_region_map(axes[0, 0], result, norm, cmap, key="b",
                                 show_xlabel=False)
    axes[0, 0].set_title(f"(a)  classic b   median {centre:.3f}",
                         fontsize=11, loc="left")

    plus_median = float(np.nanmedian(result["b_positive"]))
    bplot.draw_region_map(axes[0, 1], result, norm, cmap, key="b_positive",
                          show_xlabel=False, show_ylabel=False)
    axes[0, 1].set_title(f"(b)  b-positive   median {plus_median:.3f}",
                         fontsize=11, loc="left")
    bar = fig.colorbar(mesh, ax=axes[0, :].tolist(), fraction=.030, pad=.008,
                       extend="both")
    bar.set_label(f"b-value   (rule = mapped median {centre:.2f})", fontsize=9)
    bar.ax.axhline(centre, color="#1a1a1a", lw=1.4)

    sigma = np.ma.masked_invalid(result["sigma"])
    mesh = axes[1, 0].pcolormesh(result["x"], result["y"], sigma,
                                 cmap=bplot.sequential_cmap("sigma"),
                                 shading="nearest", rasterized=True)
    axes[1, 0].set_aspect("equal")
    finite_sigma = result["sigma"][np.isfinite(result["sigma"])]
    axes[1, 0].set_title(f"(c)  $\\sigma(b)$   median "
                         f"{np.median(finite_sigma):.3f}" if finite_sigma.size
                         else "(c)  $\\sigma(b)$", fontsize=11, loc="left")
    bar = fig.colorbar(mesh, ax=axes[1, 0], fraction=.046, pad=.008)
    bar.set_label("$\\sigma(b)$", fontsize=9)

    node_mc = np.ma.masked_invalid(result["mc_node"])
    mesh = axes[1, 1].pcolormesh(result["x"], result["y"], node_mc,
                                 cmap=bplot.sequential_cmap("probability"),
                                 shading="nearest", rasterized=True)
    axes[1, 1].set_aspect("equal")
    axes[1, 1].set_title(f"(d)  magnitude of completeness per node "
                         f"({args.mc_method})", fontsize=11, loc="left")
    bar = fig.colorbar(mesh, ax=axes[1, 1], fraction=.046, pad=.008)
    bar.set_label("$M_c$", fontsize=9)

    any_wells = False
    for index, ax in enumerate(axes.ravel()):
        any_wells = bplot.overlay_wells(ax, region, wells) or any_wells
        bplot.set_lonlat_ticks(ax, region)
        if index >= 2:
            ax.set_xlabel("longitude")
        if index % 2 == 0:
            ax.set_ylabel("latitude")
    if any_wells:
        axes[0, 0].legend(loc="upper right", fontsize=8.5, framealpha=.92)

    fig.suptitle(f"{name}: b-value in map view   r = {args.radius:.0f} km, "
                 f"Nmin = {args.nmin}, nodes {args.spacing} km   "
                 f"({result['coverage']*100:.0f}% resolved)", fontsize=13)
    path = os.path.join(OUTPUT_DIR, f"{name}_bvalueMap.png")
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    return path


def _plot_diagnostic_panels(name, result, region, wells, args):
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), sharex=True, sharey=True,
                             constrained_layout=True)

    difference = np.ma.masked_invalid(result["b_difference"])
    reach = float(np.nanpercentile(np.abs(result["b_difference"]), 98)) or 0.1
    mesh = axes[0].pcolormesh(result["x"], result["y"], difference,
                              cmap=bplot.delta_cmap(),
                              norm=TwoSlopeNorm(vcenter=0, vmin=-reach, vmax=reach),
                              shading="nearest", rasterized=True)
    axes[0].set_aspect("equal")
    axes[0].set_title(f"(a)  b-positive $-$ b   median "
                      f"{result['median_b_difference']:+.3f}  "
                      f"(large = incomplete, not tectonic)",
                      fontsize=11, loc="left")
    bar = fig.colorbar(mesh, ax=axes[0], fraction=.046, pad=.008, extend="both")
    bar.set_label("$b^+ - b$", fontsize=9)

    counts = np.ma.masked_where(result["n"] == 0, result["n"])
    mesh = axes[1].pcolormesh(result["x"], result["y"], counts,
                              cmap=bplot.sequential_cmap("count"),
                              shading="nearest", rasterized=True)
    axes[1].set_aspect("equal")
    axes[1].set_title(f"(b)  N events above $M_c$   Nmin = {args.nmin}",
                      fontsize=11, loc="left")
    bar = fig.colorbar(mesh, ax=axes[1], fraction=.046, pad=.008)
    bar.set_label("N", fontsize=9)

    for ax in axes:
        bplot.overlay_wells(ax, region, wells)
        bplot.set_lonlat_ticks(ax, region)
        ax.set_xlabel("longitude")
    axes[0].set_ylabel("latitude")

    fig.suptitle(f"{name}: b-value diagnostics   r = {args.radius:.0f} km, "
                 f"Nmin = {args.nmin}, nodes {args.spacing} km", fontsize=13)
    path = os.path.join(OUTPUT_DIR, f"{name}_bvalueMap_diagnostics.png")
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    return path


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
