"""
Step 2 - cut a catalog down to the study volume, and check it is usable.

Takes the canonical catalog from step 1, projects it onto a cross section (or a
map region), applies the time / depth / magnitude cuts, and writes the prepared
subset to output/<study>/2_prepare/ (CONFIG['study'] must match the study step
1 wrote the catalog under). Along the way it produces the diagnostics that
justify the completeness choice made in step 3 - the equivalents of Figures 1-3
of Schorlemmer et al. (2004):

    <name>_section.csv     the prepared catalog, with along/perp coordinates
    <name>_map.png         map view: box, section trace, selected vs rejected
    <name>_section.png     the section itself, along-distance against depth
    <name>_completeness.png  Mc against time, both methods (max curvature and
                            KS-distance), and cumulative event count
    <name>_fmd.png         frequency-magnitude distribution at Mc used
    <name>_ks.png           the KS-distance completeness search: FMD fit at
                            the KS-chosen Mc, and KS distance against every
                            candidate tried (compare Goebel's FMD.plotFit /
                            FMD.plotKS)
    <name>_prepare.json    every cut applied, and what it removed

No command-line flags: edit CONFIG below, then

    python scripts/2_prepareCatalog.py

Examples of what CONFIG should look like:

    # Parkfield, reproducing the paper's window and cuts
    CONFIG["catalog"] = "parkfield"
    CONFIG["preset"] = "parkfield"

    # same catalog, no time cut, to look at everything since the paper
    CONFIG["catalog"] = "parkfield"
    CONFIG["preset"] = "parkfield"
    CONFIG["end"] = ""     # "" clears a preset's end cut; None leaves it alone
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

from src import fmd, geometry                                      # noqa: E402

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# INPUT_DIR / OUTPUT_DIR are computed in main() from CONFIG['study']:
# output/<study>/1_catalog, output/<study>/2_prepare

# KS-distance completeness (fmd.ks_distance_curve) is run alongside max
# curvature on every catalog now, as a second opinion, not a replacement -
# see the KS_STEP/KS_MAX_ERR-controlled diagnostic in make_figures(). Same
# defaults as bmap.py's per-node search, for consistency.
KS_STEP = 0.1
KS_MAX_ERR = 0.25
INPUT_DIR = OUTPUT_DIR = None

# Presets mirror the ones in step 1. Values come from refs/METHOD.md.
PRESETS = {
    "parkfield": {
        "geometry": "section",
        "p1": (36.40, -121.00), "p2": (35.64, -120.20), "width_km": 5.0,
        "start": "1981-01-01", "end": "2004-01-01",
        "depth_range": (-5.0, 16.0), "mc": 1.3,
        "note": "Schorlemmer et al. (2004) paper 1. The window stops at 2004 "
                "on purpose: the paper predates the Sept 2004 M6.0.",
    },
    "salton": {
        "geometry": "map",
        "bbox": (32.0, 33.94, -117.16, -114.38),
        "start": "1981-01-01", "end": None,
        "depth_range": (-5.0, 20.0), "mc": None,
        "note": "Salton Trough / Imperial Valley, map view, 20 km depth cut.",
    },
    "elsalvador": {
        "geometry": "map",
        "bbox": (13.65, 14.20, -90.05, -89.30),
        "start": None, "end": None,   # keep the whole ~2.5-year catalog
        "depth_range": (-5.0, 25.0), "mc": None,
        "note": "El Salvador local network, hypoDD/cross-correlation "
                "relocated, 2023-12 to 2026-06. Not a Schorlemmer-style "
                "validation target (no fixed window/depth/Mc from the paper "
                "to reproduce) - map view only, for spatial b mapping. "
                "Swarm-dominated (see step 1's magnitude-through-time panel), "
                "so check the completeness/FMD diagnostics before trusting a "
                "single regional Mc.",
    },
}


# Edit this, then just run the script - no command-line flags. None means
# "use the preset" (or, for CONFIG['name'], "use the catalog name"); an explicit
# overrides the preset. start/end are special: "" (empty string) means
# "clear this cut", where None means "leave whatever the preset set alone".
CONFIG = {
    # output/<study>/... - must match CONFIG['study'] used in step 1 for
    # this catalog.
    "study": "elsalvador",

    "catalog": 'elsalvador',     # required: basename in output/<study>/1_catalog, e.g. "parkfield"
    "preset": 'elsalvador',       # a name from PRESETS above
    "p1": None,            # (lat, lon) - overrides preset geometry, use with p2
    "p2": None,            # (lat, lon)
    "width": None,         # section width, km
    "bbox": None,          # (minlat, maxlat, minlon, maxlon) - map view instead
    "start": None,
    "end": None,
    "depth": None,         # (min, max) km
    "mc": None,            # completeness to report against
    "name": None,          # output basename (default: catalog name)
}


def resolve(args):
    """Merge preset defaults with explicit overrides."""
    cfg = dict(PRESETS.get(args.preset, {}))
    if args.p1 and args.p2:
        cfg.update(geometry="section", p1=tuple(args.p1), p2=tuple(args.p2))
    if args.bbox:
        cfg.update(geometry="map", bbox=tuple(args.bbox))
    if args.width is not None:
        cfg["width_km"] = args.width
    if args.depth:
        cfg["depth_range"] = tuple(args.depth)
    if args.mc is not None:
        cfg["mc"] = args.mc
    for key, value in (("start", args.start), ("end", args.end)):
        if value is not None:
            cfg[key] = value or None      # end="" clears the cut
    if "geometry" not in cfg:
        raise SystemExit("error: set CONFIG['preset'], or CONFIG['p1']/['p2'], "
                         "or CONFIG['bbox']")
    return cfg


def main():
    global INPUT_DIR, OUTPUT_DIR
    args = SimpleNamespace(**CONFIG)
    if not args.catalog:
        raise SystemExit("error: set CONFIG['catalog'] at the top of this script.")
    INPUT_DIR = os.path.join(PROJECT_DIR, "output", args.study, "1_catalog")
    OUTPUT_DIR = os.path.join(PROJECT_DIR, "output", args.study, "2_prepare")
    cfg = resolve(args)
    name = args.name or args.catalog
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    source = os.path.join(INPUT_DIR, f"{args.catalog}.csv")
    if not os.path.exists(source):
        raise SystemExit(f"error: {source} not found. Run step 1 first.")
    events = pd.read_csv(source, parse_dates=["time"])
    print(f"loaded {len(events):,} events from {os.path.relpath(source, PROJECT_DIR)}")

    if cfg.get("note"):
        print(f"  {cfg['note']}")

    # ---- geometry ---------------------------------------------------------
    if cfg["geometry"] == "section":
        region = geometry.CrossSection(cfg["p1"], cfg["p2"], cfg.get("width_km", 5.0))
        print(f"  section {region.length_km:.1f} km long, {region.width_km} km wide, "
              f"strike {region.strike_deg:.1f} deg")
    else:
        region = geometry.MapRegion(cfg["bbox"])
        print(f"  map region {region.bbox}")

    projected = region.project(events)

    # ---- cuts, counted one at a time so nothing vanishes unexplained ------
    cuts, kept = [], projected
    n0 = len(kept)

    def apply(mask, label):
        nonlocal kept
        before = len(kept)
        kept = kept[mask.reindex(kept.index, fill_value=False)]
        cuts.append({"cut": label, "removed": before - len(kept), "kept": len(kept)})
        print(f"  {label:32s} -{before - len(kept):7,}  -> {len(kept):7,}")

    print(f"\ncuts (from {n0:,} events):")
    if cfg.get("start"):
        apply(projected["time"] >= pd.Timestamp(cfg["start"]),
              f"time >= {cfg['start']}")
    if cfg.get("end"):
        apply(projected["time"] < pd.Timestamp(cfg["end"]), f"time <  {cfg['end']}")

    if cfg["geometry"] == "section":
        apply(projected["perp"].abs() <= region.width_km / 2,
              f"|perp| <= {region.width_km / 2} km")
        apply(projected["along"].between(0, region.length_km),
              f"0 <= along <= {region.length_km:.0f} km")
    else:
        apply(projected["latitude"].between(region.bbox[0], region.bbox[1])
              & projected["longitude"].between(region.bbox[2], region.bbox[3]),
              "inside box")

    if cfg.get("depth_range"):
        low, high = cfg["depth_range"]
        apply(projected["depth"].between(low, high), f"depth {low} .. {high} km")

    if kept.empty:
        raise SystemExit("error: no events survived the cuts")

    # ---- completeness -----------------------------------------------------
    binsize = fmd.detect_binsize(kept["magnitude"])
    grid = f"{binsize}" if binsize else "continuous"
    mc_auto = fmd.mc_maxcurv(kept["magnitude"], binsize or 0.1)
    ks_candidates = fmd.ks_candidate_grid(kept["magnitude"], KS_STEP)
    ks_curve = fmd.ks_distance_curve(kept["magnitude"], candidates=ks_candidates,
                                     binsize=binsize, maxErr_b=KS_MAX_ERR)
    mc_ks_auto = ks_curve["best_mc"]
    mc = cfg.get("mc") if cfg.get("mc") is not None else mc_auto
    fit = fmd.fit_gr(kept["magnitude"], mc, binsize=binsize)

    print(f"\nmagnitude grid           {grid}  -> bin correction "
          f"{fmd.bin_correction(kept['magnitude'], binsize):.3f}")
    print(f"Mc (max curvature)       {mc_auto:.2f}")
    print(f"Mc (min KS-distance)     "
          + (f"{mc_ks_auto:.2f}" if np.isfinite(mc_ks_auto) else "none qualified"))
    print(f"Mc used                  {mc:.2f}"
          + ("  (from preset)" if cfg.get("mc") is not None else "  (auto, max curvature)"))
    print(f"N(M >= Mc)               {fit['n']:,}")
    print(f"regional b               {fit['b']:.3f} +- {fit['sigma']:.3f}")
    print(f"regional a               {fit['a']:.2f}")

    # ---- save -------------------------------------------------------------
    prepared = os.path.join(OUTPUT_DIR, f"{name}_section.csv")
    kept.to_csv(prepared, index=False)

    record = {
        "source_catalog": os.path.relpath(source, PROJECT_DIR),
        "geometry": cfg["geometry"], "n_input": n0, "n_prepared": len(kept),
        "cuts": cuts, "magnitude_binsize": binsize,
        "bin_correction": fmd.bin_correction(kept["magnitude"], binsize),
        "mc_maxcurvature": None if np.isnan(mc_auto) else round(float(mc_auto), 3),
        "mc_ks": None if np.isnan(mc_ks_auto) else round(float(mc_ks_auto), 3),
        "mc_used": float(mc), "regional_fit": fit,
    }
    if cfg["geometry"] == "section":
        record.update(p1=list(cfg["p1"]), p2=list(cfg["p2"]),
                      width_km=region.width_km,
                      length_km=round(region.length_km, 3),
                      strike_deg=round(region.strike_deg, 2))
    else:
        record["bbox"] = list(cfg["bbox"])
    for key in ("start", "end", "depth_range"):
        if cfg.get(key):
            record[key] = list(cfg[key]) if key == "depth_range" else str(cfg[key])

    with open(os.path.join(OUTPUT_DIR, f"{name}_prepare.json"), "w") as handle:
        json.dump(record, handle, indent=2, default=str)

    # ---- figures ----------------------------------------------------------
    figures = make_figures(name, events, projected, kept, region, cfg, mc, binsize, fit,
                           mc_ks_auto, ks_curve)

    print(f"\nwrote {os.path.relpath(prepared, PROJECT_DIR)}")
    for path in figures:
        print(f"      {os.path.relpath(path, PROJECT_DIR)}")
    return 0


def make_figures(name, events, projected, kept, region, cfg, mc, binsize, fit,
                 mc_ks_auto, ks_curve):
    """Write the QC figures; returns the paths."""
    written = []
    is_section = cfg["geometry"] == "section"

    def save(fig, suffix):
        path = os.path.join(OUTPUT_DIR, f"{name}_{suffix}.png")
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        written.append(path)

    # --- map view ----------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 8))
    ax.scatter(projected.longitude, projected.latitude, s=1, c="0.8",
               lw=0, label=f"rejected ({len(projected) - len(kept):,})")
    ax.scatter(kept.longitude, kept.latitude, s=2, c="crimson", lw=0,
               label=f"selected ({len(kept):,})")
    if is_section:
        ax.plot([cfg["p1"][1], cfg["p2"][1]], [cfg["p1"][0], cfg["p2"][0]],
                "k-", lw=2, zorder=5)
        for point, label in ((cfg["p1"], "P1"), (cfg["p2"], "P2")):
            ax.plot(point[1], point[0], "ko", ms=7, zorder=6)
            ax.annotate(label, (point[1], point[0]), textcoords="offset points",
                        xytext=(8, 6), fontweight="bold")
    ax.set_xlabel("longitude"), ax.set_ylabel("latitude")
    ax.set_aspect(1 / np.cos(np.radians(float(events.latitude.mean()))))
    ax.set_title(f"{name}: map view")
    ax.legend(loc="upper right", markerscale=6, framealpha=.9)
    ax.grid(alpha=.3)
    save(fig, "map")

    # --- the section itself -------------------------------------------------
    if is_section:
        fig, ax = plt.subplots(figsize=(12, 4))
        scatter = ax.scatter(kept.along, kept.depth, s=3,
                             c=kept.magnitude, cmap="viridis", lw=0)
        plt.colorbar(scatter, ax=ax, label="magnitude", pad=.01)
        ax.set_xlim(0, region.length_km)
        ax.invert_yaxis()
        ax.set_xlabel("distance along section from P1 (km)")
        ax.set_ylabel("depth (km)")
        ax.set_title(f"{name}: cross section, {len(kept):,} events")
        ax.grid(alpha=.3)
        save(fig, "section")

    # --- completeness through time -----------------------------------------
    # Both methods, same window/step: agreement is a check that a swing in
    # Mc is real rather than an artefact of one method's own coarseness (see
    # fmd.mc_vs_time's docstring).
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    centres, mc_series = fmd.mc_vs_time(kept.time.values, kept.magnitude.values,
                                        window=500, step=10, binsize=binsize or 0.1,
                                        method="maxcurv")
    _, mc_ks_series = fmd.mc_vs_time(kept.time.values, kept.magnitude.values,
                                     window=500, step=10, binsize=binsize or 0.1,
                                     method="ks", ks_step=KS_STEP, maxErr_b=KS_MAX_ERR)
    if centres.size:
        top.plot(centres, mc_series, color="0.75", lw=.8)
        smooth = pd.Series(mc_series).rolling(50, center=True, min_periods=1).mean()
        top.plot(centres, smooth, "k-", lw=1.8, label="max curvature (50-window avg)")
        top.plot(centres, mc_ks_series, color="#a8c8dd", lw=.8)
        smooth_ks = pd.Series(mc_ks_series).rolling(50, center=True, min_periods=1).mean()
        top.plot(centres, smooth_ks, color="#2b6a99", lw=1.8, label="KS-distance (50-window avg)")
    top.axhline(mc, color="crimson", ls="--", lw=1.5, label=f"Mc used = {mc:.2f}")
    top.set_ylabel("magnitude of completeness")
    top.legend(loc="upper left", framealpha=.9)
    top.grid(alpha=.3)
    top.set_title(f"{name}: completeness through time, both methods "
                 f"(500-event window, Wiemer & Wyss 2000)")

    bottom.plot(kept.time, np.arange(1, len(kept) + 1), "k-", lw=1.5,
                label="all selected events")
    above = kept[kept.magnitude >= mc]
    twin = bottom.twinx()
    twin.plot(above.time, np.arange(1, len(above) + 1), color="crimson",
              lw=1.5, ls="--")
    twin.set_ylabel(f"cumulative N(M >= {mc:.1f})", color="crimson")
    twin.tick_params(axis="y", colors="crimson")
    bottom.set_ylabel("cumulative N (all)")
    bottom.set_xlabel("time")
    bottom.grid(alpha=.3)
    save(fig, "completeness")

    # --- frequency-magnitude distribution -----------------------------------
    fig, ax = plt.subplots(figsize=(7, 6))
    centres, counts, cumulative = fmd.cumulative_fmd(kept.magnitude, binsize or 0.1)
    ax.semilogy(centres, cumulative, "o", ms=4, color="steelblue",
                label="cumulative  N(>=M)")
    ax.semilogy(centres[counts > 0], counts[counts > 0], "s", ms=4,
                color="0.45", label="binned")
    if np.isfinite(fit["b"]):
        line = np.linspace(mc, kept.magnitude.max() + .2, 20)
        ax.semilogy(line, 10 ** (fit["a"] - fit["b"] * line), "r--", lw=1.8,
                    label=f"$\\log N = {fit['a']:.2f} - {fit['b']:.3f}M$")
    ax.axvline(mc, color="crimson", ls=":", lw=1.5)
    ax.annotate(f"$M_c={mc:.1f}$", (mc, ax.get_ylim()[1]), color="crimson",
                textcoords="offset points", xytext=(4, -14))
    ax.set_xlabel("magnitude"), ax.set_ylabel("number of events")
    ax.set_title(f"{name}: FMD   b = {fit['b']:.3f} $\\pm$ {fit['sigma']:.3f}  "
                 f"(N = {fit['n']:,})")
    ax.legend(loc="upper right", framealpha=.9)
    ax.grid(alpha=.3, which="both")
    save(fig, "fmd")

    # --- KS-distance completeness search ------------------------------------
    # Compare Goebel's FMD.plotFit + FMD.plotKS: the fit at the KS-chosen Mc
    # (top), and the search that chose it (bottom) - not just the winner, so
    # a narrow or noisy minimum is visible rather than hidden behind a single
    # number, and disqualified candidates (sigma(b) >= KS_MAX_ERR) are marked
    # rather than silently absent.
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(6, 10))
    if np.isfinite(mc_ks_auto):
        fit_ks = fmd.fit_gr(kept["magnitude"], mc_ks_auto, binsize=binsize)
        centres, counts, cumulative = fmd.cumulative_fmd(kept.magnitude, binsize or 0.1)
        top.semilogy(centres, cumulative, "o", ms=4, color="steelblue",
                     label="cumulative  N(>=M)")
        top.semilogy(centres[counts > 0], counts[counts > 0], "s", ms=4,
                     color="0.45", label="binned")
        line = np.linspace(mc_ks_auto, kept.magnitude.max() + .2, 20)
        top.semilogy(line, 10 ** (fit_ks["a"] - fit_ks["b"] * line), "r--", lw=1.8,
                     label=f"$\\log N = {fit_ks['a']:.2f} - {fit_ks['b']:.3f}M$")
        top.axvline(mc_ks_auto, color="crimson", ls=":", lw=1.5)
        top.annotate(f"$M_c={mc_ks_auto:.2f}$", (mc_ks_auto, top.get_ylim()[1]),
                    color="crimson", textcoords="offset points", xytext=(4, -14))
        top.set_title(f"{name}: FMD at KS-chosen Mc   b = {fit_ks['b']:.3f} "
                     f"$\\pm$ {fit_ks['sigma']:.3f}  (N = {fit_ks['n']:,})")
        top.legend(loc="upper right", framealpha=.9)
    else:
        top.set_title(f"{name}: no candidate qualified "
                     f"(sigma(b) < {KS_MAX_ERR} throughout)")
    top.set_xlabel("magnitude"), top.set_ylabel("number of events")
    top.grid(alpha=.3, which="both")

    cand, ks_d = ks_curve["candidates"], ks_curve["ks_distance"]
    qualifies = ks_curve["qualifies"]
    if cand.size:
        bottom.plot(cand[~qualifies], ks_d[~qualifies], "x", ms=5, color="0.75",
                   label=f"sigma(b) >= {KS_MAX_ERR} (disqualified)")
        bottom.plot(cand[qualifies], ks_d[qualifies], "o", ms=5, color="steelblue",
                   label="qualifies")
        if np.isfinite(mc_ks_auto):
            bottom.axvline(mc_ks_auto, color="crimson", ls="--", lw=1.5,
                          label=f"chosen Mc = {mc_ks_auto:.2f}")
    bottom.set_xlabel("candidate magnitude of completeness")
    bottom.set_ylabel("KS distance to power law")
    bottom.set_title("search behind the KS-chosen Mc "
                     "(Clauset, Shalizi & Newman 2009)")
    bottom.legend(loc="upper right", framealpha=.9)
    bottom.grid(alpha=.3)
    fig.tight_layout()
    save(fig, "ks")

    return written


if __name__ == "__main__":
    sys.exit(main())
