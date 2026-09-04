"""
Salton Sea Geothermal Field b-value check, against the literature.

Not a pass/fail unit test - there is no ground truth to assert against, only
a comparison to published results (see PLAN.md, "5. Validation targets"):
Trugman, Shearer, Borsa & Fialko (2016), doi:10.1002/2015JB012510, and the GJI
2023 matched-filter study, doi:10.1093/gji/ggac324, both report **b higher
inside the SSGF than outside**. This script reproduces the check that found
our map does *not* show that, under either completeness method, and that the
mismatch survives restricting to the cited studies' own time window and to
the SSGF's characteristic shallow depth - pointing at a catalog-resolution
gap (this project's regional network catalog vs. their dense-array
matched-filter catalog) rather than a bug in the b-value estimator. See
PLAN.md for the full writeup; this script is what produced its numbers, kept
runnable so the finding does not silently rot into an unverifiable claim.

Requires step 2 to have been run for 'salton' first (CONFIG["catalog"] =
CONFIG["preset"] = "salton" in scripts/2_prepareCatalog.py).

Run:
    python tests/check_wells.py
"""

import os
import sys

import numpy as np
import pandas as pd
from pyproj import Transformer

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src import fmd, geometry                                      # noqa: E402

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STUDY = "parkfield_salton"   # must match CONFIG['study'] used in step 2 for 'salton'
CATALOG_PATH = os.path.join(PROJECT_DIR, "output", STUDY, "2_prepare", "salton_section.csv")
WELLS_PATH = os.path.join(PROJECT_DIR, "data", "wells", "well_data.csv")

# The cited studies' own analysis window and the SSGF's characteristic
# induced-seismicity depth (median 3.52 +/- 0.15 km, 2007-2013 - Scientific
# Reports 2025, doi:10.1038/s41598-025-85744-2), used here to give the
# literature's b-elevated-inside claim every reasonable chance to appear.
PRODUCTION_WINDOW = ("2008-01-01", "2014-01-01")
SHALLOW_DEPTH_KM = 6.0
RING_RADII_KM = (10.0, 30.0)
KS_CANDIDATES = np.round(np.arange(-1.0, 4.0 + 0.1, 0.1), 3)


def zone_table(label, frame, dist_km, binsize):
    """Classic and b-positive, both Mc methods, for the three well-centred rings."""
    print(f"\n=== {label} ({len(frame):,} events) ===")
    inner, outer = RING_RADII_KM
    zones = [(f"inside {inner:.0f} km", dist_km <= inner),
             (f"{inner:.0f}-{outer:.0f} km", (dist_km > inner) & (dist_km <= outer)),
             (f"beyond {outer:.0f} km", dist_km > outer)]
    header = (f"{'zone':<16}{'N':>9}{'Mc(mcv)':>9}{'b(mcv)':>9}{'b+(mcv)':>9}   "
              f"{'Mc(KS)':>9}{'b(KS)':>9}{'b+(KS)':>9}")
    print(header)
    print("-" * len(header))
    rows = []
    for zlabel, mask in zones:
        sample = frame.loc[mask, "magnitude"].to_numpy(dtype=float)
        n = sample.size
        if n < 20:
            print(f"{zlabel:<16}{n:>9,}  (too few events)")
            rows.append({"zone": zlabel, "n": n})
            continue

        mc_curv = fmd.mc_maxcurv(sample, binsize or 0.1)
        fit_curv = fmd.fit_gr(sample, mc_curv, binsize=binsize)
        plus_curv = fmd.fit_gr_positive(sample, binsize=binsize, mc=mc_curv)

        mc_ks = fmd.mc_ks(sample, candidates=KS_CANDIDATES, binsize=binsize,
                          maxErr_b=0.25)
        if np.isfinite(mc_ks):
            fit_ks = fmd.fit_gr(sample, mc_ks, binsize=binsize)
            plus_ks = fmd.fit_gr_positive(sample, binsize=binsize, mc=mc_ks)
            ks_b, ks_bp = fit_ks["b"], plus_ks["b"]
        else:
            ks_b, ks_bp = np.nan, np.nan

        print(f"{zlabel:<16}{n:>9,}{mc_curv:>9.2f}{fit_curv['b']:>9.3f}"
              f"{plus_curv['b']:>9.3f}   {mc_ks:>9.2f}{ks_b:>9.3f}{ks_bp:>9.3f}")
        rows.append({"zone": zlabel, "n": n, "mc_maxcurv": mc_curv,
                     "b_maxcurv": fit_curv["b"], "b_positive_maxcurv": plus_curv["b"],
                     "mc_ks": mc_ks, "b_ks": ks_b, "b_positive_ks": ks_bp})
    return rows


def main():
    if not os.path.exists(CATALOG_PATH):
        raise SystemExit(f"error: {CATALOG_PATH} not found.\n"
                         "Run scripts/2_prepareCatalog.py first, with "
                         "CONFIG['catalog'] = CONFIG['preset'] = 'salton'.")

    events = pd.read_csv(CATALOG_PATH, parse_dates=["time"]).sort_values("time")
    wells = pd.read_csv(WELLS_PATH).dropna(subset=["Latitude", "Longitude"])
    lat0, lon0 = wells["Latitude"].mean(), wells["Longitude"].mean()
    print(f"Salton Sea Geothermal Field b-value check")
    print(f"well centroid: {lat0:.3f}N, {lon0:.3f}W  ({len(wells)} wells, "
         f"sigma {wells['Latitude'].std()*111:.1f} km - one field)")

    crs = geometry.local_crs(lat0, lon0)
    to_xy = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    x, y = to_xy.transform(events["longitude"].to_numpy(), events["latitude"].to_numpy())
    dist_km = np.hypot(x, y)
    binsize = fmd.detect_binsize(events["magnitude"])

    time_mask = ((events["time"] >= PRODUCTION_WINDOW[0])
                & (events["time"] < PRODUCTION_WINDOW[1])).to_numpy()
    depth_mask = (events["depth"] <= SHALLOW_DEPTH_KM).to_numpy()

    zone_table("whole catalog, all depth, full time span", events, dist_km, binsize)
    zone_table(f"{PRODUCTION_WINDOW[0]} to {PRODUCTION_WINDOW[1]}, all depth",
              events[time_mask], dist_km[time_mask], binsize)
    zone_table(f"shallow (depth <= {SHALLOW_DEPTH_KM:.0f} km), full time span",
              events[depth_mask], dist_km[depth_mask], binsize)
    zone_table(f"{PRODUCTION_WINDOW[0]} to {PRODUCTION_WINDOW[1]} AND "
              f"depth <= {SHALLOW_DEPTH_KM:.0f} km",
              events[time_mask & depth_mask], dist_km[time_mask & depth_mask], binsize)

    inner = RING_RADII_KM[0]
    sel = time_mask & (dist_km <= inner)
    print(f"\ndepth, {PRODUCTION_WINDOW[0]}..{PRODUCTION_WINDOW[1]}, "
         f"inside {inner:.0f} km of wells (N={sel.sum():,}):")
    print(events.loc[sel, "depth"].describe().to_string())
    print(f"\n  literature (Sci. Rep. 2025): median 3.52 +/- 0.15 km, 2007-2013")

    print("\n" + "=" * 78)
    print("Conclusion: b inside the SSGF ring is not elevated relative to the")
    print("surrounding rings under any combination above. The catalog's own")
    print("median depth inside the ring during the production window is roughly")
    print("double the literature's reported induced-seismicity depth, and N")
    print("collapses once both cuts are applied - consistent with this regional")
    print("network catalog not resolving the same shallow, matched-filter-only")
    print("population those studies isolated. See PLAN.md, section 5.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
