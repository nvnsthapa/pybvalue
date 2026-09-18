"""
Step 2a - how many events does a reliable b-value actually need?

Nmin = 50 is Schorlemmer et al. (2004)'s own choice, made on Parkfield's own
catalog: refs/METHOD.md quotes their reasoning as "below that the uncertainty
in b grows rapidly" - a judgement call, not a formula, and one this project
has been carrying over to every new region unexamined. This script re-derives
it for whichever catalog CONFIG points at, two ways:

  1. Bootstrap stability: draw many random subsamples of increasing size N,
     without replacement, from the whole (Mc-cut) catalog - the direct
     empirical analogue of "how uncertain would a real sample of size N be",
     not an assumption that the Shi & Bolt formula's premises hold for this
     catalog. Report both the empirical spread of b across draws and the
     mean analytical sigma(b) at each N, so the two can be checked against
     each other rather than one taken on faith.
  2. Map-view impact: at the radius already chosen in step 3, re-run the
     actual spatial map at a few candidate Nmin values and report what
     coverage costs or buys at each - so the choice is not made on the
     bootstrap curve in isolation, disconnected from what it does to the
     map anyone will actually read.

Outputs into output/<study>/2a_nMinStability/:
    <name>_bootstrap.png   b, sigma(b) (empirical vs analytical, log-log)
                           and relative uncertainty against sample size N
    <name>_bootstrap.json  the numbers behind that figure, plus the
                           recommended N
    <name>_mapImpact.png   coverage and contrast against Nmin, at the
                           chosen radius (only when step 2 used map view)

No command-line flags: edit CONFIG below, then

    python scripts/2a_nMinStability.py

Example - what CONFIG should look like:
    CONFIG["prepared"] = "elsalvador"
    CONFIG["mc"] = 1.42
    CONFIG["radius"] = 12.0
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

from src import bmap, fmd, geometry                                # noqa: E402

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# INPUT_DIR / OUTPUT_DIR are computed in main() from CONFIG['study']:
# output/<study>/2_prepare, output/<study>/2a_nMinStability
INPUT_DIR = OUTPUT_DIR = None

DEFAULT_CANDIDATE_N = (10, 15, 20, 30, 40, 50, 75, 100, 150, 200, 300, 400,
                       500, 750, 1000)
DEFAULT_MAP_NMIN = (20, 30, 50, 75, 100, 150, 200)


# Edit this, then just run the script - no command-line flags.
CONFIG = {
    # output/<study>/... - must match CONFIG['study'] used in step 2 for
    # this prepared catalog.
    "study": "elsalvador",

    "prepared": 'elsalvador',   # required: basename in output/<study>/2_prepare
    "mc": 1.42,                  # completeness; None = whatever step 2 used
    "candidate_n": list(DEFAULT_CANDIDATE_N),   # sample sizes to bootstrap
    "trials": 500,                # bootstrap draws per candidate N
    "target_relative_error": 0.15,  # recommend the smallest N at/below
                                    # sigma(b)/b = this (a common convention;
                                    # not a law - read the curve, not just
                                    # the number it points at)
    "seed": 0,                     # reproducible bootstrap draws

    # map-view impact panel only (skipped for section geometry's own map,
    # since a section has no per-node coverage/contrast trade-off distinct
    # from what step 3 already showed):
    "radius": 12.0,                # the radius already chosen in step 3
    "mc_method": "fixed",          # 'fixed', 'maxcurv' or 'ks'
    "spacing": 0.5,                 # node spacing, km
    "map_nmin_candidates": list(DEFAULT_MAP_NMIN),
}


def main():
    global INPUT_DIR, OUTPUT_DIR
    args = SimpleNamespace(**CONFIG)
    if not args.prepared:
        raise SystemExit("error: set CONFIG['prepared'] at the top of this script.")
    INPUT_DIR = os.path.join(PROJECT_DIR, "output", args.study, "2_prepare")
    OUTPUT_DIR = os.path.join(PROJECT_DIR, "output", args.study, "2a_nMinStability")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    name = args.prepared

    catalog_path = os.path.join(INPUT_DIR, f"{name}_section.csv")
    meta_path = os.path.join(INPUT_DIR, f"{name}_prepare.json")
    if not os.path.exists(catalog_path):
        raise SystemExit(f"error: {catalog_path} not found. Run step 2 first.")

    events = pd.read_csv(catalog_path, parse_dates=["time"])
    with open(meta_path) as handle:
        meta = json.load(handle)
    is_section = meta.get("geometry") == "section"
    mc = args.mc if args.mc is not None else meta["mc_used"]
    binsize = meta.get("magnitude_binsize") or None

    pool = events["magnitude"][events["magnitude"] >= mc].to_numpy(dtype=float)
    print(f"{name}: {len(events):,} events, {pool.size:,} at or above Mc={mc:.2f}")

    full_fit = fmd.fit_gr(pool, mc, binsize=binsize)
    print(f"  full-catalog fit   b = {full_fit['b']:.3f} +- {full_fit['sigma']:.3f}  "
          f"(N = {full_fit['n']:,})")

    bootstrap = run_bootstrap(pool, mc, binsize, args)
    recommended = recommend_n(bootstrap, args.target_relative_error)

    print(f"\n{'N':>6} {'trials':>7} {'median b':>9} {'empirical sig':>14} "
          f"{'analytical sig':>15} {'rel. error':>11}")
    for row in bootstrap:
        mark = "  <- recommended" if row["n"] == recommended else ""
        print(f"{row['n']:>6} {row['trials']:>7} {row['b_median']:>9.3f} "
              f"{row['sigma_empirical']:>14.3f} {row['sigma_analytical_mean']:>15.3f} "
              f"{row['relative_error']:>10.1%}{mark}")

    if recommended is not None:
        print(f"\nrecommended Nmin = {recommended} "
              f"(smallest tested N with sigma(b)/b <= {args.target_relative_error:.0%})")
    else:
        print(f"\nno tested N reached sigma(b)/b <= {args.target_relative_error:.0%}; "
              f"the largest tried was {bootstrap[-1]['n']} at "
              f"{bootstrap[-1]['relative_error']:.1%}")
    print(f"  the paper's Nmin = 50 sits at "
          + (f"{next((r['relative_error'] for r in bootstrap if r['n'] == 50), float('nan')):.1%}"
             if any(r["n"] == 50 for r in bootstrap) else "(not in candidate_n)")
          + " relative error on this catalog")

    boot_png = plot_bootstrap(name, bootstrap, full_fit, recommended, args)
    print(f"\nwrote {os.path.relpath(boot_png, PROJECT_DIR)}")

    payload = {
        "prepared": name, "mc": float(mc), "binsize": binsize,
        "pool_size": int(pool.size), "full_catalog_fit": full_fit,
        "trials_per_n": args.trials, "target_relative_error": args.target_relative_error,
        "recommended_nmin": recommended, "bootstrap": bootstrap,
    }
    json_path = os.path.join(OUTPUT_DIR, f"{name}_bootstrap.json")
    with open(json_path, "w") as handle:
        json.dump(payload, handle, indent=2, default=float)
    print(f"      {os.path.relpath(json_path, PROJECT_DIR)}")

    if not is_section:
        map_png = run_map_impact(name, events, meta, mc, binsize, recommended, args)
        if map_png:
            print(f"      {os.path.relpath(map_png, PROJECT_DIR)}")
    return 0


def run_bootstrap(pool, mc, binsize, args):
    """
    For each candidate N, draw `trials` random subsamples (without
    replacement) of size N from `pool` and fit b in each. Returns a list of
    per-N summary dicts, in ascending N, skipping any N larger than the pool
    (without-replacement sampling cannot exceed it).
    """
    rng = np.random.default_rng(args.seed)
    rows = []
    for n in sorted(set(int(v) for v in args.candidate_n)):
        if n > pool.size:
            print(f"  ! skipping N={n}: only {pool.size:,} events at or above Mc")
            continue
        if n < 2:
            continue
        b_values, sigma_values = [], []
        for _ in range(args.trials):
            sample = rng.choice(pool, size=n, replace=False)
            fit = fmd.fit_gr(sample, mc, binsize=binsize, min_events=2)
            if np.isfinite(fit["b"]):
                b_values.append(fit["b"])
                sigma_values.append(fit["sigma"])
        if len(b_values) < 2:
            continue
        b_values = np.asarray(b_values)
        b_median = float(np.median(b_values))
        sigma_empirical = float(np.std(b_values, ddof=1))
        rows.append({
            "n": n, "trials": len(b_values),
            "b_median": b_median,
            "b_p16": float(np.percentile(b_values, 16)),
            "b_p84": float(np.percentile(b_values, 84)),
            "sigma_empirical": sigma_empirical,
            "sigma_analytical_mean": float(np.mean(sigma_values)),
            "relative_error": sigma_empirical / b_median if b_median else np.nan,
        })
    return rows


def recommend_n(bootstrap, target):
    """Smallest tested N whose empirical relative error is at or below `target`."""
    for row in bootstrap:
        if np.isfinite(row["relative_error"]) and row["relative_error"] <= target:
            return row["n"]
    return None


def plot_bootstrap(name, bootstrap, full_fit, recommended, args):
    """
    Three panels: b against N (with empirical 16-84th spread band), sigma(b)
    against N on log-log axes (empirical vs analytical, plus the b/sqrt(N)
    reference every b-value uncertainty formula reduces to), and relative
    error against N with the target threshold and recommended N marked.
    """
    n = np.array([r["n"] for r in bootstrap], dtype=float)
    b_median = np.array([r["b_median"] for r in bootstrap])
    b_p16 = np.array([r["b_p16"] for r in bootstrap])
    b_p84 = np.array([r["b_p84"] for r in bootstrap])
    sigma_emp = np.array([r["sigma_empirical"] for r in bootstrap])
    sigma_ana = np.array([r["sigma_analytical_mean"] for r in bootstrap])
    rel_err = np.array([r["relative_error"] for r in bootstrap])

    fig, (left, mid, right) = plt.subplots(1, 3, figsize=(18, 5.5))

    left.fill_between(n, b_p16, b_p84, color="steelblue", alpha=.2,
                      label="16th-84th percentile")
    left.plot(n, b_median, "o-", color="steelblue", ms=4, label="median b")
    left.axhline(full_fit["b"], color="crimson", ls="--", lw=1.5,
                 label=f"full catalog b = {full_fit['b']:.3f}")
    left.set_xscale("log")
    left.set_xlabel("N (events per sample)"), left.set_ylabel("b")
    left.set_title("(a)  b-value against sample size", fontsize=11, loc="left")
    left.legend(loc="best", fontsize=8.5, framealpha=.9)
    left.grid(alpha=.3, which="both")

    mid.plot(n, sigma_emp, "o-", color="steelblue", ms=4, label="empirical (bootstrap std)")
    mid.plot(n, sigma_ana, "s--", color="0.5", ms=4, label="analytical (Shi & Bolt, mean)")
    reference = full_fit["b"] / np.sqrt(n)
    mid.plot(n, reference, ":", color="crimson", lw=1.5,
            label=r"$b / \sqrt{N}$ reference")
    mid.set_xscale("log"), mid.set_yscale("log")
    mid.set_xlabel("N (events per sample)"), mid.set_ylabel(r"$\sigma(b)$")
    mid.set_title("(b)  uncertainty against sample size", fontsize=11, loc="left")
    mid.legend(loc="best", fontsize=8.5, framealpha=.9)
    mid.grid(alpha=.3, which="both")

    right.plot(n, rel_err * 100, "o-", color="steelblue", ms=4)
    right.axhline(args.target_relative_error * 100, color="crimson", ls="--", lw=1.5,
                 label=f"target {args.target_relative_error:.0%}")
    if 50 in [r["n"] for r in bootstrap]:
        right.axvline(50, color="0.4", ls=":", lw=1.5, label="paper's Nmin = 50")
    if recommended is not None:
        right.axvline(recommended, color="darkgreen", ls="-", lw=1.5,
                      label=f"recommended Nmin = {recommended}")
    right.set_xscale("log")
    right.set_xlabel("N (events per sample)"), right.set_ylabel(r"$\sigma(b) / b$  (%)")
    right.set_title("(c)  relative error against sample size", fontsize=11, loc="left")
    right.legend(loc="best", fontsize=8.5, framealpha=.9)
    right.grid(alpha=.3, which="both")

    fig.suptitle(f"{name}: how many events does a reliable b-value need?   "
                f"({args.trials} bootstrap draws per N, without replacement)",
                fontsize=12.5)
    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, f"{name}_bootstrap.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def run_map_impact(name, events, meta, mc, binsize, recommended, args):
    """
    Re-run the actual map-view b-value map at a few candidate Nmin values,
    at the radius already chosen in step 3, and report coverage/contrast at
    each - connecting the bootstrap curve back to what it costs on the map
    anyone will actually read, rather than leaving Nmin a number chosen in
    isolation from its own consequence.
    """
    region = geometry.MapRegion(meta["bbox"])
    candidates = sorted(set(int(v) for v in args.map_nmin_candidates) | {50})
    if recommended is not None:
        candidates = sorted(set(candidates) | {recommended})

    print(f"\nmap-view impact at r={args.radius:.0f} km "
         f"(Mc method '{args.mc_method}'):")
    rows = []
    for nmin in candidates:
        result = bmap.map_region(events, region, args.radius, mc=mc, nmin=nmin,
                                 binsize=binsize, spacing_km=args.spacing,
                                 mc_method=args.mc_method, positive=False)
        print(f"  Nmin={nmin:4d}  coverage={result['coverage']*100:5.1f}%  "
             f"contrast={result['b_contrast']:.3f}  b={result['b_min']:.2f}..{result['b_max']:.2f}")
        rows.append({"nmin": nmin, "coverage": result["coverage"],
                    "b_contrast": result["b_contrast"],
                    "n_resolved": result["n_resolved"]})

    fig, (left, right) = plt.subplots(1, 2, figsize=(12, 5.5))
    nmin_vals = [r["nmin"] for r in rows]
    coverage = [r["coverage"] * 100 for r in rows]
    contrast = [r["b_contrast"] for r in rows]

    left.plot(nmin_vals, coverage, "o-", color="steelblue", ms=5)
    left.axvline(50, color="0.4", ls=":", lw=1.5, label="paper's Nmin = 50")
    if recommended is not None:
        left.axvline(recommended, color="darkgreen", ls="-", lw=1.5,
                     label=f"bootstrap-recommended Nmin = {recommended}")
    left.set_xlabel("Nmin"), left.set_ylabel("nodes resolved (%)")
    left.set_title("(a)  coverage against Nmin", fontsize=11, loc="left")
    left.legend(loc="best", fontsize=8.5, framealpha=.9)
    left.grid(alpha=.3)

    right.plot(nmin_vals, contrast, "s-", color="#8c2d16", ms=5)
    right.axvline(50, color="0.4", ls=":", lw=1.5, label="paper's Nmin = 50")
    if recommended is not None:
        right.axvline(recommended, color="darkgreen", ls="-", lw=1.5,
                      label=f"bootstrap-recommended Nmin = {recommended}")
    right.set_xlabel("Nmin"), right.set_ylabel("contrast (p95 - p5 of b)")
    right.set_title("(b)  structure retained against Nmin", fontsize=11, loc="left")
    right.legend(loc="best", fontsize=8.5, framealpha=.9)
    right.grid(alpha=.3)

    fig.suptitle(f"{name}: Nmin's cost on the actual map   r = {args.radius:.0f} km",
                fontsize=12.5)
    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, f"{name}_mapImpact.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


if __name__ == "__main__":
    sys.exit(main())
