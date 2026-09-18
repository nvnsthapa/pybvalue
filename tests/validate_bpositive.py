"""
Validate b-positive against synthetic catalogs with a known b-value.

Three tests, in increasing order of how much they matter:

  1. Complete catalog       - both estimators must recover the true b.
  2. Step change in Mc      - the case that broke our Parkfield stationarity
                              test. Classic b should be biased; b-positive
                              should not.
  3. Bias and scatter       - repeated trials, to check b-positive is unbiased
                              rather than merely lucky, and to measure the
                              precision it costs.

Also cross-checks `fit_gr` (Aki with the dM/2 correction, what Schorlemmer et
al. use) against `beta_tinti` (the exact binned MLE) so we know the
approximation is not costing anything at dM = 0.1.

Run:
    python tests/validate_bpositive.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib                                                  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                    # noqa: E402

from src import fmd                                                # noqa: E402

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STUDY = "parkfield_salton"   # grouping folder only - this test is synthetic, not region-specific
OUTPUT_DIR = os.path.join(PROJECT_DIR, "output", STUDY, "validation")

TRUE_B = 1.0
BINSIZE = 0.1


def synth(n, b, mc_true=0.0, binsize=BINSIZE, rng=None):
    """Magnitudes from a Gutenberg-Richter distribution, binned, in time order."""
    rng = rng or np.random.default_rng()
    raw = mc_true - np.log10(rng.random(n)) / b
    return np.round(raw / binsize) * binsize


def apply_threshold(magnitudes, threshold):
    """Delete events below a per-event detection threshold, keeping time order."""
    return magnitudes[magnitudes >= threshold - BINSIZE / 2]


def report(label, magnitudes, mc, expect=TRUE_B):
    classic = fmd.fit_gr(magnitudes, mc, binsize=BINSIZE)
    positive = fmd.fit_gr_positive(magnitudes, binsize=BINSIZE, mc=mc)
    print(f"  {label:34s} classic b = {classic['b']:.3f} "
          f"({classic['b'] - expect:+.3f})   "
          f"b+ = {positive['b']:.3f} ({positive['b'] - expect:+.3f})   "
          f"N = {classic['n']:,} / diffs {positive['n']:,}")
    return classic, positive


def test_complete(rng):
    print("\n1. Complete catalog, true b = 1.0, Mc = 1.0")
    magnitudes = apply_threshold(synth(60_000, TRUE_B, rng=rng), 1.0)
    classic, positive = report("complete", magnitudes, 1.0)
    ok = abs(classic["b"] - TRUE_B) < 0.03 and abs(positive["b"] - TRUE_B) < 0.03
    print(f"   -> {'PASS' if ok else 'FAIL'}: both recover the true b")
    return ok


def test_varying_threshold(rng):
    """The case that matters: detection threshold changes partway through."""
    print("\n2. Detection threshold steps 1.0 -> 1.8 halfway through")
    print("   (analysis Mc held at 1.0 - i.e. the analyst has not noticed)")
    first = apply_threshold(synth(40_000, TRUE_B, rng=rng), 1.0)
    second = apply_threshold(synth(40_000, TRUE_B, rng=rng), 1.8)
    magnitudes = np.concatenate([first, second])

    classic, positive = report("stepped threshold", magnitudes, 1.0)
    classic_err = abs(classic["b"] - TRUE_B)
    positive_err = abs(positive["b"] - TRUE_B)
    ok = positive_err < classic_err / 2 and positive_err < 0.10
    print(f"   classic error {classic_err:.3f}, b+ error {positive_err:.3f}"
          f"  ({classic_err / max(positive_err, 1e-9):.1f}x better)")
    print(f"   -> {'PASS' if ok else 'FAIL'}: b+ resists the threshold change")
    return ok, magnitudes


def test_bias(rng, trials=200):
    print(f"\n3. Bias and scatter over {trials} trials "
          f"(threshold steps 1.0 -> 1.6, analysis Mc = 1.0)")
    classic_values, positive_values = [], []
    for _ in range(trials):
        a = apply_threshold(synth(6_000, TRUE_B, rng=rng), 1.0)
        b = apply_threshold(synth(6_000, TRUE_B, rng=rng), 1.6)
        magnitudes = np.concatenate([a, b])
        classic_values.append(fmd.fit_gr(magnitudes, 1.0, binsize=BINSIZE)["b"])
        positive_values.append(
            fmd.fit_gr_positive(magnitudes, binsize=BINSIZE, mc=1.0)["b"])

    classic_values = np.array(classic_values)
    positive_values = np.array(positive_values)
    print(f"  classic    mean {classic_values.mean():.3f}  "
          f"bias {classic_values.mean() - TRUE_B:+.3f}  "
          f"sd {classic_values.std():.3f}")
    print(f"  b-positive mean {positive_values.mean():.3f}  "
          f"bias {positive_values.mean() - TRUE_B:+.3f}  "
          f"sd {positive_values.std():.3f}")
    ok = abs(positive_values.mean() - TRUE_B) < 0.02
    print(f"   -> {'PASS' if ok else 'FAIL'}: b+ is unbiased under a moving "
          f"threshold")
    print(f"   note: b+ scatter is {positive_values.std()/classic_values.std():.1f}x "
          f"the classic one - it uses roughly half the events, so robustness "
          f"costs precision")
    return ok, classic_values, positive_values


def test_tinti(rng):
    print("\n4. Aki + dM/2 (fit_gr) against the exact binned MLE (Tinti & Mulargia)")
    magnitudes = apply_threshold(synth(200_000, TRUE_B, rng=rng), 1.0)
    aki = fmd.fit_gr(magnitudes, 1.0, binsize=BINSIZE)["b"]
    exact = fmd.beta_tinti(magnitudes, 1.0, BINSIZE)
    print(f"  Aki + dM/2 {aki:.4f}   Tinti & Mulargia {exact:.4f}   "
          f"difference {abs(aki - exact):.4f}")
    ok = abs(aki - exact) < 0.005
    print(f"   -> {'PASS' if ok else 'FAIL'}: the approximation costs nothing "
          f"at dM = 0.1")
    return ok


def test_time_order_matters(rng):
    print("\n5. Guard: b-positive on a magnitude-sorted array is meaningless")
    magnitudes = apply_threshold(synth(40_000, TRUE_B, rng=rng), 1.0)
    in_time = fmd.fit_gr_positive(magnitudes, binsize=BINSIZE, mc=1.0)["b"]
    sorted_up = fmd.fit_gr_positive(np.sort(magnitudes), binsize=BINSIZE, mc=1.0)["b"]
    print(f"  time order {in_time:.3f}   magnitude-sorted {sorted_up:.3f}")
    ok = abs(sorted_up - TRUE_B) > 0.2
    print(f"   -> {'PASS' if ok else 'FAIL'}: sorting really does break it, so "
          f"callers must pass time-ordered magnitudes")
    return ok


def make_figure(magnitudes, classic_values, positive_values):
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    fig, (left, right) = plt.subplots(1, 2, figsize=(12.5, 4.6))

    centres, counts, cumulative = fmd.cumulative_fmd(magnitudes, BINSIZE)
    left.semilogy(centres, cumulative, "o", ms=3.5, color="#3d7ea6",
                  label="cumulative  N($\\geq$M)")
    left.semilogy(centres[counts > 0], counts[counts > 0], "s", ms=3.5,
                  color="#8a857d", label="binned")
    left.axvline(1.0, color="#8c2d16", ls=":", lw=1.4)
    left.annotate("analysis $M_c$ = 1.0", (1.0, left.get_ylim()[1]),
                  color="#8c2d16", fontsize=9, textcoords="offset points",
                  xytext=(5, -14))
    left.set_xlabel("magnitude"), left.set_ylabel("number of events")
    left.set_title("Synthetic catalog, threshold stepped mid-way\n"
                   "the kink is the artefact classic b trips on", fontsize=10.5)
    left.legend(fontsize=8.5, framealpha=.9)
    left.grid(alpha=.3, which="both")

    bins = np.linspace(min(classic_values.min(), positive_values.min()) - .02,
                       max(classic_values.max(), positive_values.max()) + .02, 45)
    right.hist(classic_values, bins=bins, color="#c1613c", alpha=.85,
               label=f"classic  (bias {classic_values.mean()-TRUE_B:+.3f})")
    right.hist(positive_values, bins=bins, color="#3d7ea6", alpha=.85,
               label=f"b-positive  (bias {positive_values.mean()-TRUE_B:+.3f})")
    right.axvline(TRUE_B, color="#1a1a1a", lw=2)
    right.annotate("true b = 1.0", (TRUE_B, right.get_ylim()[1] * .96),
                   fontsize=9, textcoords="offset points", xytext=(6, -4))
    right.set_xlabel("estimated b"), right.set_ylabel("trials")
    right.set_title("200 trials under a moving detection threshold\n"
                    "classic is biased low; b-positive is centred", fontsize=10.5)
    right.legend(fontsize=8.5, framealpha=.9, loc="upper left")
    right.grid(alpha=.3)

    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, "bpositive_validation.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def main():
    rng = np.random.default_rng(20260902)
    print("b-positive validation  (van der Elst 2021)")
    print("=" * 66)

    results = [test_complete(rng)]
    ok, magnitudes = test_varying_threshold(rng)
    results.append(ok)
    ok, classic_values, positive_values = test_bias(rng)
    results.append(ok)
    results.append(test_tinti(rng))
    results.append(test_time_order_matters(rng))

    path = make_figure(magnitudes, classic_values, positive_values)
    print("\n" + "=" * 66)
    print(f"{sum(results)}/{len(results)} checks passed")
    print(f"wrote {os.path.relpath(path, PROJECT_DIR)}")
    return 0 if all(results) else 1


if __name__ == "__main__":
    sys.exit(main())
