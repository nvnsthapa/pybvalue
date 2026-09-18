# Findings — reproducing Schorlemmer et al. (2004) at Parkfield
Status of the `pybvalue` pipeline as of 2026-09-02. All five steps run end to
end and reproduce the source paper's published results. Three methodological
traps were found along the way; two of them would silently corrupt results in
any new region, so they are documented here as much as the reproduction is.

**Source paper.** Schorlemmer, D., S. Wiemer, and M. Wyss (2004), *Earthquake
statistics at Parkfield: 1. Stationarity of b values*, JGR 109, B12307,
doi:10.1029/2004JB003234.

> **Citation correction.** The project description cites paper **2**
> (doi:10.1029/2004JB003235, Schorlemmer, Wiemer, Wyss **& Jackson**), the
> probabilistic-forecasting follow-up. The b-value mapping method, the radius
> figure, and the three-step workflow all come from paper **1**, which has three
> authors and no Jackson. Cite paper 1 for the method.

---
## 1. Reproduction scorecard


Everything below is from a catalog **freshly downloaded from NCEDC**, not from
any file shipped with the project.

| Quantity | This pipeline | Paper | |
|---|---|---|---|
| Cross-section length | 110.96 km | 110 km | ✅ |
| Events, M ≥ 1.3, D ≤ 16 km, 1981–2003 | 3,618 | 3,780 | ✅ −4.3% |
| Regional b-value | **0.903 ± 0.014** | 0.92 | ✅ −1.8% |
| Nodes compared, 1992 split | 1,938–2,799 | 2,950 | ✅ |
| b at Middle Mountain asperity | ≈ 0.43 | ≈ 0.5 | ✅ |
| b in creeping section | up to 1.88 | "up to ≈ 2" | ✅ |
| Optimal sampling radius | 5–6 km | 4–5 km, chose 5 | ✅ |
| r ≥ 10 km destroys contrast | contrast 0.42 → 0.27 | "obscures any contrast" | ✅ |
| Stationarity, safe Mc | 0.07–3.0% significant | < 1.2% | ✅ |
| Timing of the one real change | splits peak 1991–92, collapse by 1993 | "initiates around 1993" | ✅ |

The residual gap in event count and b is **catalog vintage**: NCEDC in 2026 is
not NCSN as it stood in 2004, after two decades of magnitude recomputation,
relocation and event revision. Nothing was tuned to close it.

### The paper's central claim, reproduced

Schorlemmer et al. argue that assuming one regional b puts the M ≥ 6 hazard in
the wrong place. That is falsifiable, and it holds:

| Model | Peak annual P(M ≥ 6) | Location | Interpretation |
|---|---|---|---|
| constant b = 0.90, a varying | 0.00258 /yr | 55.5 km along, 4 km deep | the **creeping** section |
| a and b both varying | **0.03008 /yr** | **71.0 km along**, 13 km deep | **Middle Mountain asperity** |

The paper places the asperity at 70 km from P1; we get 71.0 km. Letting b vary
moves the peak 15.5 km along-section and raises it **11.7×**. Independent check:
both M ≥ 4.5 events in the catalog fall at 67.1 and 72.0 km along — inside the
low-b patch, exactly as the paper reports.

---

## 2. Three traps

### 2.1 A lat/lon box is 63% the wrong earthquakes

The natural first move — a rectangle around Parkfield — gives 13,146 events
against the paper's 3,780 and **b = 0.70 against 0.92**.

Cause: a rectangle over this area is **63% Coalinga/Kettleman events**. The 1983
M6.7 is a blind-thrust sequence in the Coast Ranges, tectonically unrelated to
the San Andreas, and it drags b down to 0.65.

The fix is the paper's actual geometry, which is easy to miss because the
endpoints appear only in running text: **P1 = 36.40°N, 121.00°W → P2 = 35.64°N,
120.20°W**, 5 km wide. Note P1 lies *outside* a plausible guessed box.

Three independent checks confirm the geometry:

| Check | Result |
|---|---|
| Section length | 110.96 km (paper: 110 km) |
| Middle Mountain position | 67 km along (paper: asperity at 70 km) |
| Coalinga M6.7 | 35 km off-section → correctly excluded |

The section cut removes 22,689 events; the depth cut removes **6**. At
Parkfield the depth cut is not a meaningful lever — the geometry does all the
work. That will not be true at the Salton Trough.

### 2.2 The bin-correction trap

USGS and NCEDC return **the same events with the same magnitude types**, yet
give b = 0.831 and 0.903. The reason:

```
NCEDC magnitude == floor(USGS magnitude to 0.1)   for 100.0% of 4,366 matched events
```

NCEDC floors magnitudes onto a 0.1 grid; USGS/ComCat reports full precision.
The `Mc − ΔM/2` correction exists *only* to undo that rounding, so it has to
match the catalog's convention. On one identical set of 2,393 events:

| Magnitudes | binCorrection | b | |
|---|---|---|---|
| NCEDC, binned 0.1 | 0.05 | 0.807 | ✅ |
| USGS, continuous | 0.05 | 0.746 | ❌ 7.5% low |
| USGS, continuous | 0 | 0.817 | ✅ |
| USGS, rebinned to 0.1 | 0.05 | 0.807 | ✅ recovers NCEDC |

A 7.5% systematic error shifts the whole map by more than the per-node σ(b), and
looks exactly like signal. `src/fmd.py` now **detects the magnitude grid from
the data** and sets the correction automatically.

### 2.3 The Utsu test cannot see a completeness change

This is the most consequential finding, and it is a limitation of the published
method applied to a modern catalog — not an implementation bug.

The Utsu test asks whether *one node* changed. It cannot distinguish local
tectonics from a network-wide change in what the catalog records. If
completeness improves between periods, every node's b moves together and the
test flags them all.

At the paper's Mc = 1.3 on today's NCEDC catalog, split at 1992:

| Mc | regional b 1981–92 | regional b 1992–03 | catalog-wide shift | median node Δb | significant |
|---|---|---|---|---|---|
| **1.3** | 0.865 | 0.932 | **+0.068** | **+0.069** | **19.5%** |
| 1.4 | 0.906 | 0.901 | −0.005 | −0.012 | 11.6% |
| 1.5 | 0.898 | 0.891 | −0.007 | −0.036 | 11.0% |
| 1.6 | 0.916 | 0.889 | −0.027 | +0.008 | 3.0% |

At Mc = 1.3 the **median node change equals the catalog-wide shift** — nearly
none of the apparent change is spatial. The 1981–1992 period is not reliably
complete at 1.3 in the current catalog, so it is missing small events, its mean
magnitude is too high, and its b comes out too low. The paper documents this
completeness trend itself (Mc ≈ 1.7 in 1970 → ≈ 1.0 by 1990); it simply was not
a problem for the catalog as it stood in 2004.

> **Read the significance drop carefully.** Raising Mc removes the artefact
> *and* discards events, so part of the 19.5% → 3.0% fall is lost statistical
> power, not artefact removal. The clean, power-independent diagnostic is the
> **regional shift** going from +0.068 to ≈ 0. Diagnose with the shift, never
> with the count.

`scripts/5_stationarity.py` checks this automatically and warns when the
catalog-wide shift exceeds half the median absolute node change.

At a completeness-safe Mc = 1.6 the paper's conclusion reproduces cleanly:

| Split | Compared | Significant | |
|---|---|---|---|
| 1992 | 1,938 | 58 | 2.99% |
| 1996 | 1,442 | 1 | 0.07% |
| paper, 1992 | 2,950 | 34 | < 1.2% |

---

## 3. An unlooked-for result: dating the change

Sweeping the division date year by year was added to implement the paper's
"numerous possible catalog divisions". Most splits sit at 0–1.2% significant —
stationary — with **a sharp isolated peak at 1991–1992 (3.3%, 3.0%) that
collapses to 0.15% by 1993**.

A split placed just *before* a change separates before from after and detects
it; placed *after*, the change falls inside one period and washes out. The peak
therefore brackets a real change to ≈ 1992–93.

The paper's abstract: *"only one patch of radius 5 km showed a significant
increase in b, from below average to above, as a function of time. This change
in b **initiates around 1993** and thus correlates in space and time with a
well-documented episode of creep at depth."*

Recovered independently, from a different catalog, without looking for it.

---

## 4. Radius scan

The paper is explicit that the optimal radius depends on local seismotectonic
fabric and data availability, so it must be re-derived for every region.

| r (km) | coverage | contrast (p95−p5 of b) | corr vs finest | b range |
|---|---|---|---|---|
| 2 | 15.1% | 0.814 | 1.000 | 0.44 – 2.10 |
| 3 | 31.5% | 0.775 | 0.795 | 0.42 – 1.77 |
| 4 | 46.9% | 0.636 | 0.707 | 0.42 – 1.87 |
| **5** | **59.0%** | **0.567** | **0.674** | **0.43 – 1.88** |
| 6 | 69.0% | 0.521 | 0.615 | 0.43 – 1.55 |
| 8 | 84.1% | 0.474 | 0.547 | 0.50 – 1.42 |
| 10 | 93.1% | 0.419 | 0.496 | 0.52 – 1.26 |
| 20 | 100.0% | 0.272 | 0.294 | 0.69 – 1.25 |

**Automatic selection does not recover the paper's 5 km, and that is honest.**
The paper's stated reasoning is that radii 2–5 km give "a nearly identical
pattern". Measured, they do not: correlation with the finest map decays
smoothly (1.000 → 0.795 → 0.707 → 0.674) with no break anywhere, so no
threshold isolates 2–5 km as a plateau. A first attempt using "retain 90% of
max contrast" picked 3 km; a 70% threshold would give exactly 5 km, but that is
reverse-engineering the answer.

The shipped rule is knee-detection on the contrast curve, which gives **6 km**.
Schorlemmer et al. chose 5 km by reading the panels and say so. The suggestion
is where to start looking, not an answer.

---

## 5. Deliberate deviations from the paper

| Topic | Paper | Here | Why |
|---|---|---|---|
| σ(b) | bootstrap | Shi & Bolt | bootstrap not yet implemented; **paper notes Shi & Bolt is biased low** |
| Colour ramp | rainbow | diverging, neutral grey midpoint | rainbow invents boundaries and fails for CVD readers |
| Mc for stationarity | 1.3 | 1.6 recommended | 1.3 is no longer complete for 1981–92 (§2.3) |

The colour scale is diverging because the meaningful quantity is *which side of
the regional b* a node sits on. The midpoint is pinned to the regional value so
"average" reads as absent, and every panel of a small-multiple set shares one
scale.

---

## 6. Known limitations

1. **σ(b) is Shi & Bolt, not bootstrap.** The paper explicitly says Shi & Bolt
   "tends to underestimate the true standard deviation" because it assumes a
   complete catalog and a correctly determined Mc. Our σ(b) ≈ 0.061 against
   their 0.116 is partly this and partly the full-period vs half-period
   comparison. **This is the highest-value thing left to implement.**
2. **Output filenames do not encode Mc.** Re-running step 5 at the same split
   with a different Mc silently overwrites the previous JSON and figure. The
   `_1992` files currently on disk are the Mc = 1.6 run.
3. **The section is a straight chord**, but the San Andreas curves; `along` is
   measured on the chord, not the fault trace. Absorbed at ±2.5 km half-width,
   but anomaly positions are not exactly comparable to the paper's distances.
4. **Nodes are not independent.** At 0.5 km spacing with r = 5 km, neighbouring
   samples share nearly all their events, so significant-node *counts* overstate
   the number of independent findings. The paper has the same property; it
   matters for interpreting percentages, not for the maps.
5. **Steps 3–5 assume a cross section.** `MapRegion` exists in `src/geometry.py`
   and is untested downstream.

---

## 7. Next

- **Bootstrap σ(b)**, resampling with Mc recomputed per draw (limitation 1).
- **Wire `MapRegion` through steps 3–5** for the Salton Trough, which needs map
  view rather than a fault-parallel section.
- **Salton Trough catalog.** Use SCEDC or the Hauksson–Yang–Shearer relocated
  catalog — HYS is dense in the south (123,445 events at M ≥ 1.0, D ≤ 20 km) and
  is the better choice there, though it is unusable at Parkfield, where it
  carries only 7% of the on-section events NCEDC has.
- **Re-derive the radius** for the Salton Trough. Do not carry 5 km over.
- **Run the completeness check before trusting any stationarity result**, for
  the same reason it was needed here.
- **Post-2004 extension.** The catalog runs to 2026-09-01 but the reproduction
  stops at 2003 by design — the paper predates the September 2004 Parkfield
  M6.0. Whether the b-value structure survived the earthquake it forecast is a
  real open question the data now supports.
