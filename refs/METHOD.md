# Schorlemmer et al. (2004) — spatial b-value mapping

Method notes extracted from the source paper, to be implemented in an
**area-agnostic** form (no Salton Trough / Parkfield assumptions in the core).

## Provenance — citation correction

The project description cites **Schorlemmer, Wiemer, Wyss & Jackson (2004),
"Earthquake statistics at Parkfield: 2. Probabilistic forecasting and testing"**,
JGR 109, B12308, doi:10.1029/2004JB003235.

The b-value mapping method, the radius scan figure, and the three-step workflow
quoted in the project description all come from the **companion paper 1**:

> Schorlemmer, D., S. Wiemer, and M. Wyss (2004), *Earthquake statistics at
> Parkfield: 1. Stationarity of b values*, J. Geophys. Res., 109, B12307,
> doi:10.1029/2004JB003234.

Note paper 1 has **three** authors (no Jackson). Paper 2 is the forecasting/RELM
testing follow-up. Cite paper 1 for everything below.

Local copy: `refs/Schorlemmer2004_Parkfield1_stationarity_bvalues.pdf`
The gridding technique itself is credited to Wiemer & Wyss (2002), implemented
in ZMAP (Wiemer, 2001).

## The three-step workflow (paper 1, abstract)

1. **Determine the optimal sampling volume** by mapping b over a wide range of
   radii and selecting the *largest* radius that still resolves the b-value
   heterogeneity. Larger radii smooth anomalies away; smaller radii only reduce
   coverage without revealing new detail.
2. **Map the difference in b between two periods**, trying numerous catalog
   divisions in time.
3. **Identify significant changes** with the Utsu (1992) test.

## Estimators

### b-value — maximum likelihood (Aki 1965; Utsu 1965; Bender 1983)

```
b = log10(e) / (M̄ − M_min)
```

`M̄` is the mean magnitude of the sample. Critically, **`M_min = Mc − ΔM/2`** —
the completeness must be corrected by half a bin width to compensate for the
bias from rounding magnitudes into ΔM bins. With ΔM = 0.1 the correction is 0.05.

### Uncertainty — Shi & Bolt (1982)

```
σ(b) = 2.30 · b² · sqrt( Σᵢ (Mᵢ − M̄)² / (n(n−1)) )
```

⚠️ The paper explicitly notes this **underestimates** the true standard
deviation, because it assumes a complete catalog and a correctly determined Mc.
Schorlemmer et al. instead compute σ(b) by **bootstrap** (Schorlemmer et al.,
2003), resampling with Mc recomputed per bootstrap sample. The project
description only specifies Shi & Bolt — see "Deviations" below.

### Utsu (1992) test for two b-values

```
ΔAIC = −2(N₁+N₂)·ln(N₁+N₂)
       + 2N₁·ln(N₁ + N₂·b₁/b₂)
       + 2N₂·ln(N₁·b₂/b₁ + N₂)
       − 2

Pb = exp(−ΔAIC/2 − 2)
```

`Pb` is the probability the two b-values come from the same population.

| ΔAIC | Pb | Interpretation |
|---|---|---|
| < 2 | > 0.05 | not significant (stationary) |
| > 2 | ≤ 0.05 | significant — `log₁₀ Pb ≤ −1.3` |
| > 5 | ≤ 0.01 | highly significant — `log₁₀ Pb ≤ −1.9` |

## Sampling geometry

- Nodes on a regular grid, spacing **0.5 km × 0.5 km**.
- At each node, take all events inside a **cylinder of radius r** centred on the
  node. This is a **fixed-radius** neighbourhood, *not* fixed-N nearest
  neighbours. (Wiemer & Wyss's Mc mapping uses fixed-N; the b-value maps in this
  paper use fixed-r.)
- Parkfield is mapped **in cross section** along the fault, so the cylinders lie
  perpendicular to the cross-sectional plane, with length equal to the section
  width (**5 km**). Concretely, an event belongs to the sample at node
  `(along₀, depth₀)` when both hold:

  ```
  |perpendicular distance from the P1-P2 plane|  ≤ 2.5 km
  sqrt((along − along₀)² + (depth − depth₀)²)    ≤ r
  ```

  The radius acts in the **(along-section, depth)** plane; the 5 km width acts
  across it. Treating r as a 3-D epicentral radius is a different sample and
  gives different maps.

  Verified against the paper: with P1/P2 as given, the section is 111.0 km long
  (paper: 110 km) and Middle Mountain falls at 67 km along (paper places the
  asperity at 70 km from P1). The 1983 Coalinga M6.7 sequence lands 35 km
  off-section and is correctly excluded — a plain lat/lon rectangle over the
  same area is 63% Coalinga events and drops regional b from 0.92 to 0.70.
- Require **N ≥ Nmin events with M ≥ Mc** in the sample. Below Nmin, leave the
  node undefined (masked) — do not compute b. Nmin = **50**, chosen because
  below that the uncertainty in b grows rapidly.

### Radius selection

Radii scanned: **r = 2 km to 20 km**.

- r = 2–5 km give a nearly identical pattern → smaller than 5 km reveals no new
  detail, only reduces coverage.
- r ≥ 6 km begins to blur contrasts and displace anomalies toward the edge of
  the active volume.
- r ≥ 10 km obscures the b-value contrast entirely.
- **Optimal for Parkfield: 4–5 km; the study used r = 5 km.**

> "The selection of the optimal radius for Parkfield is not applicable to other
> areas because it depends on the local seismotectonic fabric and data
> availability."

This sentence is the reason the radius scan must be step 1 of *any* new region —
it is not a parameter to copy from the paper.

## Magnitude of completeness

For the b-value maps, Schorlemmer et al. **do not compute Mc per node**. They
first confirm Mc has no strong spatial variability, then assume a spatially
**homogeneous Mc = 1.3** for stability. Mc is only computed per node inside the
bootstrap, to propagate Mc uncertainty into σ(b); for that they cut the catalog
at M = 0.8 so Mc is free to drop that low in resampled draws.

Time-varying Mc was assessed with sample sizes of 500 events, step 10 events
(Wiemer & Wyss, 2000), plus the GENAS algorithm (Habermann, 1983) for
magnitude-dependent rate changes.

## Parkfield reference values (for reproduction)

| Parameter | Value |
|---|---|
| Catalog | NCSN, 1981–2003 (primary) |
| Secondary catalog | 1967–1981, cut at Mc = 1.7 |
| Independent check | HRSN borehole network, 1987–1998.5 |
| Magnitude binning | ΔM = 0.1 |
| Assumed Mc | 1.3 (homogeneous) |
| Events above Mc | 3780 |
| Overall regional b | 0.92 |
| Depth range | D ≤ 16 km |
| Grid spacing | 0.5 × 0.5 km |
| **Cross section** | **P1 = 36.40°N, 121.00°W → P2 = 35.64°N, 120.20°W** |
| Cross-section width | 5 km (±2.5 km either side of the plane) |
| Radius used | 5 km |
| Nmin | 50 |
| Fault segment length | 110 km (P1→P2) |
| Stationarity divisions | 1981–1992 / 1992–2003, and 1981–1996 / 1996–2003 |
| Nodes computed | 2950 |
| Significant changes | 34 of 2950 (< 1.2%) |
| Typical σ(b) | 0.116 (period 1), 0.122 (period 2) |

Expected reproduction targets: regional b = 0.92; a low-b patch (b ≈ 0.5) at the
Middle Mountain asperity ~70 km from P1; high b (1–2) at the southern end of the
creeping zone, 30–60 km from P1.

## Deviations between the paper and the project description

Worth resolving explicitly before implementing:

| Topic | Paper 1 | Project description |
|---|---|---|
| Mc | homogeneous, assumed (1.3) | per-sample, maximum curvature |
| σ(b) | bootstrap (Shi & Bolt noted as biased low) | Shi & Bolt only |
| Projection | cross section along fault | map view (lat/lon) |
| Citation | paper 1 (003234) | paper 2 (003235) |

Neither project choice is wrong — per-node max-curvature Mc is standard practice
and map view is the natural product for a basin rather than a single fault
strand — but they are departures, and the code should support both sides.

### Decision (2026-09-02)

**Follow the paper.** For the Parkfield reproduction we adopt Schorlemmer's
choices exactly — fixed-radius cylinders, homogeneous assumed Mc, cross-section
geometry, Shi & Bolt σ(b) — and drop the project-description variants for now.
Reproducing published numbers requires reproducing published method; departures
get revisited **before map creation** for the new region, once the code is
validated against a known answer.

Deferred, not discarded: per-node max-curvature Mc, bootstrap σ(b), map-view
geometry, fixed-N neighbourhoods.

## Catalog source selection (measured, Parkfield section, 1981–2003)

| Source | On-section | n(M≥1.3, D≤16) | b | Verdict |
|---|---|---|---|---|
| **NCEDC / NCSN** | 8,621 | **3,618** | **0.903** | correct — the network Schorlemmer used |
| USGS / ComCat | 8,591 | 3,620 | 0.831 | same events, see bin-correction trap below |
| SCEDC | 537 | 461 | — | unusable here: wrong network |
| Hauksson relocated (.mat) | 641 | 566 | — | unusable here: *southern* California catalog |

Reference: n = 3,780, b = 0.92.

SCEDC and the Hauksson–Yang–Shearer relocated catalog are both built on the
*southern* California network. Parkfield sits at their northern edge, so both
carry roughly 7% of the on-section events NCEDC has, with a much higher true Mc.
Neither is wrong as a catalog — HYS is the better choice for the Salton Trough,
where it is dense — but neither can be used at Parkfield.

### The bin-correction trap

USGS and NCEDC return **the same events with the same magnitude types**, but:

```
NCEDC magnitude == floor(USGS magnitude to 0.1)   for 100.0% of 4,366 matched events
```

NCEDC floors magnitudes onto a 0.1 grid; USGS/ComCat reports full precision
(`1.1, 1.2, 1.3…` vs `0.8, 0.81, 0.82, 0.83…`). The `Mc − ΔM/2` correction
exists *only* to compensate for that rounding, so it must match the catalog's
convention. On one identical set of 2,393 events:

| Magnitudes | binCorrection | b |
|---|---|---|
| NCEDC, binned 0.1 | 0.05 | **0.807** ✔ |
| USGS, continuous | 0.05 | 0.746 ✘ — 7.5% low, nothing to correct |
| USGS, continuous | 0 | 0.817 ✔ |
| USGS, rebinned to 0.1 | 0.05 | **0.807** ✔ — recovers NCEDC exactly |

**Rule: check the magnitude grid before choosing binCorrection.** Applying 0.05
to an unbinned catalog biases b by ~7% — a systematic shift of the whole map,
comparable to the σ(b) ≈ 0.116 the paper reports per node, and easily mistaken
for signal. This is why `FMD.fit_GR()` takes `binCorrection` as a parameter and
why `mag_distr.py` passes `binCorrection=0` for its continuous synthetic data.
The paper does the same thing from the other direction: it *rebins* the HRSN
catalog to ΔM = 0.1 before use.

## The stationarity test is confounded by changing completeness

The Utsu test asks whether *one node* changed between two periods. It cannot
distinguish a local tectonic change from a network-wide change in what the
catalog records. If completeness improves between the periods, every node's b
moves the same way and the test flags large numbers of them — correctly, in its
own terms, and misleadingly for the science.

Measured on the modern NCEDC catalog at Parkfield, split at 1992, r = 5 km:

| Mc | regional b, 1981–92 | regional b, 1992–03 | catalog-wide shift | nodes significant |
|---|---|---|---|---|
| **1.3** (the paper's) | 0.865 | 0.932 | **+0.068** | **19.5%** |
| 1.4 | 0.906 | 0.901 | −0.005 | 11.6% |
| 1.5 | 0.898 | 0.891 | −0.007 | 11.0% |
| 1.6 | 0.916 | 0.889 | −0.027 | 3.0% |

At Mc = 1.3 the **median node change (+0.069) is essentially identical to the
catalog-wide shift (+0.068)** — almost all of the apparent "change" is a single
global offset, not spatial structure. It disappears by Mc = 1.4.

The cause is the completeness improvement the paper itself documents (Mc ≈ 1.7
in 1970 falling to ≈ 1.0 by 1990): the 1981–1992 period is not reliably
complete at 1.3 in the current catalog, so it is missing small events, its mean
magnitude is too high, and its b comes out too low.

**Read the drop in significance carefully.** Raising Mc removes the artefact
*and* discards events, so it costs statistical power — some of the fall from
19.5% to 3.0% is simply a smaller sample. The artefact-free evidence is the
**regional shift**, which is power-independent: it goes from +0.068 to ≈0.
Diagnose with the shift, not with the significance count.

`scripts/5_stationarity.py` checks this automatically and warns when the
catalog-wide shift is more than half the median absolute node change.

## Implications for an area-agnostic design

- **Nothing region-specific in the core.** Region bounds, depth cut, time
  window, ΔM, r, Nmin and grid spacing are all inputs.
- **Sampling geometry must be pluggable**: map-view (vertical cylinder, epicentral
  radius) vs. cross-section (horizontal cylinder perpendicular to a
  user-defined P1→P2 plane of given width). Parkfield reproduction needs the
  second; the Salton Trough / Arkansas products need the first.
- **Distances in km, not degrees.** A degree of longitude is ~85 km at 35°N and
  ~94 km at 33°N — a fixed-radius neighbourhood in degrees is not a circle on
  the ground and would silently distort every map. Project to a local metric CRS
  (UTM or azimuthal equidistant centred on the region) before neighbour search.
- **Mc strategy must be a choice**, not a hardcoded assumption: fixed float,
  per-node max curvature, or per-node KS — all three appear across the paper and
  the project description.
- **Nmin masking is part of the result.** Coverage (fraction of nodes resolved)
  is what the radius scan trades against detail, so it must be reported, not
  hidden.

## References

- Aki, K. (1965), *Maximum likelihood estimate of b in the formula log N = a−bM
  and its confidence limits*, Bull. Earthq. Res. Inst. Univ. Tokyo, 43, 237–239.
- Bender, B. (1983), *Maximum likelihood estimation of b values for magnitude
  grouped data*, BSSA, 73(3), 831–851.
- Habermann, R. E. (1983), GENAS — magnitude-dependent rate change detection.
- Schorlemmer, D., S. Wiemer, and M. Wyss (2004), *Earthquake statistics at
  Parkfield: 1. Stationarity of b values*, JGR, 109, B12307,
  doi:10.1029/2004JB003234.
- Schorlemmer, D., S. Wiemer, M. Wyss, and D. D. Jackson (2004), *Earthquake
  statistics at Parkfield: 2. Probabilistic forecasting and testing*, JGR, 109,
  B12308, doi:10.1029/2004JB003235.
- Shi, Y., and B. A. Bolt (1982), *The standard error of the magnitude-frequency
  b value*, BSSA, 72(5), 1677–1687.
- Utsu, T. (1992/1999), Utsu test for difference in b values (AIC-based).
- Wiemer, S. (2001), *A software package to analyze seismicity: ZMAP*, SRL.
- Wiemer, S., and M. Wyss (2000), *Minimum magnitude of completeness in
  earthquake catalogs*, BSSA.
- Wiemer, S., and M. Wyss (2002), *Mapping spatial variability of the
  frequency-magnitude distribution of earthquakes*, Adv. Geophys., 45, 259–302.
