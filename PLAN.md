# Plan — spatial b-value map of the Salton Trough / Imperial Valley

Written 2026-09-02, after the Parkfield reproduction (see `output_findings.md`)
and a literature review. This is the project's science target; Parkfield was the
validation case.

## Objective

A spatial map of the Gutenberg-Richter b-value across the Salton Trough and
Imperial Valley, with defensible uncertainty and defensible significance, using
the Schorlemmer et al. (2004) framework updated where the last two decades of
literature say it needs updating.

## Decisions taken

| Decision | Choice |
|---|---|
| Estimator | **Classic Aki MLE *and* b-positive at every node**, compared |
| Sampling | **Fixed radius**, as validated at Parkfield; radius re-derived by scan |
| Catalog | **Hauksson relocated (GrowClust), 1981–2025** |
| Geometry | Map view (`MapRegion`), depth 0–20 km per the project description |
| Declustering | **No** — see §3.3 |
| Completeness (Mc) method | **No single default** — max curvature and KS-distance both computed (`CONFIG['mc_method'] = "maxcurv"\|"ks"`); neither wins outright (§Phase B), read the b⁺−b map to judge which to trust per node |

**Output layout (2026-09-03):** every script's outputs now live under
`output/<study>/<step>/`, where `CONFIG['study']` (default
`"parkfield_salton"`) groups a set of related regions under one tree. This
made room for a second, independent study (El Salvador, catalog already in
`data/Elsalvador/`) without its outputs mixing into these Parkfield/Salton
ones — set `CONFIG['study'] = "elsalvador"` (or similar) in every step of
that run. All filenames and reference numbers elsewhere in this document are
otherwise unchanged, just nested one level deeper on disk.

## 1. The catalog is in hand and already parses

`sc_1981_2025_1d_3d_gc_soda_noqb_10_1_SCSN.gc` from the SCEDC alternative
catalogs (152 MB, updated 2026-02-27) is the current Hauksson–Yang–Shearer
relocated catalog. It is the 29-column GrowClust format that
`src/catalog.py::_read_growclust` already handles — it parsed unmodified.

```
https://scedc.caltech.edu/ftp/catalogs/hauksson/Socal_DD/sc_1981_2025_1d_3d_gc_soda_noqb_10_1_SCSN.gc
```

**Implementation note:** the SCEDC file server returns an empty body to a
default curl/urllib request. It needs a browser `User-Agent`. Step 1's
downloader must set one when a fetch-by-URL path is added.

| | |
|---|---|
| Whole catalog | 808,517 events, 1981-01-01 → 2025-12-31 |
| Salton Trough box (32.0–33.94 N, −117.16 – −114.38 E) | 312,501 |
| … and D ≤ 20 km | **311,756** |

That is **86× the Parkfield on-section sample** and 2.5× what the 1981–2011
`.mat` gives. Sample size is not a constraint here; it is the reason small
radii become viable and the reason significance tests need care.

## 2. Two things the data said immediately

### 2.1 The magnitude grid is 0.01, not 0.1

`detect_binsize` returns **0.01** for this catalog — the relocated catalog
carries recomputed magnitudes at higher precision than the network catalog.
The bin correction is therefore **0.005, not 0.05**.

Carrying Parkfield's 0.05 over would have introduced exactly the error
documented in `refs/METHOD.md` §"bin-correction trap". The automatic detection
catches it. This is the second time that check has paid for itself.

### 2.2 b drifts strongly with Mc — the central problem

| Mc | N | b |
|---|---|---|
| 0.5 | 275,238 | 0.497 ± 0.001 |
| 1.0 | 186,220 | 0.644 ± 0.001 |
| 1.5 | 102,675 | 0.815 ± 0.002 |

Max-curvature returns Mc = 1.00, but **b is still climbing steeply at Mc = 1.5**.
A stable b above Mc is the whole premise of the estimator; this is not stable.
Formal σ is ±0.001–0.002, i.e. the drift is ~400σ. Two candidate causes, not
mutually exclusive:

1. **The catalog is not complete at Mc = 1.0** over this whole box. Network
   density varies enormously — dense around the Salton Sea Geothermal Field,
   sparse at the basin margins — so one regional Mc cannot be right.
2. **Mixed populations.** The box contains high-b geothermal volumes and
   lower-b tectonic volumes. Superposing FMDs with different b produces
   curvature that mimics incompleteness.

Either way the conclusion is the same and it is the core of this plan:
**a single regional Mc is not defensible here.** At Parkfield we could follow
the paper and assume one; that assumption is what produced the 19.5%
false-positive rate in step 5, and here it would be worse.

## 3. What the literature says to do differently

### 3.1 b-positive (van der Elst 2021, JGR, doi:10.1029/2020JB021027)

Uses only positive differences between successive magnitudes. By conditional
probability the estimator is near-immune to a varying detection threshold —
`b⁺ ≈ b` even for incomplete catalogs.

Run it **alongside** classic Aki at every node. Where the two agree, the node is
complete and the result is trustworthy. Where they diverge, completeness is the
suspect. This converts our worst trap into a routine per-node diagnostic, and it
is the single highest-value addition to the pipeline.

Caveat: Lippiello & Petrillo (2024, JGR doi:10.1029/2023JB027849) show b-positive
loses efficiency when events below threshold are reported, and propose
`b-more-incomplete`. Worth knowing; not worth implementing first.

### 3.2 Significance must be corrected (Marzocchi et al. 2020, GJI doi:10.1093/gji/ggz541)

*How to be fooled searching for significant variations of the b-value* names our
step-5 limitations precisely: the look-elsewhere effect when scanning many
space-time windows, and the fact that **overlapping windows violate the
independence assumption the Utsu test rests on**. At 0.5 km nodes with r = 5 km,
neighbouring samples share nearly all their events.

Their guidance, adopted here:
- Report the **effective number of independent samples**, not the node count.
- Correct for multiple testing on any headline claim.
- Prefer **non-overlapping** samples for a claim; keep the dense grid for the
  picture.
- Pre-specify hypotheses rather than reporting the best of many scans.

### 3.3 Do not decluster (Mizrahi, Nandan & Wiemer 2021, SRL doi:10.1785/0220200231)

Declustering lowers b by **up to 30%**, and they demonstrate on synthetic ETAS
catalogs with known b that the reduction is an artifact of the algorithm rather
than a real difference between mainshocks and aftershocks.

The Salton Trough is swarm-dominated — the Brawley Seismic Zone is driven by
aseismic creep and fluid diffusion, not mainshock-aftershock sequences, so
declustering algorithms are mis-specified here on top of being biased. Default
to **no declustering**, and treat it as a sensitivity test only.

### 3.4 Deferred: Kamer & Hiemer (2015, JGR doi:10.1002/2014JB011510)

Voronoi tessellation + penalized likelihood + BIC selects the spatial partition
from the data — no grid spacing, no radius, and non-overlapping by construction,
which fixes §3.2 as a side effect. This is the principled answer to the radius
problem that knee-detection could only approximate at Parkfield. Substantial
build; revisit once the fixed-radius map exists and if this goes to publication.

## 4. Phases

### Phase A — port to map view
- Wire `MapRegion` through steps 3–5. The estimator (`src/bmap.py`) and plotting
  (`src/bplot.py`) layers are already geometry-agnostic; the work is in the
  node-sampling call and the figure axes.
- `map_region()` alongside `map_section()`, sampling by epicentral distance in
  the projected km frame.
- **Check:** re-run Parkfield through the map-view path and confirm the section
  result is recovered where the two are comparable.

### Phase B — completeness, done properly
- Per-node Mc rather than one regional value.
- ~~Add b-positive to `src/fmd.py`~~ **done** — `fmd.fit_gr_positive`, validated
  in `tests/validate_bpositive.py` (5/5). Under a stepped detection threshold
  the classic estimator is biased −0.218 while b-positive is biased −0.002, at
  ~2.2× the scatter. **Caveat found on real data:** b⁺ is far more stable than
  classic across Mc in the 1.3–1.6 range, but it still drifts at a very
  permissive cut (Mc = 0.5), so it is a cross-check, not a licence to ignore
  completeness entirely. Compare b and b⁺ over a *range* of Mc, not at one value.
- Compute both estimators at every node.
- New diagnostic output: a **b − b⁺ difference map**. Structure in that map is a
  completeness artifact map, and it should be read before the b map is.
- Produce a regional Mc(x, y) map as a first-class product, not an assumption.
- ~~Per-node Mc method: maximum curvature~~ **also implemented: KS-distance**
  (`fmd.mc_ks`, ported from Goebel's `FMD.Mc_KS`/`KS_D_value_PL`, Clauset,
  Shalizi & Newman 2009). Minimises the KS distance between the observed and
  power-law-modelled cumulative distribution among candidates whose Shi & Bolt
  σ(b) is below 0.25, rather than just taking the modal histogram bin.
  Verified bit-for-bit against Goebel's own `FMD` class on synthetic data.
  Wired in as `mc_method="ks"` in `bmap.map_region` (per-node candidates on a
  0.1-mag grid — Mc is never asserted finer than that, and searching every
  unique magnitude per node, Goebel's whole-catalog default, does not scale to
  tens of thousands of nodes).

  **Mixed verdict, not a clean win:**
  | | Parkfield (whole section, paper Mc=1.3 → b=0.92) | Salton (whole box) |
  |---|---|---|
  | maximum curvature | Mc=1.15 → b=1.004 | Mc=1.00 → b=0.644 |
  | KS-distance | Mc=0.90 → b=0.792 | Mc=2.13 → b=1.013 |

  At Parkfield, KS *undershoots* the paper's Mc and misses b by 0.13 — the raw
  Clauset statistic is not very sensitive to a gradual detection rolloff (as
  opposed to the sharp cutoff it was designed for), so a sample that is still
  missing small events can look enough like a power law to win on KS distance.
  Tightening `maxErr_b` from 0.25 to 0.03 does not move it. At Salton, KS goes
  the other way and picks a much higher Mc than max curvature.

  On the map itself (r=8 km, 1 km nodes), KS **does** shrink the completeness
  diagnostic: median(b⁺−b) falls from **+0.170** (max curvature, per-node Mc
  median 1.30) to **+0.106** (KS, per-node Mc median 1.80), at a small coverage
  cost (44.6% → 42.1%). That is the direction "more honest completeness" should
  move it. Both methods are kept (`CONFIG['mc_method'] = "fixed"|"maxcurv"|"ks"`);
  neither is assumed correct on its own — see the b⁺−b map before trusting
  either b map.

### Phase C — the b-value map
- ~~Radius scan re-derived for this region~~ **done, per Mc method** — the
  knee moves with the completeness method, so it had to be re-derived twice,
  not once:

  | mc_method | knee | coverage | contrast at knee |
  |---|---|---|---|
  | maxcurv | **15 km** | 66% | 0.491 |
  | ks | **20 km** | 71% | 0.460 |

  KS's per-node candidates are noisier (a coarser search than max curvature's
  histogram mode), so it needs a larger volume to stabilise — consistent with
  it landing on a higher regional Mc in every check so far (Phase B). Neither
  radius is close to Parkfield's 5 km: contrast itself is nearly flat with
  radius here (maxcurv: 0.556→0.486 across the whole 2–30 km scan, an 11%
  drop, against Parkfield's 67%), so the "knee" is doing less work than it did
  at Parkfield and coverage is the more decisive factor in practice.
  `mc_method` is now baked into step 3/4 output filenames
  (`<name>_<mc_method>_radiusScan.*` etc.) — before this fix, scanning a
  second method silently overwrote the first method's figures under the same
  `CONFIG['prepared']` name.
- b, σ(b), N, and b⁺ maps at the chosen radius. **Done** — `4_bvalueMap.py`'s
  map-view figure is now six panels (b, b⁺, b⁺−b / σ(b), N, per-node Mc), not
  four; previously σ(b) and N were computed and saved but never plotted for
  map view (only for a section).
- Overlay the geothermal wells from `data/wells/` — they stop being decoration
  and become the interpretive frame. **Done**, on all six panels.

**The two final maps**, each at its own knee radius
(`output/parkfield_salton/4_bvalueMap/salton_maxcurv_bvalueMap.png`,
`salton_ks_bvalueMap.png`):

| | maxcurv, r=15 km | KS, r=20 km |
|---|---|---|
| coverage | 65.8% | 71.3% |
| classic b (median) | 0.887 | 1.070 |
| b⁺ (median) | 1.059 | 1.174 |
| median(b⁺−b) | +0.181 | **+0.101** |
| per-node Mc (median) | 1.40 | 2.00 |

KS's smaller completeness gap holds even at its own, larger, more-conservative
radius — the same direction as every other comparison in Phase B. Read
alongside the maps themselves: both show the same qualitative structure (a
low-b band along the main NW-SE trend flanked by higher-b patches), which is
some reassurance the pattern is not an artifact of the completeness method,
only its absolute level and how much of the map that method manages to
resolve. Neither map's reading at the wells themselves is validated — see §5.

### Phase D — significance, honestly
- Effective-independent-sample count and multiple-testing correction. **Done**
  — `bmap.bonferroni`, `bmap.benjamini_hochberg`, `bmap.effective_samples`,
  wired into `5_stationarity.py` as `significance_report()`, run automatically
  after every split (not opt-in).
- Non-overlapping comparison for any headline claim. **Done** — the same call
  re-samples on a grid spaced 2×radius apart (genuinely non-overlapping
  circles/spans, not overlapping-node subsampling) and reports its count as
  the one to quote; the dense grid's raw count is now explicitly labelled a
  picture, not a test.
- Temporal analysis only after Phase B; the completeness check from step 5
  runs first, every time. **Already was** (`completeness_check`, step 5 §2.3
  of the earlier findings), now runs before the new correction too.
- **Scope note:** this only exercises `5_stationarity.py`, which is
  section-geometry only (`load()` refuses a map-view `_prepare.json`). Salton
  was prepared in map view, so **Salton has no temporal-stationarity product
  yet** — Phase D corrected the significance machinery that exists (Parkfield)
  and made it non-optional; a map-view equivalent of step 5 (period-vs-period
  b-difference over the Salton grid, corrected the same way) is not built.
  Add to a future phase if a temporal claim about Salton is wanted.

**The finding, at Parkfield, is not what the paper's own framing suggests.**
Running the existing 1992 split (r = 5 km, the paper's own choice) through the
new correction:

| method | significant | of |
|---|---|---|
| raw (dAIC > 2, no correction) | 58 | 1,938  (2.99%) |
| Bonferroni (α=0.05) | **0** | 1,938 |
| Benjamini-Hochberg FDR (α=0.05) | **0** | 1,938 |
| non-overlapping grid (spacing = 2r, the honest re-run) | **0** | 4 |

The effective-sample estimate explains why: a 111 km section sampled at r=5 km
tiles into only **~11 truly independent spans**. Both the Bonferroni/FDR
correction of the dense grid's own p-values *and* the actual non-overlapping
re-run agree: at the standard α=0.05, **nothing in this section survives
multiple-testing correction** — not even the two dAIC>5 "highly significant"
nodes the raw count contains (Pb ≈ 0.004 and 0.012, both far above the
Bonferroni threshold of 2.6×10⁻⁵ at n=1,938 tests). This does not mean the
paper's underlying physical claim is wrong — a 111 km section was never going
to give many independent samples at 5 km radius, small-N problems are common
in seismology, and the paper's own significance framing (Utsu 1992, no
multiple-testing correction) was standard practice for 2004 — but it does mean
**the specific numeric claim "34 of 2,950 nodes, <1.2%, is stationary" is not
statistically defensible against a corrected null**, and neither is this
pipeline's own uncorrected reproduction of it in `output_findings.md` §5.
That document needs a correction alongside the Salton write-up in Phase E.

### Phase E — write-up
- **Paused (2026-09-03), by decision** — not doing write-ups for now. Deferred:
  extending `output_findings.md` with the Salton Trough results, and correcting
  its stale §5 stationarity claim (see Phase D finding above). Pick back up
  when write-ups are wanted again.
- The other open item, also not started: a Salton temporal-stationarity
  product (map-view equivalent of `5_stationarity.py`) — offered and declined
  for now, stays on the list for later.

## 5. Validation targets

Unlike Parkfield there is no single paper to reproduce, but there are published
results to check against:

- **b is higher inside the Salton Sea Geothermal Field than outside**, with a
  rising trend attributed to poroelastic contraction from production (Zhai,
  Shearer & Fialko-lineage GJI 2023, doi:10.1093/gji/ggac324, "Temporal changes
  of seismicity in Salton Sea Geothermal Field due to distant earthquakes and
  geothermal productions" — mean b ≈ 0.99 inside vs ≈ 0.96 outside, 2008–2013,
  a matched-filter-detected catalog of >70,000 events; Trugman, Shearer, Borsa
  & Fialko 2016, doi:10.1002/2015JB012510). If our map does not show elevated
  b over the SSGF, something is wrong.

  **Checked, and it does not hold in this coarse form.** `data/wells/well_data.csv`
  is 90 wells tightly clustered at 33.015°N, −115.531°W (σ ≈ 2 km) — one field,
  the SSGF. Splitting the whole prepared catalog into rings around that centroid:

  | zone | N | Mc (maxcurv) | b | b⁺ | Mc (KS) | b | b⁺ |
  |---|---|---|---|---|---|---|---|
  | inside 10 km | 7,288 | 1.50 | 0.886 | 1.054 | 2.40 | 0.895 | 0.910 |
  | 10–30 km | 33,920 | 0.50 | 0.405 | 0.926 | 2.70 | 1.026 | 1.132 |
  | beyond 30 km | 270,548 | 1.00 | 0.647 | 0.842 | 3.20 | 1.098 | 1.150 |

  Two things stand out. First, KS is clearly the more trustworthy column here:
  b and b⁺ agree to 0.015 inside 10 km, against a 0.168 gap for max curvature
  — direct evidence for the same conclusion as the map-wide diagnostic above.
  Second, and unresolved: **under KS, b inside the field (0.895) is the
  *lowest* of the three zones**, not the highest — the opposite of the cited
  literature. Before treating that as a real disagreement, note the confounds
  this coarse a check does not control for: a 10 km ring is much larger than
  the SSGF itself and the depth range is 0–20 km with no ring-specific
  Mc-vs-radius check, so "inside" is diluted by non-SSGF crust and possibly
  by mixing shallow (field) and deep (regional) seismicity that the cited
  studies would have separated.

  **Checked at node level too, same conclusion.** Nearest map node to the well
  centroid (r=8 km, 1 km spacing): KS gives b=0.909 (7×7-node window median
  0.896), max curvature gives b=0.883 (window median 0.882) — both sit at or
  slightly *below* the map's regional median (1.012 KS / 0.875 maxcurv), not
  elevated above it. The discrepancy survives moving from a coarse ring to the
  actual fine-grained map product, so it is not a resolution artifact of the
  zone check above.

  **Re-checked against the GJI (2023) study's own window (2008–2013) and the
  SSGF's characteristic depth (shallow, reservoir-associated seismicity is
  reported at median 3.52 ± 0.15 km, 2007–2013 — Scientific Reports 2025,
  doi:10.1038/s41598-025-85744-2) — the mismatch does not go away, and it
  points at the catalog rather than the estimator:**

  | subset | inside 10 km: N, Mc(KS), b(KS), b⁺(KS) | 10–30 km | beyond 30 km |
  |---|---|---|---|
  | whole catalog, all depth | 7,288 / 2.40 / 0.895 / 0.910 | 1.026 / 1.132 | 1.098 / 1.150 |
  | 2008–2013, all depth | 2,019 / 1.40 / 0.860 / 1.026 | 1.104 / 1.212 | 0.913 / 1.000 |
  | shallow (≤6 km), all time | 1,279 / 1.70 / 0.742 / 0.951 | 1.015 / 1.173 | 1.066 / 1.181 |
  | 2008–2013 **and** ≤6 km | 267 / 1.80 / 0.742 / 0.813 | 1.155 / 1.276 | 0.847 / 0.924 |

  In every combination, b inside 10 km is at or below the other two zones —
  never the highest, as the literature reports. That rules out the
  time-period and depth explanations as the *whole* story. What it turned up
  instead: restricted to 2008–2013 within 10 km of the wells, **this catalog's
  events have a median depth of 7.8 km** (IQR 6.9–8.4 km) — more than double
  the 3.52 km the specialized studies report for the same period, and it also
  leaves very few events (267) once both cuts are applied, well below what
  either method needs for a stable per-zone estimate. The most likely
  explanation is a **catalog mismatch, not a b-value mismatch**: the cited
  studies build a matched-filter-detected catalog from dense local
  instrumentation specifically to resolve the SSGF's own shallow, low-magnitude
  induced population; the Hauksson relocated catalog used here is a standard
  regional network product and appears to be seeing a different — deeper,
  more regionally-tectonic — population in that same 10 km ring rather than
  the induced swarms themselves. Building or sourcing a matched-filter SSGF
  catalog is out of scope here; until then, **treat this map's b-value near
  the wells as not validated against the geothermal literature**, and read
  the rest of the map (large-scale structure, the b⁺−b diagnostic) with more
  confidence than the specific SSGF reading.
- Dynamically triggered events cluster near high-b locations (same source).
- The 2012 Brawley swarm and the 1979 Imperial Valley M6.5 are known features
  that should appear sensibly in the maps.

## 6. Risks

| Risk | Mitigation |
|---|---|
| Mc varies so strongly that per-node Mc is itself noisy | b-positive cross-check; report Mc map and b−b⁺ map as products |
| Swarms violate the independence the estimator assumes | no declustering; report sensitivity |
| 312k events over a large box makes the node loop slow | already KD-tree based; profile before optimising |
| Geothermal and tectonic populations mix inside one sampling volume | this is partly the science; small radii and the well overlay address it |
| Real anomalies are confused with network-geometry artifacts | compare against the Mc map before interpreting any anomaly |

## 7. First actions

1. Download the catalog into `data/` and add a `--url` path (with User-Agent) to
   step 1.
2. Run steps 1–2 for the Salton Trough box to get the prepared subset and the
   completeness diagnostics.
3. Implement b-positive in `src/fmd.py` and validate it against a synthetic
   catalog with known b and an imposed time-varying detection threshold.
4. Then Phase A.
