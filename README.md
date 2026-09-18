[![DOI](https://zenodo.org/badge/811449976.svg)](https://zenodo.org/doi/10.5281/zenodo.11508510)

# pybvalue

Spatial mapping of earthquake b-values, following the method of
**Schorlemmer, Wiemer & Wyss (2004)**, *Earthquake statistics at Parkfield:
1. Stationarity of b values*, JGR 109, B12307, doi:10.1029/2004JB003234.

The pipeline is **area-agnostic**: region, time span, depth cut, magnitude
binning, sampling radius and `Nmin` are all inputs. Parkfield is the
reproduction target used to validate the code; the Salton Trough / Imperial
Valley is the science target.

The extracted method — every formula, parameter and reference value needed to
reproduce the paper — is in **[`refs/METHOD.md`](refs/METHOD.md)**. Read that
before changing anything in `src/`.

## Layout

```
data/                  input catalogs (large files are gitignored)
docs/                  project description and report
refs/                  METHOD.md - the method, extracted from the source paper
scripts/               numbered workflow steps
src/                   importable modules
output/<study>/<step>/ everything produced, mirroring the script numbering
```

Scripts are numbered by workflow step, and each writes into the matching
`output/<study>/<step>/` folder, where `<study>` is `CONFIG['study']` -
a name grouping a set of related regions under one output tree (default
`parkfield_salton`, covering both regions below; a new study, e.g.
`elsalvador`, keeps its outputs separate - set `CONFIG['study']` the same
way in every step of that run). `output/` contents are gitignored; the
folder skeleton is not.

## Workflow

| Step | Script | Produces |
|---|---|---|
| 1 | `1_downloadCatalog.py` | a canonical catalog + provenance sidecar |
| 2 | `2_prepareCatalog.py` | study-volume subset + completeness/FMD diagnostics |
| 3 | `3_radiusScan.py` | b mapped over r = 2…20 km; pick the radius |
| 4 | `4_bvalueMap.py` | b-value map at the chosen radius, and the hazard it implies |
| 5 | `5_stationarity.py` | Δb between periods + Utsu significance |

Steps 3–5 are Schorlemmer's three-step workflow. (The regional Mc and b-value
baseline that would have been a separate step is already produced by step 2.)
Step 3 is not optional boilerplate: the paper is explicit that the optimal
radius depends on the local seismotectonic fabric and data availability, so it
must be re-derived for every new region rather than copied from Parkfield.

**No command-line flags.** Every script has a `CONFIG = {...}` dict near the
top instead of an argparse parser — open the file, edit the values, run
`python scripts/N_something.py`. Each script's default `CONFIG` is a working
example (Parkfield, reproducing the numbers below), so running any step with
no edits at all reproduces the reference case; switching region or parameter
is a one-line edit, not a new set of flags to remember.

## Step 1 — acquire a catalog

One step, two paths in, one format out. Download fresh from an FDSN event
service, or normalise a catalog you already have. Both write the same canonical
CSV to `output/<study>/1_catalog/` alongside a `.json` sidecar recording
exactly where it came from, and a `<name>_catalog.png`/`.pdf` reconnaissance
figure (epicenters on a map, magnitude through time) - there's no geometry yet
at this stage, so it's deliberately just the two things you need to choose one.
The map has a real basemap under it - shaded relief, political boundaries,
place names (Esri World Topo, via `contextily`; `pip install contextily` if
it's not already there) - and degrades gracefully to a plain background if
that package, or the network, isn't available.

```python
# in scripts/1_downloadCatalog.py, CONFIG:

CONFIG["region"] = "parkfield"          # fresh download, named region (the default)

# any other area:
CONFIG["region"] = None
CONFIG["datacenter"] = "NCEDC"
CONFIG["bbox"] = (35.6, 36.3, -120.8, -120.0)
CONFIG["start"] = "1981-01-01"
CONFIG["name"] = "myregion"

# a catalog you already have (format auto-detected):
CONFIG["region"] = None
CONFIG["file"] = "data/1977_20240611.txt"
CONFIG["name"] = "socal"
```
then `python scripts/1_downloadCatalog.py`. Named regions are in the
`REGIONS` dict at the top of the file, next to `CONFIG`.

Canonical columns: `time, latitude, longitude, depth, magnitude, magtype, evid`
— depth in km positive down, time UTC. Readers exist for FDSN text, SCEDC
`date_mag_loc`, GrowClust `.gc`, the Hauksson–Yang–Shearer relocated `.mat`, and
generic CSV.

Two things the loader handles that are easy to get wrong:

- **FDSN servers truncate at 10,000 rows silently** rather than erroring. Any
  window returning exactly that many is split and re-requested until every
  piece is under the cap.
- **Network catalogs round origin seconds to `60.00`** (really 59.995), which
  pandas cannot parse. Those events are recovered rather than dropped — 17 of
  them in the bundled SCEDC catalog. Any event that *is* dropped is reported
  loudly, because losing events silently biases every b-value downstream.

## Step 2 — prepare the study volume

Projects the catalog onto a cross section (or a map region), applies the time,
depth and geometry cuts, and writes the diagnostics that justify the
completeness choice — the equivalents of Figures 1–3 of the paper.

```python
CONFIG["catalog"] = "parkfield"
CONFIG["preset"] = "parkfield"
```
then `python scripts/2_prepareCatalog.py` (this is the file's default `CONFIG`,
so running it unedited does the same thing).

Produces `<name>_section.csv` plus map, section, completeness and FMD figures,
and a `_prepare.json` recording **every cut and how many events it removed** —
so a surprising event count can always be traced to the cut that caused it.

Geometry lives in `src/geometry.py` as `CrossSection` and `MapRegion`, sharing a
`project` / `select` interface. All distances are kilometres on a locally
projected plane (azimuthal equidistant, centred on the region), never degrees.

Parkfield reproduction, against Schorlemmer et al.:

| | this pipeline | paper |
|---|---|---|
| section length | 111.0 km | 110 km |
| N(M ≥ 1.3, D ≤ 16 km) | 3,618 | 3,780 |
| regional b | 0.903 ± 0.014 | 0.92 |

The residual is catalog vintage: NCEDC today is not NCSN as it stood in 2004.

## Step 3 — choose the sampling radius

```python
CONFIG["prepared"] = "parkfield"
CONFIG["features"] = True
```
then `python scripts/3_radiusScan.py`.

Maps b at every node for each radius and reports the two quantities the choice
trades off — **coverage** (fraction of nodes with ≥ Nmin events) and
**contrast** (p95 − p5 of b), plus the correlation of each map with the finest
one. Parkfield:

| r (km) | coverage | contrast | corr vs finest |
|---|---|---|---|
| 2 | 15% | 0.814 | 1.000 |
| 4 | 47% | 0.636 | 0.707 |
| 5 | 59% | 0.567 | 0.674 |
| 6 | 69% | 0.521 | 0.615 |
| 10 | 93% | 0.419 | 0.496 |
| 20 | 100% | 0.272 | 0.294 |

The script suggests a radius from the knee of the contrast curve, which lands
at 6 km against the paper's 5. **That gap is honest, not a bug**: the
correlation between maps decays smoothly with no sharp break, so no automatic
rule recovers the paper's choice — Schorlemmer et al. picked 5 km by reading
the panels, and say so. The suggestion is where to start looking.

`src/bmap.py` holds the mapping core (`map_section`, `radius_scan`,
`compare_maps`); `src/bplot.py` holds the shared colour scale.

**Map view** (a region rather than a single fault section, e.g. the Salton
Trough) samples epicentral cylinders instead of a section plane, and adds a
completeness-method choice per node, `CONFIG["mc_method"]`:

```python
CONFIG["prepared"] = "salton"
CONFIG["mc_method"] = "maxcurv"   # or "ks"
```

| `CONFIG["mc_method"]` | what it does |
|---|---|
| `"fixed"` (default) | one Mc for the whole map, from step 2 or `CONFIG["mc"]` |
| `"maxcurv"` | per-node Mc, the modal bin of that node's histogram |
| `"ks"` | per-node Mc, minimum KS distance to a power law (Clauset et al. 2009) among candidates with Shi & Bolt σ(b) < 0.25 |

Map-view output filenames carry the method (`<name>_<mc_method>_radiusScan.png`
etc.) — a section always uses one fixed Mc, so its filenames are untouched.
Neither per-node method won outright at the Salton Trough: max curvature keeps
a Parkfield-consistent Mc≈1.15, but KS's diagnostic (b⁺ − b, below) is smaller
there. Both are kept rather than picked, and the two knees landed at different
radii (15 km for maxcurv, 20 km for KS) — see `PLAN.md`, Phase B/C, for the
numbers behind that.

## Step 4 — the b-value map, and what it implies

```python
CONFIG["prepared"] = "parkfield"
CONFIG["radius"] = 5
CONFIG["features"] = True
```
then `python scripts/4_bvalueMap.py`.

Maps b, σ(b) and sample size on the section, then converts a and b into the
annual probability of one or more M ≥ M₀ (paper equations 3–4), computed twice:
once holding b at the regional average, once letting it vary in space.

That second figure is the paper's argument, and it reproduces:

| | peak annual P(M ≥ 6) | where |
|---|---|---|
| constant b = 0.90 | 0.0026 /yr | 55.5 km along, 4 km deep — the **creeping** section |
| a and b both varying | 0.0301 /yr | 71.0 km along, 13 km deep — the **Middle Mountain asperity** |

Schorlemmer et al. place the asperity at 70 km from P1. Assuming a single
regional b puts the highest M ≥ 6 hazard in the creeping segment, which
contradicts the observations; letting b vary moves the peak onto the asperity
where the 1966 Parkfield event nucleated. Both M ≥ 4.5 events in the catalog
(67.1 and 72.0 km along) sit inside that low-b patch.

The map itself matches Figure 5: a b ≈ 0.45 patch at the asperity, a low-b zone
running ~25 km southeast into the locked section, and b up to 1.88 in the
creeping section. σ(b) averages 0.061 with 90% below 0.088 — lower than the
paper's 0.116/0.122, as expected, since those are for half-length periods.

**Note on `constant_b_avalue`:** the a value absorbs b through its
normalisation (a = log₁₀N + b·Mc), so the constant-b panel recomputes a rather
than reusing the values fitted with a varying b. Skipping that step silently
mixes the two models and flatters the comparison.

**Map view** produces a six-panel figure instead of the section's stacked b /
σ(b) / N: classic b, b-positive and their difference on top; σ(b), N and
per-node Mc underneath. Geothermal wells (`data/wells/`, if present) are
overlaid on every panel.

```python
CONFIG["prepared"] = "salton"
CONFIG["radius"] = 20
CONFIG["mc_method"] = "ks"
```

Read panel (c), b⁺ − b, **before** panel (a): a large value there means the
node is incomplete, not that it is tectonically distinctive, and panel (a)'s
regional structure should be trusted in proportion to how flat panel (c) is
nearby.

## Step 5 — is the pattern stationary?

```python
CONFIG["prepared"] = "parkfield"
CONFIG["split"] = "1992-01-01"    # or CONFIG["scan"] = True to sweep every year
CONFIG["mc"] = 1.6
```
then `python scripts/5_stationarity.py`.

Splits the catalog in two, maps b in each period on the same grid, differences
them, and runs the Utsu test at every node. Six panels (b₁, b₂, Δb, σ₁, σ₂,
log₁₀ Pb) plus an optional sweep of the division date.

**The Utsu test cannot tell a local change from a network-wide one.** If
completeness improves between the periods, every node's b moves together and
the test flags them all. At the paper's Mc = 1.3 on the modern NCEDC catalog
that is exactly what happens:

| Mc | catalog-wide b shift | median node Δb | nodes significant |
|---|---|---|---|
| 1.3 | **+0.068** | **+0.069** | 19.5% |
| 1.4 | −0.005 | −0.012 | 11.6% |
| 1.6 | −0.027 | +0.008 | 3.0% |

At Mc = 1.3 the median node change *is* the global shift — almost none of it is
spatial. The 1981–1992 period is not reliably complete at 1.3 in today's
catalog, so it is missing small events and its b comes out too low. The script
detects this and warns. Diagnose with the **shift**, which is power-independent,
not with the significance count — raising Mc also discards events and costs
power.

At a completeness-safe Mc = 1.6 the paper's conclusion reproduces:

| split | nodes compared | significant | |
|---|---|---|---|
| 1992 | 1,938 | 58 | 2.99% |
| 1996 | 1,442 | 1 | 0.07% |
| paper (1992) | 2,950 | 34 | <1.2% |

and the division scan puts most splits at 0–1.2%, i.e. **stationary**, with a
sharp isolated peak at 1991–1992 that collapses by 1993. That brackets a change
to ≈1992–93 — the paper reports one patch whose b rose significantly, "initiating
around 1993", correlated with a documented episode of creep at depth.

### On the colour scale

b-value is rendered as **diverging**, not sequential: what carries meaning is
which side of the regional average a node sits on. The ramp is warm ↔ cool with
a **neutral grey midpoint pinned to the regional b**, so "average" reads as
absent rather than as a value of its own, and every panel of a scan shares one
scale so the panels can actually be compared. Unresolved nodes get a flat
neutral and never blend into the low end of the ramp.

The published figures use a rainbow ramp. That is not reproduced — a multi-hue
ramp for a magnitude invents boundaries the data does not have and is
unreadable to colour-vision-deficient readers.

## Environment

Developed against the `mypy` conda environment: Python 3.12, pygmt 0.15 with
GMT 6.5, numpy, pandas, scipy, matplotlib, obspy, geopandas, pyproj.

The b-value / Mc estimator itself comes from
[tgoebel/magnitude-distribution](https://github.com/tgoebel/magnitude-distribution)
(`src/FMD_GR.py`: maximum-curvature and KS completeness, Aki MLE, Shi & Bolt
uncertainty, and the Utsu test). Clone it into the project root — it is
gitignored:

```bash
git clone https://github.com/tgoebel/magnitude-distribution
```

## References

Aki, K., 1965, Maximum likelihood estimate of b in the formula log N = a − bM
and its confidence limits: Bull. Earthquake Res. Inst., Tokyo Univ., v. 43,
p. 237–239.

Clauset, A., Shalizi, C.R., and Newman, M.E.J., 2009, Power-law distributions in
empirical data: SIAM Review, v. 51, no. 4, p. 661–703.

Hauksson, E., Yang, W., and Shearer, P.M., 2012, Waveform relocated earthquake
catalog for southern California (1981 to June 2011): BSSA, v. 102, no. 5,
p. 2239–2244.

Schorlemmer, D., Wiemer, S., and Wyss, M., 2004, Earthquake statistics at
Parkfield: 1. Stationarity of b values: JGR Solid Earth, v. 109, B12307,
doi:10.1029/2004JB003234.

Schorlemmer, D., Wiemer, S., Wyss, M., and Jackson, D.D., 2004, Earthquake
statistics at Parkfield: 2. Probabilistic forecasting and testing: JGR Solid
Earth, v. 109, B12308, doi:10.1029/2004JB003235.

Shi, Y., and Bolt, B.A., 1982, The standard error of the magnitude-frequency
b value: BSSA, v. 72, no. 5, p. 1677–1687.

Utsu, T., 1992, Certain aspects of magnitude frequency relation of earthquakes:
Zisin, v. 45, p. 65–74.

Wiemer, S., and Wyss, M., 2002, Mapping spatial variability of the
frequency-magnitude distribution of earthquakes: Advances in Geophysics, v. 45,
p. 259–302.
