# Candidate methods beyond Schorlemmer et al. (2004)

Not implemented. A survey of papers that advance spatial b-value mapping past
the fixed-radius method this pipeline currently uses (`refs/METHOD.md`), kept
here so a future decision to adopt one of these has the trade-offs written
down rather than re-researched from scratch. No code changes here.

## Tormann, Wiemer & Mignan (2014) — distance-weighted sampling

> Tormann, T., S. Wiemer, and A. Mignan (2014), *Systematic survey of
> high-resolution b-value imaging along Californian faults: Inference on
> asperities*, J. Geophys. Res. Solid Earth, 119, 2029–2054,
> doi:10.1002/2013JB010867.

**Approach.** Proposes "distance exponential weighted" (DEW) sampling as a
replacement for the two prior conventions — constant radius (what this
pipeline uses) and fixed nearest-neighbour count. Instead of a hard cutoff,
every event within a search radius contributes to a node's b-value fit with a
weight that decays exponentially with distance (`exp(-λ·d)`), so nearby
events dominate without a sharp in/out boundary. Two new parameters beyond
Nmin/radius: `Nmax` (a cap on events considered) and `λ` (the decay rate).
Both are tuned once via a synthetic Monte Carlo recovery test — three known
b-value regions (0.5, 1.3, 1.8) embedded in a b=1 background — by minimizing
a score function that measures how well the recovered map matches the known
input. The calibrated method (λ=0.7, effectively Nmax=∞ in the real
application, so the real limiting parameter is Rmax=7.5 km) is then run across
every major Californian fault to build a systematic, high-resolution survey
looking for **asperities**: the paper's central claim is that low b marks
locked, high-differential-stress patches likely to nucleate large ruptures,
and high b marks creeping sections — extending the Parkfield-specific
asperity argument (already the interpretive frame this pipeline's step 4
probability panel uses) to fault systems generally.

**Why not adopted here (yet).** Kamer's comment, below, found the calibration
step circular in a way that matters in practice, not just in principle - see
the two side-by-side Parkfield maps in Kamer's Figure 3 (b range 0.50-1.86 vs
0.80-0.96, same catalog).

## Kamer (2014) — comment on Tormann et al., and why calibration is circular

> Kamer, Y. (2014), Comment on "Systematic survey of high-resolution b-value
> imaging along Californian faults: Inference on asperities" by T. Tormann
> et al., J. Geophys. Res. Solid Earth, doi:10.1002/2014JB011147.
> (Preprint: arXiv:1403.7423)

**The critique.** Tormann et al.'s score function picks λ/Nmax by minimizing
recovery error against a *synthetic* catalog whose true b-value has to be
assumed in advance - but the whole point of the method is to discover the
real spatial b-value distribution, which is exactly what is unknown. Kamer
shows the optimal parameters are not stable across that assumption: a
synthetic input of b=0.5 (what Tormann et al. actually calibrated on, since
it was their dominant synthetic region) wants a small, tightly-localized
sample (λ=0.7, ~150 events), while a synthetic input of b=1 - the single most
well-established regional average in seismology - wants a much larger,
flatter sample (λ=0.01, Rmax≈40 km). Calibrating on the wrong assumed b
undersamples regions whose true b differs from it, and undersampling a
Gutenberg-Richter distribution is known to manufacture spurious high/low-b
artefacts (Felzer 2006, cited therein) rather than reveal real ones -
demonstrated directly by rerunning Tormann et al.'s own Parkfield catalog
under both parameter sets (Kamer's Figure 3): the same data reads as
0.50-1.86 or 0.80-0.96 depending only on which synthetic b was assumed at
calibration time. Kamer's conclusion extends past DEW specifically to
constant-radius and nearest-neighbour sampling too: any method whose
resolution (sample size) is fixed by hand, rather than chosen by the data,
risks confusing an artefact of undersampling for a real spatial anomaly.
Points to likelihood-based alternatives with built-in model-complexity
control (Imoto 1987; Ogata & Katsura 1993) as the way out — which is the
direction Kamer's own follow-up paper below takes.

**Relevance to this project.** This is a direct caution on the method
already in use here (fixed radius, Nmin=50 - itself Schorlemmer's own
judgement call, not re-derived, though `scripts/2a_nMinStability.py` now at
least checks it empirically per catalog). It does not invalidate fixed-radius
sampling, but it is the reason any *contrast* number this pipeline reports
(step 3's radius scan, `bmap.heterogeneity`) should be read as an upper bound
on real spatial structure, not a clean measurement of it - some of it is
provably attributable to sampling choice rather than the ground, and the
scan can only show how contrast falls off with radius, not what fraction of
the value at any one radius is real.

## Kamer & Hiemer (2015) — data-driven partitioning, no radius or Nmin at all

> Kamer, Y., and S. Hiemer (2015), Data-driven spatial b value estimation
> with applications to California seismicity: To b or not to b, J. Geophys.
> Res. Solid Earth, 120, 5191–5214, doi:10.1002/2014JB011510.

(Already named as a deferred candidate in `PLAN.md`; this is the first full
summary of what it actually does.)

**Approach.** Removes the radius/Nmin/λ choice entirely. The study area is
partitioned into Voronoi cells around a set of randomly placed nodes; b is
fit once per cell (plain Aki MLE, no distance weighting), and the whole
map's joint likelihood is scored. The number of nodes (the map's
*complexity* - more nodes means more free parameters, each fit on fewer
events) is itself swept and selected by a penalized-likelihood criterion, so
the data decide how much spatial resolution the catalog can actually support
rather than an analyst picking a radius by eye or by a synthetic test. Mc is
handled per partition with the Clauset et al. (2009) KS-distance method -
the same estimator already implemented here as `fmd.mc_ks` / `fmd.ks_distance_curve`.

**Main finding.** Applied to California statewide, the data-driven model
selects far less spatial complexity than prior fixed-radius/DEW studies
implied: b varies only across a narrow 0.94±0.04 to 1.15±0.06 range at most
locations, markedly tighter than Tormann et al.'s reported 0.5-1.86. The
paper's reading is that most of the dramatic spatial b-value contrast in the
earlier literature is a parameter-choice artefact (exactly what Kamer's 2014
comment predicted), not evidence that California's crust is that
heterogeneous.

**Why this is the strongest deferred candidate.** It directly answers this
project's own step 3 problem (Kamer's earlier critique of the fixed-radius
convention this pipeline also uses) with an actual alternative rather than
just a caution, and it already uses this project's own KS-Mc estimator as a
component. The cost is implementation complexity: Voronoi tessellation
search plus a penalized-likelihood model-selection loop is a materially
larger piece of code than `bmap.map_region`'s fixed-radius grid, and would
need its own validation pass against Parkfield before being trusted on
Salton or El Salvador.

## Gulia & Wiemer (2019) — b-value as a real-time forecasting signal

> Gulia, L., and S. Wiemer (2019), Real-time discrimination of earthquake
> foreshocks and aftershocks, Nature, 574, 193–199,
> doi:10.1038/s41586-019-1606-4.

**Approach.** Not a new spatial-mapping method - it reuses the same
machinery (regional Mc, Aki b, a background/reference b-value) but turns it
into a real-time monitor rather than a static map. As a sequence unfolds
after a moderate earthquake, b is recomputed in a shrinking time window and
compared to the pre-sequence background: a b-value drop of more than ~10%
below background is read as a foreshock warning (more large events
relatively likely, i.e. the sequence may culminate in something bigger), a
rise of more than ~10% as reassurance the sequence is a decaying aftershock
series, and anything in between as ambiguous - a three-state "traffic light"
built directly on the b-value-as-stress-proxy idea this whole pipeline
already uses (low b -> closer to failure).

**Relevance to this project.** Orthogonal to the other three papers here -
those are about *where* to trust a b-value estimate spatially; this is about
using the *same* estimate as a live signal through time. It is the natural
next step if a temporal-monitoring product is ever wanted for either study
area (this pipeline's step 5 already does periodwise b comparison with
proper significance correction, per Phase D in `PLAN.md`, but as a
retrospective stationarity check, not a live monitor), and it is a concrete
reason to keep the b-positive/completeness-correction machinery this
pipeline already has, rather than a simple Aki fit: a live monitor is
exactly the setting where a moving detection threshold does the most damage
to a naive b-value comparison.

## Not covered here

This is not an exhaustive literature survey - it is the four papers directly
upstream of a "should this pipeline change its sampling method" decision,
found by following citations from the paper named in the request (Tormann et
al. 2014) forward and backward. A broader sweep (Bayesian Mc estimation,
ETAS-based declustering-aware b-value work, ML-based magnitude estimation)
would be a separate pass if one of these four turns into an actual adoption
decision.
