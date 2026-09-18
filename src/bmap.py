"""
Spatial b-value mapping.

At every node of a regular grid, take the events inside a sampling volume of
radius r, require at least Nmin of them above Mc, and fit the Gutenberg-Richter
b-value. Nodes with too few events stay undefined - the masked area is part of
the result, not a rendering detail, because coverage is exactly what the radius
trades against resolution.

Following Schorlemmer et al. (2004): fixed radius (not fixed-N nearest
neighbours), and for a cross section the radius acts in the (along, depth)
plane while the section width acts across it. See refs/METHOD.md.
"""

import numpy as np
from scipy.spatial import cKDTree

from . import fmd


def map_section(events, section, radius_km, mc, nmin=50, binsize=None,
                spacing_km=0.5, depth_range=(0.0, 16.0), verbose=False):
    """
    Map b on the plane of a cross section.

    Events must already be restricted to the section volume (step 2) - the
    section width has done its work by then, so sampling here is 2-D in
    (along, depth).

    Args:
        events (pd.DataFrame): prepared catalog with `along`, `depth`, `magnitude`
        section (geometry.CrossSection): the section being mapped
        radius_km (float): sampling radius in the section plane
        mc (float): magnitude of completeness
        nmin (int): minimum events above Mc for a node to be computed
        binsize (float): magnitude bin width; detected from the data when None
        spacing_km (float): node spacing
        depth_range (tuple): (min, max) depth for the node grid, km
        verbose (bool): print a one-line summary
    Returns:
        dict: along, depth (1-D node axes); b, sigma, a, n (2-D, NaN where
              unresolved); plus radius_km, mc, nmin, binsize, coverage, n_nodes
    """
    if binsize is None:
        binsize = fmd.detect_binsize(events["magnitude"])
    correction = binsize / 2.0

    complete = events[events["magnitude"] >= mc]
    along_axis, depth_axis, grid_along, grid_depth = section.node_grid(
        spacing_km=spacing_km, depth_range=depth_range)
    shape = grid_along.shape

    b = np.full(shape, np.nan)
    sigma = np.full(shape, np.nan)
    a_value = np.full(shape, np.nan)
    counts = np.zeros(shape, dtype=int)

    if len(complete) >= nmin:
        tree = cKDTree(np.column_stack([complete["along"].values,
                                        complete["depth"].values]))
        magnitudes = complete["magnitude"].values
        nodes = np.column_stack([grid_along.ravel(), grid_depth.ravel()])
        neighbours = tree.query_ball_point(nodes, r=radius_km)

        flat_b = b.ravel()
        flat_sigma = sigma.ravel()
        flat_a = a_value.ravel()
        flat_n = counts.ravel()
        log10e = np.log10(np.e)

        for i, index in enumerate(neighbours):
            n = len(index)
            flat_n[i] = n
            if n < nmin:
                continue
            sample = magnitudes[index]
            mean_mag = sample.mean()
            denominator = mean_mag - (mc - correction)
            if denominator <= 0:
                continue
            value = log10e / denominator
            flat_b[i] = value
            flat_sigma[i] = (2.30 * value**2
                             * np.sqrt(((sample - mean_mag) ** 2).sum() / (n * (n - 1))))
            flat_a[i] = np.log10(n) + value * mc

        b, sigma, a_value = (flat_b.reshape(shape), flat_sigma.reshape(shape),
                             flat_a.reshape(shape))
        counts = flat_n.reshape(shape)

    resolved = np.isfinite(b)
    result = {
        "along": along_axis, "depth": depth_axis,
        "grid_along": grid_along, "grid_depth": grid_depth,
        "b": b, "sigma": sigma, "a": a_value, "n": counts,
        "radius_km": float(radius_km), "mc": float(mc), "nmin": int(nmin),
        "binsize": float(binsize), "spacing_km": float(spacing_km),
        "n_nodes": int(b.size), "n_resolved": int(resolved.sum()),
        "coverage": float(resolved.mean()),
    }
    result.update(heterogeneity(b))
    if verbose:
        print(f"  r={radius_km:4.1f} km  coverage={result['coverage']*100:5.1f}%  "
              f"b={result['b_min']:.2f}..{result['b_max']:.2f}  "
              f"contrast={result['b_contrast']:.3f}")
    return result


def map_region(events, region, radius_km, mc=None, nmin=50, binsize=None,
               spacing_km=1.0, mc_method="fixed", positive=True,
               ks_step=0.1, ks_max_err=0.25, verbose=False):
    """
    Map b in map view: vertical cylinders of radius r about each grid node.

    Differs from `map_section` in three ways that matter for a basin rather
    than a single fault strand:

      * sampling is epicentral, in a locally projected km frame;
      * Mc can be estimated **per node** instead of assumed constant, because
        network density varies across a region like the Salton Trough far more
        than it does along one fault;
      * b-positive is computed alongside the classic estimator, so the two can
        be differenced into a map of where completeness is suspect.

    Args:
        events (pd.DataFrame): needs `latitude`, `longitude`, `magnitude`,
            `time`. Sorted by time internally - b-positive requires it.
        region (geometry.MapRegion): the area being mapped
        radius_km (float): sampling radius
        mc (float): completeness when `mc_method` is 'fixed'; also the floor
            below which a per-node estimate is not allowed to fall
        nmin (int): minimum events above Mc for a node to be computed
        binsize (float): magnitude bin width; detected when None
        spacing_km (float): node spacing. 0.5 km over a whole basin is a very
            large grid; 1 km is a saner default here than at Parkfield
        mc_method (str): 'fixed', 'maxcurv' (per-node maximum curvature) or
            'ks' (per-node minimum-KS-distance completeness, Clauset et al.
            2009 / Goebel's FMD.Mc_KS - looks at the whole magnitude range
            above each candidate rather than just the modal bin, so it is
            less easily fooled by a catalog that is locally complete but
            still rolling off just above its mode)
        positive (bool): also compute b-positive at every node
        ks_step (float): candidate spacing for mc_method='ks'. Mc is never
            asserted finer than this - searching every unique magnitude per
            node (Goebel's whole-catalog default) is both wasted resolution
            and, over tens of thousands of nodes, far too slow
        ks_max_err (float): mc_method='ks' rejects candidates with Shi & Bolt
            sigma(b) at or above this, to keep the search out of the tail
        verbose (bool): print a one-line summary
    Returns:
        dict: x, y (1-D node axes, km); grid_x, grid_y; b, sigma, a, n, mc_node
              (2-D, NaN where unresolved); b_positive, sigma_positive,
              n_positive when `positive`; plus the usual metadata and coverage
    """
    frame = events.sort_values("time")
    magnitudes = frame["magnitude"].to_numpy(dtype=float)
    if binsize is None:
        binsize = fmd.detect_binsize(magnitudes)

    projected = region.project(frame)
    coords = np.column_stack([projected["x"].to_numpy(), projected["y"].to_numpy()])
    x_axis, y_axis, grid_x, grid_y = region.node_grid(spacing_km=spacing_km)
    shape = grid_x.shape

    tree = cKDTree(coords)
    nodes = np.column_stack([grid_x.ravel(), grid_y.ravel()])
    neighbours = tree.query_ball_point(nodes, r=radius_km, workers=-1)

    fields = {name: np.full(shape, np.nan).ravel()
              for name in ("b", "sigma", "a", "mc_node", "b_positive",
                           "sigma_positive")}
    counts = np.zeros(shape, dtype=int).ravel()
    counts_positive = np.zeros(shape, dtype=int).ravel()
    floor = mc if mc is not None else -np.inf

    ks_candidates = fmd.ks_candidate_grid(magnitudes, ks_step) if mc_method == "ks" else None

    for i, index in enumerate(neighbours):
        if len(index) < nmin:
            continue
        # Events were sorted by time before the tree was built, so sorting the
        # neighbour indices restores time order - which b-positive requires and
        # query_ball_point does not preserve.
        index = np.sort(np.asarray(index))
        sample = magnitudes[index]

        if mc_method == "maxcurv":
            node_mc = fmd.mc_maxcurv(sample, binsize or 0.1)
            if not np.isfinite(node_mc):
                continue
            node_mc = max(node_mc, floor) if np.isfinite(floor) else node_mc
        elif mc_method == "ks":
            local_candidates = ks_candidates[ks_candidates >= sample.min()]
            node_mc = fmd.mc_ks(sample, candidates=local_candidates,
                                binsize=binsize, maxErr_b=ks_max_err)
            if not np.isfinite(node_mc):
                continue
            node_mc = max(node_mc, floor) if np.isfinite(floor) else node_mc
        else:
            node_mc = mc
        fields["mc_node"][i] = node_mc

        classic = fmd.fit_gr(sample, node_mc, binsize=binsize, min_events=nmin)
        counts[i] = classic["n"]
        if classic["n"] < nmin or not np.isfinite(classic["b"]):
            continue
        fields["b"][i] = classic["b"]
        fields["sigma"][i] = classic["sigma"]
        fields["a"][i] = classic["a"]

        if positive:
            plus = fmd.fit_gr_positive(sample, binsize=binsize, mc=node_mc)
            if np.isfinite(plus["b"]) and plus["n"] >= max(2, nmin // 4):
                fields["b_positive"][i] = plus["b"]
                fields["sigma_positive"][i] = plus["sigma"]
                counts_positive[i] = plus["n"]

    result = {"x": x_axis, "y": y_axis, "grid_x": grid_x, "grid_y": grid_y,
              "n": counts.reshape(shape), "n_positive": counts_positive.reshape(shape),
              "radius_km": float(radius_km), "mc": mc, "mc_method": mc_method,
              "nmin": int(nmin), "binsize": float(binsize),
              "spacing_km": float(spacing_km), "n_nodes": int(counts.size)}
    result.update({name: field.reshape(shape) for name, field in fields.items()})

    resolved = np.isfinite(result["b"])
    result["n_resolved"] = int(resolved.sum())
    result["coverage"] = float(resolved.mean())
    result.update(heterogeneity(result["b"]))

    both = resolved & np.isfinite(result["b_positive"])
    difference = np.where(both, result["b_positive"] - result["b"], np.nan)
    result["b_difference"] = difference
    result["median_b_difference"] = (float(np.nanmedian(difference))
                                     if both.any() else np.nan)

    if verbose:
        line = (f"  r={radius_km:4.1f} km  coverage={result['coverage']*100:5.1f}%  "
               f"b={result['b_min']:.2f}..{result['b_max']:.2f}")
        # b-positive is only computed when positive=True (the radius scan
        # skips it - the choice of radius is made on coverage/contrast alone,
        # and computing it at every radius would be needless cost); printing
        # "median(b+ - b)=+nan" regardless, as if it had been tried and come
        # up empty, would misreport a field that was simply never asked for.
        if positive:
            line += f"  median(b+ - b)={result['median_b_difference']:+.3f}"
        print(line)
    return result


def heterogeneity(b):
    """
    How much b-value structure a map retains.

    The radius scan is a trade-off: larger radii cover more nodes but smooth
    anomalies away. Coverage measures the first, these measure the second, so
    the two together are what the choice of radius is actually made on.

    Args:
        b (np.ndarray): mapped b values, NaN where unresolved
    Returns:
        dict: b_min, b_max, b_median, b_std, b_contrast (p95-p5), b_iqr
    """
    values = b[np.isfinite(b)]
    if values.size == 0:
        return {k: np.nan for k in
                ("b_min", "b_max", "b_median", "b_std", "b_contrast", "b_iqr")}
    p5, p25, p75, p95 = np.percentile(values, [5, 25, 75, 95])
    return {"b_min": float(values.min()), "b_max": float(values.max()),
            "b_median": float(np.median(values)), "b_std": float(values.std()),
            "b_contrast": float(p95 - p5), "b_iqr": float(p75 - p25)}


def radius_scan(events, section, radii, mc, nmin=50, binsize=None,
                spacing_km=0.5, depth_range=(0.0, 16.0), verbose=True):
    """
    Map b over a range of sampling radii - step 1 of the Schorlemmer workflow.

    The paper's rule: take the *largest* radius that still resolves the
    heterogeneity. Smaller radii only cost coverage without revealing new
    detail; larger ones blur real contrasts away. The optimum is local and must
    be re-derived for every region.

    Args:
        events, section, mc, nmin, binsize, spacing_km, depth_range: see map_section
        radii (iterable): sampling radii, km
        verbose (bool): print a line per radius
    Returns:
        list: one map_section result per radius, in the order given
    """
    if verbose:
        print(f"radius scan over {len(list(radii))} radii "
              f"({len(events):,} events, Mc={mc}, Nmin={nmin})")
    return [map_section(events, section, radius, mc, nmin=nmin, binsize=binsize,
                        spacing_km=spacing_km, depth_range=depth_range,
                        verbose=verbose)
            for radius in radii]


def recurrence(a, b, m0, duration_years):
    """
    Probabilistic recurrence time and annual probability (paper equations 3-4).

        Tr = dT / 10^(a - b*M0)          local recurrence time, years
        lambda = 1 / Tr                  annual rate
        Pr = 1 - exp(-lambda)            P(one or more M >= M0 in a year)

    Args:
        a, b (np.ndarray | float): Gutenberg-Richter parameters, per node
        m0 (float): target magnitude
        duration_years (float): length of the period the a value came from
    Returns:
        tuple: (recurrence time in years, annual probability)
    """
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        expected = 10.0 ** (np.asarray(a, dtype=float)
                            - np.asarray(b, dtype=float) * m0)
        rate = expected / float(duration_years)
        time = np.where(rate > 0, 1.0 / rate, np.inf)
        probability = 1.0 - np.exp(-rate)
    return time, probability


def constant_b_avalue(result, b_constant):
    """
    Recompute a at every node as if b were the regional average.

    The a value absorbs b through its normalisation (a = log10(N) + b*Mc), so
    holding b fixed while reusing the a values fitted with a spatially varying
    b would mix the two models. The paper's comparison - "spatially varying a,
    constant b" against "both varying" - needs this recomputation to be honest.

    Args:
        result (dict): output of map_section
        b_constant (float): the regional b value
    Returns:
        np.ndarray: a values on the same grid, NaN where unresolved
    """
    counts = result["n"].astype(float)
    valid = np.isfinite(result["b"]) & (counts > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(valid, np.log10(counts) + b_constant * result["mc"], np.nan)


def compare_maps(first, second):
    """
    Difference two b-value maps and test each node with the Utsu test.

    Used for the stationarity step: two time periods mapped on the same grid,
    then asked whether each node's change is significant.

    Args:
        first, second (dict): results from map_section, same grid
    Returns:
        dict: delta_b, p_b, d_aic, log_p (2-D, NaN where either map is
              unresolved), plus counts of significant and highly significant
              nodes
    """
    if first["b"].shape != second["b"].shape:
        raise ValueError("maps are on different grids")

    both = (np.isfinite(first["b"]) & np.isfinite(second["b"])
            & (first["n"] >= 2) & (second["n"] >= 2)
            & (first["b"] > 0) & (second["b"] > 0))
    delta = np.where(both, second["b"] - first["b"], np.nan)

    # Vectorised form of fmd.utsu_test - the division scan calls this on tens of
    # thousands of nodes, so the per-node Python loop was the bottleneck.
    n1 = np.where(both, first["n"], np.nan).astype(float)
    n2 = np.where(both, second["n"], np.nan).astype(float)
    b1 = np.where(both, first["b"], np.nan)
    b2 = np.where(both, second["b"], np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        total = n1 + n2
        d_aic = (-2 * total * np.log(total)
                 + 2 * n1 * np.log(n1 + n2 * b1 / b2)
                 + 2 * n2 * np.log(n1 * b2 / b1 + n2)
                 - 2)
        p_b = np.exp(-d_aic / 2 - 2)
        log_p = np.log10(p_b)

    return {"delta_b": delta, "p_b": p_b, "d_aic": d_aic, "log_p": log_p,
            "n_compared": int(both.sum()),
            "n_significant": int(np.nansum(d_aic > 2)),
            "n_highly_significant": int(np.nansum(d_aic > 5))}


# ---------------------------------------------------------------------------
# significance, corrected for multiple testing (Marzocchi, Zechar & Jordan
# 2020 do not name this paper by title here; see PLAN.md Phase D /
# doi:10.1093/gji/ggz541, "How to be fooled searching for significant
# variations of the b-value")
# ---------------------------------------------------------------------------

def bonferroni(p_values, alpha=0.05):
    """
    Bonferroni family-wise correction: reject only p <= alpha / n_tested.

    Reading the raw "dAIC > 2" node count from `compare_maps` as if every node
    were an independent, pre-specified test is exactly the mistake Marzocchi
    et al. (2020) warn against: at alpha = 0.05 you expect 5% of nodes to look
    significant by chance even when nothing changed, and a dense grid tests
    thousands of them. Bonferroni controls the probability of *any* false
    positive across the whole grid - conservative, but simple to defend.

    Args:
        p_values (array-like): p-values (here, Utsu Pb from `compare_maps`),
            any shape, NaN allowed where untested
        alpha (float): family-wise error rate
    Returns:
        dict: n_tested, n_significant, threshold_p, reject (bool array, same
              shape as `p_values`, False where NaN or not significant)
    """
    values = np.asarray(p_values, dtype=float)
    finite = np.isfinite(values)
    n = int(finite.sum())
    reject = np.zeros(values.shape, dtype=bool)
    if n == 0:
        return {"n_tested": 0, "n_significant": 0, "threshold_p": np.nan,
                "reject": reject}
    threshold = alpha / n
    reject[finite] = values[finite] <= threshold
    return {"n_tested": n, "n_significant": int(reject.sum()),
            "threshold_p": float(threshold), "reject": reject}


def benjamini_hochberg(p_values, alpha=0.05):
    """
    Benjamini-Hochberg false-discovery-rate correction.

    Less conservative than `bonferroni`: controls the *expected proportion* of
    significant nodes that are false discoveries, rather than the probability
    of any false positive at all. Both are legitimate; Marzocchi et al. (2020)
    do not mandate one, only that some correction replaces the raw count.

    Args:
        p_values (array-like): p-values, any shape, NaN allowed where untested
        alpha (float): target false discovery rate
    Returns:
        dict: n_tested, n_significant, threshold_p (the largest p that still
              passes), reject (bool array, same shape as `p_values`)
    """
    values = np.asarray(p_values, dtype=float)
    finite = np.isfinite(values)
    flat = values[finite]
    n = flat.size
    reject = np.zeros(values.shape, dtype=bool)
    if n == 0:
        return {"n_tested": 0, "n_significant": 0, "threshold_p": np.nan,
                "reject": reject}

    order = np.argsort(flat)
    sorted_p = flat[order]
    ranks = np.arange(1, n + 1)
    passing = np.where(sorted_p <= ranks / n * alpha)[0]
    if passing.size == 0:
        return {"n_tested": n, "n_significant": 0, "threshold_p": np.nan,
                "reject": reject}

    cutoff = int(passing.max())
    flat_reject = np.zeros(n, dtype=bool)
    flat_reject[order[:cutoff + 1]] = True
    reject[finite] = flat_reject
    return {"n_tested": n, "n_significant": cutoff + 1,
            "threshold_p": float(sorted_p[cutoff]), "reject": reject}


def effective_samples(radius_km, region_area_km2=None, length_km=None,
                      width_km=None):
    """
    How many genuinely independent samples a dense grid's area actually holds.

    Neighbouring nodes spaced closer than the sampling radius share most of
    their events, so the node count from a dense grid vastly overstates the
    number of independent tests - the root cause of the multiple-testing
    problem `bonferroni`/`benjamini_hochberg` correct for. This estimates the
    honest count as non-overlapping disks (map view) or non-overlapping
    along-section spans (section view) tiling the same area, i.e. what you
    would get by literally re-sampling on a grid spaced 2*radius apart.

    This is an estimate, not a substitute for actually re-running the
    comparison on a non-overlapping grid (which is what step 5 also does, and
    is the more defensible number for a headline claim) - use this to sanity
    check that number, not instead of it.

    Args:
        radius_km (float): sampling radius
        region_area_km2 (float): map-view area; give this OR length_km
        length_km (float): section length (section view)
        width_km (float): section width (section view; sampling is 2-D along
            the section, not 3-D, so only length tiles - width is the sampling
            volume's other axis already, not additional independent area)
    Returns:
        int: estimated number of independent samples
    """
    if region_area_km2 is not None:
        return max(1, int(region_area_km2 / (np.pi * radius_km ** 2)))
    if length_km is not None:
        return max(1, int(length_km / (2 * radius_km)))
    raise ValueError("give region_area_km2 or length_km")
