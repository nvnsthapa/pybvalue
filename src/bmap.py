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
