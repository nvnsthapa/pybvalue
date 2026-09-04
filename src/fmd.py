"""
Frequency-magnitude distribution: completeness, b-value, and the Utsu test.

The estimators themselves come from tgoebel/magnitude-distribution
(`src/FMD_GR.py`), cloned alongside this project. That file is loaded by path
rather than imported normally, because both repositories ship a package called
`src` and a plain import would resolve to the wrong one.

What this module adds on top is the part that is easy to get wrong: choosing
the bin correction to match the catalog's own magnitude convention. See
refs/METHOD.md, "The bin-correction trap".
"""

import importlib.util
import os

import numpy as np

_FMD_CLASS = None
_SEARCH_PATHS = [
    "magnitude-distribution/src/FMD_GR.py",
    "../magnitude-distribution/src/FMD_GR.py",
]


def load_fmd():
    """
    Return the FMD class from the magnitude-distribution clone.

    Returns:
        type: the FMD class
    Raises:
        ImportError: if the clone is missing
    """
    global _FMD_CLASS
    if _FMD_CLASS is not None:
        return _FMD_CLASS

    project = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    for relative in _SEARCH_PATHS:
        candidate = os.path.normpath(os.path.join(project, relative))
        if os.path.exists(candidate):
            spec = importlib.util.spec_from_file_location("FMD_GR", candidate)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            _FMD_CLASS = module.FMD
            return _FMD_CLASS

    raise ImportError(
        "Could not find magnitude-distribution/src/FMD_GR.py.\n"
        "Clone it into the project root:\n"
        "    git clone https://github.com/tgoebel/magnitude-distribution"
    )


# ---------------------------------------------------------------------------
# magnitude binning
# ---------------------------------------------------------------------------

def detect_binsize(magnitudes, candidates=(0.1, 0.05, 0.01)):
    """
    Infer the magnitude bin width a catalog was reported on.

    Network catalogs quantise magnitudes (NCEDC floors to 0.1); aggregators
    like USGS/ComCat report full precision. Which one you have decides the bin
    correction, and getting it wrong biases b by several percent.

    Args:
        magnitudes (array-like): observed magnitudes
        candidates (tuple): bin widths to test, largest first
    Returns:
        float: the bin width, or 0.0 if the magnitudes look continuous
    """
    values = np.asarray(magnitudes, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 0.0
    for width in candidates:
        residual = np.abs(values / width - np.round(values / width))
        if (residual < 1e-6).mean() > 0.99:
            return float(width)
    return 0.0


def bin_correction(magnitudes, binsize=None):
    """
    The correction to subtract from Mc before the MLE fit.

    Schorlemmer et al. use Mmin = Mc - dM/2 to undo the bias from rounding
    magnitudes into bins. A continuous catalog has no such bias, so its
    correction is zero.

    Args:
        magnitudes (array-like): observed magnitudes
        binsize (float): override the detected bin width
    Returns:
        float: the correction in magnitude units
    """
    if binsize is None:
        binsize = detect_binsize(magnitudes)
    return binsize / 2.0


# ---------------------------------------------------------------------------
# estimators
# ---------------------------------------------------------------------------

def fit_gr(magnitudes, mc, binsize=None, min_events=2):
    """
    Aki maximum-likelihood b-value with Shi & Bolt uncertainty.

        b       = log10(e) / (mean(M) - (Mc - dM/2))
        sigma_b = 2.30 b^2 sqrt( sum (Mi - meanM)^2 / (n(n-1)) )
        a       = log10(n) + b * Mc

    Args:
        magnitudes (array-like): magnitudes of the sample (any completeness)
        mc (float): magnitude of completeness
        binsize (float): magnitude bin width; detected when None
        min_events (int): below this many events above Mc, return NaNs
    Returns:
        dict: b, sigma, a, n, mc, binsize
    """
    values = np.asarray(magnitudes, dtype=float)
    values = values[np.isfinite(values)]
    if binsize is None:
        binsize = detect_binsize(values)

    # Half a bin of slack. Binned magnitudes are bin centres, so the Mc bin
    # itself belongs in the sample; without the tolerance a value that is
    # 1.2999999999999998 after arithmetic drops the whole lowest bin and
    # biases b low. Clean data is unaffected - the next bin down is a full
    # bin away.
    above = values[values >= mc - (binsize / 2 if binsize > 0 else 0.0)]
    n = above.size
    blank = {"b": np.nan, "sigma": np.nan, "a": np.nan, "n": int(n),
             "mc": float(mc), "binsize": float(binsize)}
    if n < max(2, min_events):
        return blank

    mean_mag = above.mean()
    denominator = mean_mag - (mc - binsize / 2.0)
    if denominator <= 0:
        return blank

    b = np.log10(np.e) / denominator
    sigma = 2.30 * b**2 * np.sqrt(((above - mean_mag) ** 2).sum() / (n * (n - 1)))
    return {"b": float(b), "sigma": float(sigma),
            "a": float(np.log10(n) + b * mc), "n": int(n),
            "mc": float(mc), "binsize": float(binsize)}


def fit_gr_positive(magnitudes, binsize=None, dmc=None, mc=None, min_events=2):
    """
    b-positive: b from the positive differences between successive magnitudes.

    If magnitudes above completeness are exponentially distributed, the
    difference between two successive ones is Laplace distributed, and the
    positive half of that is exponential with the *same* decay rate. Because
    the estimator conditions on a difference rather than on an absolute
    magnitude, a detection threshold that moves around in time largely cancels:
    b+ stays close to the true b even where the catalog is incomplete.

    That is what makes it worth computing alongside the classic estimator.
    Where the two agree, the sample is complete and the classic value can be
    trusted; where they diverge, completeness is the suspect rather than
    tectonics.

    !! `magnitudes` must be in TIME ORDER. The estimator differences events
    that are consecutive in time, so a magnitude-sorted array silently returns
    a meaningless number rather than an error.

    Args:
        magnitudes (array-like): magnitudes, ordered in time
        binsize (float): magnitude bin width; detected from the data when None
        dmc (float): smallest difference to keep; defaults to binsize
        mc (float): optional pre-cut. b-positive tolerates a far lower cut than
            the classic estimator - that is the point - so leave it None or set
            it well below the apparent completeness
        min_events (int): minimum number of positive differences
    Returns:
        dict: b, sigma, a, n (number of differences used), mc, binsize,
              plus n_events and dmc
    Source:
        van der Elst (2021), JGR Solid Earth, doi:10.1029/2020JB021027
    """
    values = np.asarray(magnitudes, dtype=float)
    values = values[np.isfinite(values)]
    if binsize is None:
        binsize = detect_binsize(values)
    if dmc is None:
        dmc = binsize
    if dmc < 0:
        raise ValueError("dmc must be >= 0")

    if mc is not None:
        values = values[values >= mc - binsize / 2]

    blank = {"b": np.nan, "sigma": np.nan, "a": np.nan, "n": 0,
             "mc": float(mc) if mc is not None else np.nan,
             "binsize": float(binsize), "n_events": int(values.size),
             "dmc": float(dmc)}
    if values.size < 2:
        return blank

    differences = np.diff(values)
    if binsize > 0:
        # Subtracting two binned magnitudes lands just off the grid
        # (0.1 - 0.0 gives 0.09999999999999998), which then falls on the wrong
        # side of every later comparison. Snap back onto it.
        differences = np.round(differences / binsize) * binsize
    # Half a bin of slack so that magnitudes differing by exactly one bin are
    # kept rather than lost to a floating-point comparison.
    positive = differences[differences > dmc - binsize / 2]
    if positive.size < max(2, min_events):
        return blank

    # The positive differences are exponential above dmc, so the ordinary
    # estimator applies to them with dmc playing the role of completeness.
    fit = fit_gr(positive, dmc, binsize=binsize, min_events=min_events)
    fit["n_events"] = int(values.size)
    fit["dmc"] = float(dmc)
    fit["mc"] = float(mc) if mc is not None else np.nan
    return fit


def beta_tinti(magnitudes, mc, binsize):
    """
    Exact maximum-likelihood b for binned magnitudes (Tinti & Mulargia 1987).

        p    = 1 + dM / mean(M - Mc)
        beta = ln(p) / dM,        b = beta / ln(10)

    Aki's estimator with the Mc - dM/2 correction, which is what Schorlemmer et
    al. use and what `fit_gr` implements, is the first-order expansion of this.
    Provided for cross-checking that the approximation is not costing anything
    at a given bin width.

    Args:
        magnitudes (array-like): magnitudes at or above mc
        mc (float): completeness
        binsize (float): magnitude bin width; 0 for continuous magnitudes
    Returns:
        float: b-value, or NaN for an unusable sample
    """
    values = np.asarray(magnitudes, dtype=float)
    values = values[np.isfinite(values) & (values >= mc)]
    if values.size < 2:
        return np.nan
    mean_excess = (values - mc).mean()
    if mean_excess <= 0:
        return np.nan
    if binsize > 0:
        beta = np.log(1 + binsize / mean_excess) / binsize
    else:
        beta = 1.0 / mean_excess
    return float(beta / np.log(10))


def mc_maxcurv(magnitudes, binsize=0.1):
    """
    Magnitude of completeness by maximum curvature: the magnitude bin holding
    the most events in the non-cumulative distribution.

    Args:
        magnitudes (array-like): observed magnitudes
        binsize (float): histogram bin width
    Returns:
        float: Mc, or NaN for an empty sample
    """
    values = np.asarray(magnitudes, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan
    edges = np.arange(values.min() - binsize, values.max() + 2 * binsize, binsize)
    counts, _ = np.histogram(values, edges)
    if counts.sum() == 0:
        return np.nan
    return float(edges[:-1][counts.argmax()] + binsize / 2)


def mc_vs_time(times, magnitudes, window=500, step=10, binsize=0.1):
    """
    Completeness as a function of time, following Wiemer & Wyss (2000):
    a sliding window of fixed event count, Mc by maximum curvature in each.

    Schorlemmer et al. use window=500, step=10 to justify the choice of a
    single homogeneous Mc for the whole catalog.

    Args:
        times (array-like): origin times, sorted or not
        magnitudes (array-like): magnitudes, same order as times
        window (int): events per window
        step (int): events to advance between windows
        binsize (float): magnitude bin width
    Returns:
        tuple: (window centre times, Mc per window)
    """
    times = np.asarray(times)
    magnitudes = np.asarray(magnitudes, dtype=float)
    order = np.argsort(times)
    times, magnitudes = times[order], magnitudes[order]
    if times.size < window:
        return np.array([]), np.array([])

    centres, values = [], []
    for start in range(0, times.size - window + 1, step):
        chunk = magnitudes[start:start + window]
        centres.append(times[start + window // 2])
        values.append(mc_maxcurv(chunk, binsize))
    return np.asarray(centres), np.asarray(values)


def utsu_test(n1, n2, b1, b2):
    """
    Probability that two b-values come from the same population (Utsu, 1992).

        dAIC = -2(N1+N2) ln(N1+N2) + 2N1 ln(N1 + N2 b1/b2)
                                   + 2N2 ln(N1 b2/b1 + N2) - 2
        Pb   = exp(-dAIC/2 - 2)

    dAIC < 2 is not significant; dAIC > 2 is significant (Pb <= 0.05);
    dAIC > 5 is highly significant (Pb <= 0.01).

    Args:
        n1, n2 (int): sample sizes above Mc
        b1, b2 (float): the two b-values
    Returns:
        tuple: (Pb, dAIC)
    """
    if not all(np.isfinite([n1, n2, b1, b2])) or min(n1, n2) < 2 \
            or b1 <= 0 or b2 <= 0:
        return np.nan, np.nan
    total = n1 + n2
    d_aic = (-2 * total * np.log(total)
             + 2 * n1 * np.log(n1 + n2 * b1 / b2)
             + 2 * n2 * np.log(n1 * b2 / b1 + n2)
             - 2)
    return float(np.exp(-d_aic / 2 - 2)), float(d_aic)


def cumulative_fmd(magnitudes, binsize=0.1):
    """
    Binned and cumulative frequency-magnitude distributions.

    Args:
        magnitudes (array-like): observed magnitudes
        binsize (float): bin width
    Returns:
        tuple: (bin centres, counts per bin, cumulative counts N(>=M))
    """
    values = np.asarray(magnitudes, dtype=float)
    values = values[np.isfinite(values)]
    edges = np.arange(np.floor(values.min() * 10) / 10 - binsize / 2,
                      values.max() + 1.5 * binsize, binsize)
    counts, _ = np.histogram(values, edges)
    centres = edges[:-1] + binsize / 2
    return centres, counts, counts[::-1].cumsum()[::-1]
