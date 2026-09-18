"""
Rendering for b-value maps.

Colour choices are deliberate and shared by every figure in the project:

  * b-value is **diverging**, not sequential. The quantity that carries meaning
    is which side of the regional average a node sits on, so the ramp is two
    opposing hues with a **neutral grey midpoint** pinned to the regional b.
    A hue at the midpoint would make "average" look like a value of its own.
  * Warm = low b, cool = high b, matching the usual reading of low b as
    relatively more large events.
  * The published figures use a rainbow ramp. That is a 2004 convention and is
    not reproduced: a multi-hue ramp for magnitude invents boundaries the data
    does not have and is unreadable to colour-vision-deficient readers.
  * Every panel of a small-multiple set shares one scale, otherwise the panels
    cannot be compared - which is the entire point of the radius scan.
  * Unresolved nodes (below Nmin) get a flat neutral that reads as absent, and
    are never left to blend into the low end of the ramp.
"""

import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.ticker import FuncFormatter

# warm -> neutral grey -> cool
B_COLORS = ["#8c2d16", "#c1613c", "#e0a180", "#d8d6d2", "#8fb3c8", "#3d7ea6", "#17445f"]
NO_DATA = "#f2f1ef"
DELTA_COLORS = ["#762a5a", "#a86598", "#d9bcd0", "#d8d6d2", "#b9cbb2", "#7da26f", "#3f6b34"]


def bvalue_cmap():
    """Diverging warm-grey-cool ramp for b-values."""
    cmap = LinearSegmentedColormap.from_list("bvalue", B_COLORS)
    cmap.set_bad(NO_DATA)
    return cmap


def delta_cmap():
    """Diverging ramp for b-value differences, distinct in hue from bvalue_cmap."""
    cmap = LinearSegmentedColormap.from_list("delta_b", DELTA_COLORS)
    cmap.set_bad(NO_DATA)
    return cmap


def sequential_cmap(kind="sigma"):
    """
    Single-hue light-to-dark ramps for magnitude quantities.

    Sequential data gets one hue: a multi-hue ramp would invent category
    boundaries the numbers do not have. `sigma` (uncertainty) is a cool slate,
    `count` a neutral, `probability` a warm ramp where darker reads as more
    likely.

    Args:
        kind (str): 'sigma', 'count' or 'probability'
    Returns:
        LinearSegmentedColormap
    """
    ramps = {
        "sigma": ["#f0efed", "#c3cdd4", "#8ba3b3", "#547a91", "#2b566e", "#14384b"],
        "count": ["#f2f1ef", "#d4d1cb", "#aaa59c", "#7d776d", "#544e45", "#2e2a24"],
        "probability": ["#fbf3ec", "#f3d3bb", "#e8a982", "#d67a4f", "#b04e2a", "#7a2d13"],
    }
    if kind not in ramps:
        raise ValueError(f"unknown ramp {kind!r}; choose from {sorted(ramps)}")
    cmap = LinearSegmentedColormap.from_list(f"seq_{kind}", ramps[kind])
    cmap.set_bad(NO_DATA)
    return cmap


def bvalue_norm(centre, values=None, half_range=0.45, percentile=5):
    """
    Diverging norm centred on the regional/mapped b value.

    A fixed half_range either saturates the colour scale (real spread wider
    than guessed - every node past the edge reads as the same extreme colour)
    or wastes most of it (real spread narrower - everything crowds into a
    sliver near the centre). Passing `values` fixes both: vmin/vmax become
    the `percentile`/`100-percentile` of the data actually being mapped, the
    same robust-to-one-freak-node approach step 3's radius scan already uses
    for its own colour range.

    Args:
        centre (float): the value that should read as neutral
        values (array-like): mapped values (b and b-positive combined, if both
            share this scale) to set a data-driven vmin/vmax from; None keeps
            the old fixed-width behaviour
        half_range (float): span either side of centre, used only when
            `values` is not given (or has fewer than 2 finite values)
        percentile (float): low percentile for vmin; high is 100 - this
    Returns:
        TwoSlopeNorm
    """
    if values is not None:
        finite = np.asarray(values, dtype=float)
        finite = finite[np.isfinite(finite)]
        if finite.size >= 2:
            vmin, vmax = np.percentile(finite, [percentile, 100 - percentile])
            # A centre outside [vmin, vmax] (possible if it's the regional fit,
            # biased by the same incompleteness the map is correcting for)
            # would make TwoSlopeNorm raise - widen to include it rather than
            # silently falling back to the fixed-width range.
            vmin, vmax = min(vmin, centre), max(vmax, centre)
            if vmin < vmax:
                return TwoSlopeNorm(vmin=vmin, vcenter=centre, vmax=vmax)
    return TwoSlopeNorm(vmin=centre - half_range, vcenter=centre,
                        vmax=centre + half_range)


def set_lonlat_ticks(ax, region):
    """
    Relabel a map-view axis's km ticks as longitude/latitude.

    Nodes are sampled and plotted in a locally projected km frame (a true
    circle of radius r has to be a circle on the ground - see geometry.py),
    but km-from-centre is a much harder number to place on a real map than
    a longitude or latitude, so the tick *labels* are converted back through
    `region.to_lonlat` while the underlying data stays in km.

    The azimuthal equidistant projection `region` uses is not separable, so
    each x tick is converted at y=0 (and each y tick at x=0) rather than at
    its true row/column - checked against the full joint conversion, this is
    accurate to under 1 km even at Salton's ~200 km span, which is exact
    enough for a tick label.

    Args:
        ax (matplotlib.axes.Axes): target axis, x/y already in region km
        region (geometry.MapRegion): supplies the inverse projection
    """
    ax.xaxis.set_major_formatter(
        FuncFormatter(lambda x, _pos: f"{region.to_lonlat(x, 0)[0]:.2f}"))
    ax.yaxis.set_major_formatter(
        FuncFormatter(lambda y, _pos: f"{region.to_lonlat(0, y)[1]:.2f}"))


def draw_section_map(ax, result, norm, cmap=None, key="b", show_ylabel=True,
                     show_xlabel=True):
    """
    Draw one b-value section onto an axis.

    Args:
        ax (matplotlib.axes.Axes): target axis
        result (dict): output of bmap.map_section
        norm: colour norm, shared across panels
        cmap: colormap; defaults to bvalue_cmap()
        key (str): which field of `result` to draw
        show_xlabel, show_ylabel (bool): label only the outer panels of a grid
    Returns:
        the QuadMesh, for a colorbar
    """
    cmap = cmap or bvalue_cmap()
    field = np.ma.masked_invalid(result[key])
    mesh = ax.pcolormesh(result["along"], result["depth"], field,
                         cmap=cmap, norm=norm, shading="nearest", rasterized=True)
    ax.set_xlim(result["along"].min(), result["along"].max())
    ax.set_ylim(result["depth"].max(), result["depth"].min())   # depth downward
    if show_xlabel:
        ax.set_xlabel("distance along section from P1 (km)")
    if show_ylabel:
        ax.set_ylabel("depth (km)")
    ax.tick_params(labelsize=9)
    for spine in ax.spines.values():
        spine.set_color("#b9b6b1")
    return mesh


def draw_region_map(ax, result, norm, cmap=None, key="b", show_xlabel=True,
                    show_ylabel=True, region=None):
    """
    Draw one map-view field onto an axis, in the projected km frame.

    The aspect is locked to equal: a b-value anomaly's shape is part of the
    evidence, and a stretched axis invents elongation that is not in the data.

    Args:
        ax (matplotlib.axes.Axes): target axis
        result (dict): output of bmap.map_region
        norm: colour norm, shared across panels
        cmap: colormap; defaults to bvalue_cmap()
        key (str): which field of `result` to draw
        show_xlabel, show_ylabel (bool): label only the outer panels
        region (geometry.MapRegion): if given, tick labels read longitude/
            latitude (see set_lonlat_ticks) instead of km from centre
    Returns:
        the QuadMesh, for a colorbar
    """
    cmap = cmap or bvalue_cmap()
    field = np.ma.masked_invalid(result[key])
    mesh = ax.pcolormesh(result["x"], result["y"], field, cmap=cmap, norm=norm,
                         shading="nearest", rasterized=True)
    ax.set_aspect("equal")
    if region is not None:
        set_lonlat_ticks(ax, region)
    if show_xlabel:
        ax.set_xlabel("longitude" if region is not None else "east of region centre (km)")
    if show_ylabel:
        ax.set_ylabel("latitude" if region is not None else "north of region centre (km)")
    ax.tick_params(labelsize=9)
    for spine in ax.spines.values():
        spine.set_color("#b9b6b1")
    return mesh


def overlay_wells(ax, region, wells, size=26):
    """
    Mark geothermal wells on a map-view axis, projected into the same km frame.

    `wells` is not assumed to belong to this region - it is one fixed dataset
    (currently just the Salton Sea Geothermal Field) loaded unconditionally
    for any map-view run, so a region far away (El Salvador, say) must not
    plot California wells thousands of km off in its local frame and blow
    the axis limits out to include them. Silently keeping only the wells
    that actually fall inside `region.bbox` makes that the caller's problem
    to never worry about, rather than something every new region has to
    remember to guard against.

    Args:
        ax: target axis
        region (geometry.MapRegion): supplies the projection and bbox
        wells (pd.DataFrame): needs `Latitude`, `Longitude`; `Well Type` optional
        size (float): marker size
    Returns:
        bool: whether any well fell inside the region and was plotted - callers
              use this to skip an otherwise-empty legend entry
    """
    if wells is None or not len(wells):
        return False
    min_lat, max_lat, min_lon, max_lon = region.bbox
    inside = wells[wells["Latitude"].between(min_lat, max_lat)
                  & wells["Longitude"].between(min_lon, max_lon)]
    if not len(inside):
        return False
    placed = region.project(inside.rename(columns={"Latitude": "latitude",
                                                    "Longitude": "longitude"}))
    ax.scatter(placed["x"], placed["y"], s=size, marker="^",
               facecolor="#f5f0e8", edgecolor="#1a1a1a", linewidth=.8,
               zorder=7, label=f"geothermal wells ({len(placed)})")
    return True


def annotate_panel(ax, text, coverage=None):
    """Label a panel with its radius and, optionally, its coverage."""
    ax.text(0.012, 0.94, text, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="left", color="#1a1a1a",
            bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="#b9b6b1", alpha=.92))
    if coverage is not None:
        ax.text(0.988, 0.94, f"{coverage * 100:.0f}% covered",
                transform=ax.transAxes, fontsize=8.5, va="top", ha="right",
                color="#4a4a4a",
                bbox=dict(boxstyle="round,pad=0.24", fc="white",
                          ec="#d3d0cb", alpha=.88))


def mark_features(ax, features, ymax=None):
    """
    Mark named locations along the section (asperity, creeping section, ...).

    Args:
        ax: target axis
        features (dict): {label: (along_km, depth_km)}
        ymax (float): unused, kept for signature stability
    """
    for label, (along, depth) in features.items():
        ax.plot(along, depth, marker="o", ms=8, mfc="none", mec="#1a1a1a",
                mew=1.6, zorder=6)
        ax.annotate(label, (along, depth), textcoords="offset points",
                    xytext=(9, -3), fontsize=8.5, color="#1a1a1a", zorder=6,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white",
                              ec="none", alpha=.75))
