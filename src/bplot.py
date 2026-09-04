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


def bvalue_norm(centre, half_range=0.45):
    """
    Diverging norm centred on the regional b value.

    Args:
        centre (float): regional b - the value that should read as neutral
        half_range (float): span either side of centre
    Returns:
        TwoSlopeNorm
    """
    return TwoSlopeNorm(vmin=centre - half_range, vcenter=centre,
                        vmax=centre + half_range)


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
