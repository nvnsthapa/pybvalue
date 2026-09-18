"""
Step 1 - acquire an earthquake catalog.

Two paths, one output. Either download fresh from an FDSN event service, or
normalise a catalog you already have; both write the same canonical CSV plus a
JSON provenance sidecar into output/<study>/1_catalog/, so every later step
reads one format and never needs to know which path was taken. CONFIG['study']
groups a related set of regions under one output tree (default
"parkfield_salton"; use a new name, e.g. "elsalvador", to keep a different
study's outputs separate).

Also writes a raw reconnaissance figure, <name>_catalog.png (+ a vector
<name>_catalog.pdf for print): epicenters on a map, and magnitude through
time. There is no geometry yet at this stage (no bbox or cross section has
been chosen), so it is deliberately just the two things you need to choose
one - unlike step 2's diagnostics, which assume a geometry already exists.

The map panel has a real basemap under it (Esri World Topo: shaded relief,
political boundaries, place names - one no-auth tile source, not plain
OpenStreetMap, whose own tile server blocks scripted access per its usage
policy). Needs `pip install contextily` (pulls in rasterio); without it, or
without network, or if the tile service is down, the map still draws, just
without the background - a missing basemap is never fatal to step 1. Tiles
are cached to output/.tile_cache/ so repeat runs do not re-fetch.

No command-line flags: edit CONFIG below, then

    python scripts/1_downloadCatalog.py

Set exactly one of CONFIG["region"] (a name from REGIONS, below) or
CONFIG["file"] (a path to a catalog already on disk). Examples of each,
as what CONFIG should look like:

    # fresh download, named region
    CONFIG["region"] = "parkfield"

    # fresh download, explicit box (region left as None)
    CONFIG["datacenter"] = "NCEDC"
    CONFIG["bbox"] = (35.6, 36.3, -120.8, -120.0)
    CONFIG["start"] = "1981-01-01"
    CONFIG["name"] = "parkfield"

    # a catalog you already have (format auto-detected)
    CONFIG["file"] = "data/1977_20240611.txt"
    CONFIG["name"] = "socal_scedc"
"""

import os
import sys
from types import SimpleNamespace

import numpy as np
import pyproj

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib                                                  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                    # noqa: E402
import matplotlib.dates as mdates                                  # noqa: E402
from matplotlib.lines import Line2D                                # noqa: E402
from matplotlib.ticker import FuncFormatter                        # noqa: E402

try:
    import contextily as cx        # optional: basemap tiles for plot_catalog()
except ImportError:
    cx = None

from src import catalog  # noqa: E402

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# OUTPUT_DIR is computed in main() from CONFIG['study']: output/<study>/1_catalog

# Web Mercator (EPSG:3857) is what every web basemap tile is projected in, so
# the map panel of plot_catalog() plots in it directly rather than raw
# lon/lat - see that function for why.
_TO_MERCATOR = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
_FROM_MERCATOR = pyproj.Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)

# Named regions are a convenience only - the pipeline itself is area-agnostic
# and any box can be given via CONFIG['bbox']. 'parkfield' reproduces the study area
# of Schorlemmer et al. (2004), paper 1; see refs/METHOD.md.
REGIONS = {
    "parkfield": {
        "datacenter": "NCEDC",
        # Box wraps the paper's cross section P1 (36.4N, 121.0W) -> P2
        # (35.64N, 120.2W) with a small margin. Do not narrow it: a tighter
        # box clips P1, and a wider rectangle pulls in the 1983 Coalinga
        # sequence, which lies 35 km off-section on a different fault system
        # and pushes the regional b value from 0.92 down to ~0.70.
        "bbox": (35.55, 36.50, -121.10, -120.10),
        "cross_section": {"p1": (36.40, -121.00), "p2": (35.64, -120.20),
                          "width_km": 5.0},
        "start": "1981-01-01",
        "end": None,  # today
        "note": "Parkfield segment, San Andreas. NCEDC archives NCSN, the "
                "network Schorlemmer et al. (2004) used; SCEDC barely covers "
                "this box. Cross section P1->P2 is 111 km long, 5 km wide.",
    },
    "salton": {
        "datacenter": "SCEDC",
        "bbox": (32.0, 33.94, -117.16, -114.38),
        "start": "1981-01-01",
        "end": None,
        "note": "Salton Trough / Imperial Valley, the project target region.",
    },
}


# Edit this, then just run the script - no command-line flags.
# "source (choose one)": region OR file. Everything else is optional; None
# means "use the region preset" (download options) or the shown default
# (file options).
CONFIG = {
    # output/<study>/1_catalog/ - groups related regions under one output
    # tree; use a new name (e.g. "elsalvador") to keep a study separate.
    "study": "elsalvador",

    # -- source: exactly one of these two -----------------------------------
    "region": None,   # a name from REGIONS above, e.g. "salton"
    "file": "data/Elsalvador/slvsafull_local_filtered_cc.csv",   # already have it, no download

    # -- download options (used when "region" is set; None = preset's value)
    "datacenter": None,       # e.g. one of catalog.FDSN_SERVICES, or a full URL
    "bbox": None,              # (minlat, maxlat, minlon, maxlon), degrees
    "start": None,             # e.g. "1981-01-01"
    "end": None,               # default: now
    "minmag": None,
    "maxdepth": None,

    # -- file options (used when "file" is set) ------------------------------
    "format": "csv",          # generic CSV reader; columns matched by alias

    "name": "elsalvador",      # output basename (default: region or file stem)
    "force": False,            # overwrite an existing catalog of the same name
}


def main():
    args = SimpleNamespace(**CONFIG)
    OUTPUT_DIR = os.path.join(PROJECT_DIR, "output", args.study, "1_catalog")

    if not args.region and not args.file:
        print("error: set CONFIG['region'] or CONFIG['file'] at the top of "
              "this script.", file=sys.stderr)
        return 2
    if args.region and args.region not in REGIONS:
        print(f"error: CONFIG['region'] = {args.region!r} is not in REGIONS: "
              f"{sorted(REGIONS)}", file=sys.stderr)
        return 2

    # ---- resolve where the data is coming from ----------------------------
    if args.file:
        name = args.name or os.path.splitext(os.path.basename(args.file))[0]
        provenance = {"source": "file", "original_path": os.path.abspath(args.file),
                      "format": args.format}
    else:
        cfg = REGIONS[args.region]
        datacenter = args.datacenter or cfg["datacenter"]
        bbox = tuple(args.bbox) if args.bbox else cfg["bbox"]
        start = args.start or cfg["start"]
        end = args.end or cfg["end"] or pd_now()
        name = args.name or args.region
        provenance = {
            "source": "fdsn", "datacenter": datacenter,
            "bbox_minlat_maxlat_minlon_maxlon": list(bbox),
            "starttime": str(start), "endtime": str(end),
            "minmagnitude": args.minmag, "maxdepth": args.maxdepth,
        }
        if cfg.get("cross_section") and not args.bbox:
            provenance["cross_section"] = cfg["cross_section"]

    destination = os.path.join(OUTPUT_DIR, f"{name}.csv")
    if os.path.exists(destination) and not args.force:
        print(f"error: {os.path.relpath(destination, PROJECT_DIR)} already exists.\n"
              f"       Set CONFIG['force'] = True to overwrite, or a different CONFIG['name'].",
              file=sys.stderr)
        return 1

    # ---- fetch ------------------------------------------------------------
    if args.file:
        print(f"Reading {args.file}")
        events = catalog.read_file(args.file, fmt=args.format)
        provenance["detected_format"] = catalog.detect_format(args.file) \
            if args.format == "auto" else args.format
        # A regional catalog file usually covers far more ground than the study
        # area, so allow the same box (and magnitude/depth limits) to trim it
        # here rather than carrying every event through the pipeline.
        box = tuple(args.bbox) if args.bbox else (
            REGIONS[args.region]["bbox"] if args.region else None)
        if box:
            before = len(events)
            events = events[events["latitude"].between(box[0], box[1])
                            & events["longitude"].between(box[2], box[3])]
            print(f"  box {box}: kept {len(events):,} of {before:,}")
            provenance["bbox_minlat_maxlat_minlon_maxlon"] = list(box)
        for limit, column, op in (("minmag", "magnitude", "ge"),
                                  ("maxdepth", "depth", "le")):
            value = getattr(args, limit)
            if value is not None:
                before = len(events)
                keep = (events[column] >= value if op == "ge"
                        else events[column] <= value)
                events = events[keep]
                print(f"  {column} {'>=' if op == 'ge' else '<='} {value}: "
                      f"kept {len(events):,} of {before:,}")
    else:
        print(f"Region {name}: {bbox[0]}..{bbox[1]} N, {bbox[2]}..{bbox[3]} E")
        print(f"Time   {start} .. {end}")
        events = catalog.download_fdsn(
            datacenter=datacenter, starttime=start, endtime=end,
            minlatitude=bbox[0], maxlatitude=bbox[1],
            minlongitude=bbox[2], maxlongitude=bbox[3],
            minmagnitude=args.minmag, maxdepth=args.maxdepth,
        )

    if events.empty:
        print("No events matched.", file=sys.stderr)
        return 1

    # ---- report and save --------------------------------------------------
    print()
    catalog.summarise(events, name)
    catalog.write(events, destination, provenance)
    figure_png = os.path.join(OUTPUT_DIR, f"{name}_catalog.png")
    figure_pdf = os.path.join(OUTPUT_DIR, f"{name}_catalog.pdf")
    fig = plot_catalog(name, events)
    fig.savefig(figure_png, dpi=300, bbox_inches="tight")
    fig.savefig(figure_pdf, bbox_inches="tight")   # vector, for print/publication
    plt.close(fig)
    print(f"\nWrote {os.path.relpath(destination, PROJECT_DIR)}")
    print(f"      {os.path.relpath(os.path.splitext(destination)[0] + '.json', PROJECT_DIR)}")
    print(f"      {os.path.relpath(figure_png, PROJECT_DIR)}")
    print(f"      {os.path.relpath(figure_pdf, PROJECT_DIR)}")
    return 0


PUB_STYLE = {
    "font.size": 11, "axes.labelsize": 12, "axes.titlesize": 12,
    "xtick.labelsize": 10, "ytick.labelsize": 10, "legend.fontsize": 9,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "pdf.fonttype": 42, "ps.fonttype": 42,   # embed real text, not paths, in vector output
}


def _size_for_magnitude(magnitude, floor):
    """Marker area (points^2) for a magnitude - shared between data and legend markers."""
    return 8 + 10 * np.clip(magnitude - floor, 0, None)


def _add_scalebar(ax, latitude0):
    """
    A short km scale bar in the lower-left corner, sized to the panel's own
    span. `ax` is in Web Mercator metres, which overstate ground distance
    away from the equator by 1/cos(latitude) - correct for that so the bar
    is true to scale (~3% at El Salvador's latitude, ~20% at Parkfield's).
    """
    k = 1 / np.cos(np.radians(latitude0))
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    span_km = (x1 - x0) / k / 1000
    length_km = min((c for c in (1, 2, 5, 10, 20, 25, 50, 100, 200, 250)
                     if c >= span_km / 6), default=250)
    length_m = length_km * 1000 * k
    bar_x = x0 + 0.06 * (x1 - x0)
    bar_y = y0 + 0.07 * (y1 - y0)
    ax.plot([bar_x, bar_x + length_m], [bar_y, bar_y], color="0.1", lw=2.2,
            solid_capstyle="butt", zorder=6)
    ax.text(bar_x + length_m / 2, bar_y + 0.018 * (y1 - y0), f"{length_km:g} km",
            ha="center", va="bottom", fontsize=9, zorder=6)


def _add_basemap(ax):
    """
    Esri World Topo: shaded relief, political boundaries and place names in
    one no-auth tile source (contextily's default, plain OpenStreetMap tiles,
    actively blocks scripted/automated access per its own tile usage policy -
    this does not). Tiles are cached to disk so repeat runs do not re-fetch.

    Never fatal: missing `contextily`, no network, or a tile service outage
    should not stop step 1 from producing its CSV, so failure is reported and
    the map falls back to a plain background rather than raising.
    """
    if cx is None:
        print("  ! contextily not installed, skipping basemap "
              "(pip install contextily)", file=sys.stderr)
        return False
    try:
        cx.set_cache_dir(os.path.join(PROJECT_DIR, "output", ".tile_cache"))
        cx.add_basemap(ax, crs="EPSG:3857", source=cx.providers.Esri.WorldTopoMap,
                       attribution=False, zorder=0)
        return True
    except Exception as exc:   # noqa: BLE001 - any failure here must stay non-fatal
        print(f"  ! basemap unavailable ({type(exc).__name__}: {exc}); "
              f"plotting without one", file=sys.stderr)
        return False


def plot_catalog(name, events):
    """
    Raw reconnaissance figure, before any geometry exists: epicenters on a
    map, and magnitude through time. Events with no magnitude (kept, not
    dropped - see src/catalog.py) are drawn separately so they do not vanish
    silently or break the colour scale.
    """
    has_mag = events["magnitude"].notna()
    dated, undated = events[has_mag], events[~has_mag]
    floor = float(dated.magnitude.min())

    # Plotted in Web Mercator (what every web basemap tile is projected in),
    # not raw lon/lat, so the tile image lines up with the events. Tick
    # labels are converted back to degrees below, so this is invisible to
    # the reader - the map just gains true north-up geometry and a real
    # background (relief, political boundaries, place names).
    dated_x, dated_y = _TO_MERCATOR.transform(dated.longitude.values, dated.latitude.values)
    if len(undated):
        undated_x, undated_y = _TO_MERCATOR.transform(undated.longitude.values,
                                                       undated.latitude.values)

    # Edge-stroked, size-encoded bubbles read well for a few thousand events
    # (El Salvador) but overplot into a solid blob well before catalogs the
    # size of Parkfield or Salton (tens to hundreds of thousands) - drop the
    # edge and shrink/flatten the size ramp once there are too many to resolve
    # individually, the same way step 2's map view already does.
    dense = len(dated) > 20_000
    marker_kw = (dict(lw=0, edgecolor="none", alpha=.55) if dense
                else dict(lw=.4, edgecolor="0.15", alpha=.85))
    size_scale = .3 if dense else 1.0

    with plt.rc_context(PUB_STYLE):
        fig, (map_ax, mag_ax) = plt.subplots(1, 2, figsize=(13, 6),
                                             constrained_layout=True)

        # --- map -------------------------------------------------------
        legend_handles = []
        if len(undated):
            map_ax.scatter(undated_x, undated_y,
                           s=10 if not dense else 4,
                           c="0.9", zorder=2,
                           **(dict(lw=.5, edgecolor="0.2") if not dense else dict(lw=0)))
            legend_handles.append(Line2D([0], [0], marker="o", color="none",
                                         markerfacecolor="0.9", markeredgecolor="0.2",
                                         markeredgewidth=.5, markersize=5,
                                         label=f"no magnitude ({len(undated):,})"))
        scatter = map_ax.scatter(dated_x, dated_y,
                                 s=_size_for_magnitude(dated.magnitude, floor) * size_scale,
                                 c=dated.magnitude, cmap="viridis", zorder=3, **marker_kw)
        cbar = fig.colorbar(scatter, ax=map_ax, label="magnitude", pad=.02, shrink=.85)
        cbar.outline.set_linewidth(.6)

        legend_mags = sorted(set(np.round(np.linspace(floor, dated.magnitude.max(), 3) * 2) / 2))
        legend_handles += [
            Line2D([0], [0], marker="o", color="none", markerfacecolor="0.6",
                   markeredgecolor="0.15", markeredgewidth=.4,
                   markersize=np.sqrt(_size_for_magnitude(m, floor) * size_scale),
                   label=f"M {m:g}")
            for m in legend_mags
        ]
        map_ax.legend(handles=legend_handles, loc="upper right", framealpha=.92,
                      title="magnitude", borderpad=.8, labelspacing=.8)

        map_ax.set_xlabel("longitude"), map_ax.set_ylabel("latitude")
        map_ax.set_aspect("equal")   # Mercator x/y already share one scale
        map_ax.xaxis.set_major_formatter(
            FuncFormatter(lambda x, _pos: f"{_FROM_MERCATOR.transform(x, 0)[0]:.2f}"))
        map_ax.yaxis.set_major_formatter(
            FuncFormatter(lambda y, _pos: f"{_FROM_MERCATOR.transform(0, y)[1]:.2f}"))
        map_ax.set_title(f"{name}  ·  {len(events):,} events  ·  "
                         f"{events.time.min():%Y-%m-%d} to {events.time.max():%Y-%m-%d}",
                         fontsize=10.5, loc="left")
        _add_basemap(map_ax)
        _add_scalebar(map_ax, float(dated.latitude.median()))
        map_ax.text(-0.10, 1.0, "a", transform=map_ax.transAxes,
                    fontsize=13, fontweight="bold", va="top")

        # --- magnitude through time --------------------------------------
        mag_ax.scatter(dated.time, dated.magnitude, s=6, c="steelblue", lw=.2,
                       edgecolor="0.2", alpha=.65, zorder=3)
        mag_ax.set_xlabel("time"), mag_ax.set_ylabel("magnitude")
        mag_ax.set_title("magnitude through time", fontsize=10.5, loc="left")
        mag_ax.grid(alpha=.25, lw=.5)
        for spine in ("top", "right"):
            mag_ax.spines[spine].set_visible(False)
        mag_ax.text(-0.10, 1.0, "b", transform=mag_ax.transAxes,
                    fontsize=13, fontweight="bold", va="top")
        # fig.autofmt_xdate() fights constrained_layout (rotation gets reset
        # on the next layout pass) - rotate this axis's own ticks instead.
        mag_ax.xaxis.set_major_locator(mdates.AutoDateLocator(maxticks=8))
        plt.setp(mag_ax.get_xticklabels(), rotation=30, ha="right")

    return fig


def pd_now():
    import pandas as pd
    return pd.Timestamp.utcnow().tz_localize(None).strftime("%Y-%m-%d")


if __name__ == "__main__":
    sys.exit(main())
