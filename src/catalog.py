"""
Catalog acquisition and normalisation.

Every catalog entering the project - freshly downloaded or supplied by the
user - is converted to one canonical DataFrame so that no downstream step has
to know where the data came from:

    time       datetime64[ns]  origin time (UTC, tz-naive)
    latitude   float           degrees north
    longitude  float           degrees east
    depth      float           km, positive down
    magnitude  float
    magtype    str             Md, Ml, Mw, ... ('' if unknown)
    evid       str             event id, unique within a catalog

Two ways in:
    download_fdsn()  - fetch from an FDSN event service
    read_file()      - load a catalog already on disk (format auto-detected)
"""

import json
import os
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone

import numpy as np
import pandas as pd

COLUMNS = ["time", "latitude", "longitude", "depth", "magnitude", "magtype", "evid"]

# FDSN event services. Parkfield/NCSN lives at NCEDC; SCEDC only sees the
# southern edge of that box, so the datacenter choice is not cosmetic.
FDSN_SERVICES = {
    "NCEDC": "https://service.ncedc.org/fdsnws/event/1/query",
    "SCEDC": "https://service.scedc.caltech.edu/fdsnws/event/1/query",
    "USGS": "https://earthquake.usgs.gov/fdsnws/event/1/query",
    "EMSC": "https://www.seismicportal.eu/fdsnws/event/1/query",
    "INGV": "https://webservices.ingv.it/fdsnws/event/1/query",
}

# Most FDSN servers silently truncate at 10k rows rather than erroring, so a
# response of exactly this length is treated as "there may be more".
SERVER_ROW_CAP = 10000


# ---------------------------------------------------------------------------
# download
# ---------------------------------------------------------------------------

def _fdsn_request(url, params, timeout=300):
    """Issue one FDSN text query. Returns the body, or '' when no data."""
    query = f"{url}?{urllib.parse.urlencode(params)}"
    try:
        with urllib.request.urlopen(query, timeout=timeout) as response:
            return response.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError as err:
        if err.code in (204, 404):  # FDSN "no data matched"
            return ""
        raise


def _parse_fdsn_text(body):
    """Parse the pipe-delimited FDSN event text format into canonical columns."""
    lines = [ln for ln in body.splitlines() if ln.strip() and not ln.startswith("#")]
    if not lines:
        return pd.DataFrame(columns=COLUMNS)

    rows = [ln.split("|") for ln in lines]
    # EventID|Time|Latitude|Longitude|Depth/km|Author|Catalog|Contributor|
    # ContributorID|MagType|Magnitude|MagAuthor|EventLocationName|EventType
    frame = pd.DataFrame({
        "evid": [r[0].strip() for r in rows],
        "time": [r[1].strip() for r in rows],
        "latitude": [r[2] for r in rows],
        "longitude": [r[3] for r in rows],
        "depth": [r[4] for r in rows],
        "magtype": [r[9].strip() if len(r) > 9 else "" for r in rows],
        "magnitude": [r[10] if len(r) > 10 else np.nan for r in rows],
    })
    return _coerce(frame)


def download_fdsn(datacenter, starttime, endtime, minlatitude, maxlatitude,
                  minlongitude, maxlongitude, minmagnitude=None,
                  maxmagnitude=None, mindepth=None, maxdepth=None,
                  chunk_years=1.0, timeout=300, verbose=True):
    """
    Download an event catalog from an FDSN event service.

    FDSN servers cap a single response (10k rows at NCEDC) and truncate
    silently rather than erroring, so any window that comes back full is split
    in half and re-requested until every piece is under the cap. This is why
    the function takes a whole time span rather than making the caller chunk.

    Args:
        datacenter (str): key of FDSN_SERVICES, or a full query URL
        starttime, endtime (str | datetime): time span, UTC
        minlatitude, maxlatitude (float): degrees north
        minlongitude, maxlongitude (float): degrees east
        minmagnitude, maxmagnitude (float): optional magnitude bounds
        mindepth, maxdepth (float): optional depth bounds, km
        timeout (int): per-request timeout, seconds
        verbose (bool): report progress per chunk
    Returns:
        pd.DataFrame: canonical catalog, sorted by time ascending
    """
    url = FDSN_SERVICES.get(str(datacenter).upper(), datacenter)
    start, end = pd.Timestamp(starttime), pd.Timestamp(endtime)
    if start >= end:
        raise ValueError(f"starttime {start} is not before endtime {end}")

    box = {
        "minlatitude": minlatitude, "maxlatitude": maxlatitude,
        "minlongitude": minlongitude, "maxlongitude": maxlongitude,
        "format": "text", "orderby": "time",
    }
    for key, value in [("minmagnitude", minmagnitude), ("maxmagnitude", maxmagnitude),
                       ("mindepth", mindepth), ("maxdepth", maxdepth)]:
        if value is not None:
            box[key] = value

    collected = []

    def fetch(win_start, win_end, depth=0):
        params = dict(box)
        params["starttime"] = win_start.strftime("%Y-%m-%dT%H:%M:%S")
        params["endtime"] = win_end.strftime("%Y-%m-%dT%H:%M:%S")
        chunk = _parse_fdsn_text(_fdsn_request(url, params, timeout))

        if len(chunk) >= SERVER_ROW_CAP:
            middle = win_start + (win_end - win_start) / 2
            # A window that cannot be split further would loop forever; keep
            # the truncated result and warn rather than recursing to nothing.
            if depth > 20 or middle <= win_start or middle >= win_end:
                print(f"  ! {win_start.date()}..{win_end.date()} still at the "
                      f"{SERVER_ROW_CAP} row cap and cannot be split further; "
                      f"this window is incomplete")
                collected.append(chunk)
                return
            fetch(win_start, middle, depth + 1)
            fetch(middle, win_end, depth + 1)
            return

        if verbose and len(chunk):
            print(f"  {win_start.date()} .. {win_end.date()}  {len(chunk):6d} events")
        collected.append(chunk)

    if verbose:
        print(f"Downloading from {url}")
    # Walk fixed windows rather than recursing from the whole span: a decades
    # -long first request just hits the cap and is thrown away, and each retry
    # costs a full round-trip. Windows only subdivide when they come back full.
    step = pd.Timedelta(days=365.25 * chunk_years)
    window_start = start
    while window_start < end:
        window_end = min(window_start + step, end)
        fetch(window_start, window_end)
        window_start = window_end

    if not collected:
        return pd.DataFrame(columns=COLUMNS)
    catalog = pd.concat(collected, ignore_index=True)
    return _finalise(catalog)


# ---------------------------------------------------------------------------
# read existing files
# ---------------------------------------------------------------------------

def read_file(path, fmt="auto", **kwargs):
    """
    Load a catalog already on disk.

    Args:
        path (str): catalog file
        fmt (str): 'auto', 'fdsn_text', 'scedc', 'growclust', 'hys_mat', 'csv'
        **kwargs: passed to the chosen reader (e.g. colmap= for 'csv')
    Returns:
        pd.DataFrame: canonical catalog, sorted by time ascending
    """
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    if fmt == "auto":
        fmt = detect_format(path)
    readers = {
        "fdsn_text": _read_fdsn_text_file,
        "scedc": _read_scedc,
        "growclust": _read_growclust,
        "hys_mat": _read_hys_mat,
        "csv": _read_csv,
    }
    if fmt not in readers:
        raise ValueError(f"unknown format {fmt!r}; choose from {sorted(readers)}")
    return _finalise(readers[fmt](path, **kwargs))


def detect_format(path):
    """Guess a catalog's format from its extension and first few lines."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".mat":
        return "hys_mat"
    if ext == ".gc":
        return "growclust"

    with open(path, "r", errors="replace") as handle:
        head = [handle.readline() for _ in range(6)]
    joined = "".join(head)

    if "#EventID" in joined or joined.count("|") > 10:
        return "fdsn_text"
    if "SCEDC" in joined or "#YYY/MM/DD" in joined:
        return "scedc"
    if ext in (".csv", ".txt") and "," in joined:
        return "csv"
    return "csv"


def _read_fdsn_text_file(path):
    """A catalog saved straight from an FDSN text query."""
    with open(path, "r", errors="replace") as handle:
        return _parse_fdsn_text(handle.read())


def _read_scedc(path):
    """SCEDC 'date_mag_loc' text catalog (the format behind data/1977_*.txt)."""
    names = ["date", "time", "ET", "GT", "magnitude", "M", "latitude",
             "longitude", "depth", "Q", "EVID", "NPH", "NGRM"]
    raw = pd.read_csv(path, sep=r"\s+", header=None, names=names,
                      comment="#", index_col=False, skiprows=3,
                      on_bad_lines="skip")
    # The file ends with a literal '</PRE>' tag that parses as a junk row.
    raw = raw[pd.to_numeric(raw["latitude"], errors="coerce").notna()]
    frame = pd.DataFrame({
        "time": _combine_date_time(raw["date"], raw["time"]),
        "latitude": raw["latitude"], "longitude": raw["longitude"],
        "depth": raw["depth"], "magnitude": raw["magnitude"],
        "magtype": raw["M"].astype(str), "evid": raw["EVID"],
    })
    return _coerce(frame)


def _read_growclust(path):
    """GrowClust relocated catalog (.gc)."""
    names = ["yr", "mon", "day", "hr", "min", "sec", "eID", "latR", "lonR",
             "depR", "mag", "qID", "cID", "nbranch", "qnpair", "qndiffP",
             "qndiffS", "rmsP", "rmsS", "eh", "ez", "et", "latC", "lonC",
             "depC", "eType", "magType", "GrowClust", "relocBox"]
    raw = pd.read_csv(path, sep=r"\s+", header=None, names=names,
                      index_col=False, on_bad_lines="skip")
    raw = raw[pd.to_numeric(raw["latR"], errors="coerce").notna()]
    dates = (raw["yr"].astype(int).astype(str).str.zfill(4) + "-"
             + raw["mon"].astype(int).astype(str).str.zfill(2) + "-"
             + raw["day"].astype(int).astype(str).str.zfill(2))
    times = (raw["hr"].astype(int).astype(str).str.zfill(2) + ":"
             + raw["min"].astype(int).astype(str).str.zfill(2) + ":"
             + raw["sec"].astype(float).map("{:06.3f}".format))
    frame = pd.DataFrame({
        "time": _combine_date_time(dates, times),
        "latitude": raw["latR"], "longitude": raw["lonR"],
        "depth": raw["depR"], "magnitude": raw["mag"],
        "magtype": raw.get("magType", ""), "evid": raw["eID"],
    })
    return _coerce(frame)


def _read_hys_mat(path):
    """Hauksson, Yang & Shearer relocated southern California catalog (.mat)."""
    from scipy.io import loadmat

    mat = loadmat(path)
    get = lambda key: np.asarray(mat[key]).ravel()
    decimal_year = get("Time")
    frame = pd.DataFrame({
        "time": _decimal_year_to_datetime(decimal_year),
        "latitude": get("Lat"), "longitude": get("Lon"),
        "depth": get("Depth"), "magnitude": get("Mag"),
        "magtype": "", "evid": np.arange(decimal_year.size).astype(str),
    })
    return _coerce(frame)


def _read_csv(path, colmap=None):
    """
    Generic CSV. Column names are matched case-insensitively against common
    spellings; pass colmap={'canonical': 'your_column'} to override.
    """
    raw = pd.read_csv(path)
    aliases = {
        "time": ["time", "datetime", "origintime", "origin_time", "date"],
        "latitude": ["latitude", "lat"],
        "longitude": ["longitude", "lon", "long"],
        "depth": ["depth", "depth_km", "dep", "z"],
        "magnitude": ["magnitude", "mag", "m"],
        "magtype": ["magtype", "mag_type", "magnitudetype"],
        "evid": ["evid", "eventid", "event_id", "id"],
    }
    lowered = {str(c).lower().replace(" ", "_"): c for c in raw.columns}
    frame = pd.DataFrame()
    for canonical, options in aliases.items():
        source = (colmap or {}).get(canonical)
        if source is None:
            source = next((lowered[o] for o in options if o in lowered), None)
        if source is not None:
            frame[canonical] = raw[source]
        elif canonical in ("magtype", "evid"):
            frame[canonical] = ""
        else:
            raise ValueError(
                f"could not find a '{canonical}' column in {os.path.basename(path)}; "
                f"columns are {list(raw.columns)}. Pass colmap={{'{canonical}': ...}}."
            )
    return _coerce(frame)


# ---------------------------------------------------------------------------
# shared helpers
# ---------------------------------------------------------------------------

def _combine_date_time(dates, times):
    """
    Join separate date and time-of-day columns into timestamps.

    Network catalogs round origin seconds for display, which produces values
    like '05:47:60.00' (really 59.995) and even '23:24:60.99'. Those are valid
    events but pandas cannot parse them, so seconds are added as an offset to
    the HH:MM baseline instead of being parsed in place - 60 then simply rolls
    over into the next minute.

    Args:
        dates (pd.Series): date strings, any pandas-parseable form
        times (pd.Series): 'HH:MM:SS.ss' strings
    Returns:
        pd.Series: datetime64
    """
    parts = times.astype(str).str.strip().str.split(":", expand=True)
    if parts.shape[1] < 3:
        return pd.to_datetime(dates.astype(str).str.replace("/", "-", regex=False)
                              + " " + times.astype(str), errors="coerce")
    hours = pd.to_numeric(parts[0], errors="coerce")
    minutes = pd.to_numeric(parts[1], errors="coerce")
    seconds = pd.to_numeric(parts[2], errors="coerce")
    midnight = pd.to_datetime(dates.astype(str).str.replace("/", "-", regex=False),
                              errors="coerce")
    return (midnight
            + pd.to_timedelta(hours, unit="h")
            + pd.to_timedelta(minutes, unit="m")
            + pd.to_timedelta(seconds, unit="s"))


def _decimal_year_to_datetime(values):
    """Convert decimal years (1981.5 -> mid-1981) to datetimes."""
    values = np.asarray(values, dtype=float)
    years = np.floor(values).astype(int)
    starts = pd.to_datetime([f"{y}-01-01" for y in years])
    lengths = np.array([366 if (y % 4 == 0 and y % 100 != 0) or y % 400 == 0 else 365
                        for y in years], dtype=float)
    offsets = pd.to_timedelta((values - years) * lengths, unit="D")
    return starts + offsets


def _coerce(frame):
    """Force canonical dtypes on a partially built frame."""
    out = pd.DataFrame(index=frame.index)
    out["time"] = pd.to_datetime(frame["time"], errors="coerce", format="mixed", utc=True)
    out["time"] = out["time"].dt.tz_localize(None)
    for column in ("latitude", "longitude", "depth", "magnitude"):
        out[column] = pd.to_numeric(frame[column], errors="coerce")
    out["magtype"] = frame.get("magtype", "").astype(str).replace("nan", "")
    out["evid"] = frame.get("evid", "").astype(str)
    return out[COLUMNS]


def _finalise(frame):
    """Drop unusable rows, de-duplicate on evid, sort by time."""
    before = len(frame)
    unusable = frame[["time", "latitude", "longitude"]].isna().any(axis=1)
    # Losing events to a parsing quirk biases every b-value computed downstream,
    # so a silent drop is never acceptable - say what went and why.
    if unusable.any():
        reasons = ", ".join(
            f"{col}={frame[col].isna().sum()}"
            for col in ("time", "latitude", "longitude") if frame[col].isna().any()
        )
        print(f"  ! dropped {unusable.sum()} of {before} events with unusable "
              f"fields ({reasons})")
    frame = frame[~unusable]
    # Adaptive chunking splits on a shared boundary instant, so the same event
    # can arrive twice. Only rows carrying a real id can be de-duplicated;
    # rows without one are kept as-is rather than collapsed together.
    has_id = frame["evid"].ne("") & frame["evid"].ne("nan")
    frame = pd.concat([
        frame[has_id].drop_duplicates(subset="evid", keep="first"),
        frame[~has_id],
    ])
    return frame.sort_values("time").reset_index(drop=True)


def summarise(frame, label=""):
    """One-screen description of a catalog. Returns the text as well as printing it."""
    if frame.empty:
        text = f"{label}: empty catalog"
        print(text)
        return text
    magnitudes = frame["magnitude"].dropna()
    lines = [
        f"{label or 'catalog'}: {len(frame):,} events",
        f"  time      {frame['time'].min()}  ..  {frame['time'].max()}",
        f"  latitude  {frame['latitude'].min():.4f} .. {frame['latitude'].max():.4f}",
        f"  longitude {frame['longitude'].min():.4f} .. {frame['longitude'].max():.4f}",
        f"  depth km  {frame['depth'].min():.2f} .. {frame['depth'].max():.2f}"
        f"   ({frame['depth'].isna().sum():,} missing)",
        f"  magnitude {magnitudes.min():.2f} .. {magnitudes.max():.2f}"
        f"   ({frame['magnitude'].isna().sum():,} missing)",
    ]
    types = [t for t in frame["magtype"].value_counts().head(4).items() if t[0]]
    if types:
        lines.append("  mag types " + ", ".join(f"{k}={v:,}" for k, v in types))
    text = "\n".join(lines)
    print(text)
    return text


def write(frame, path, provenance=None):
    """
    Save a canonical catalog as CSV, alongside a .json sidecar recording where
    it came from. Without the sidecar a catalog on disk is unreproducible.

    Args:
        frame (pd.DataFrame): canonical catalog
        path (str): destination .csv
        provenance (dict): query parameters / source description
    Returns:
        str: the path written
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    frame.to_csv(path, index=False)

    record = dict(provenance or {})
    record.update({
        "retrieved_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "n_events": int(len(frame)),
        "catalog_file": os.path.basename(path),
    })
    if not frame.empty:
        record["time_range"] = [str(frame["time"].min()), str(frame["time"].max())]
    sidecar = os.path.splitext(path)[0] + ".json"
    with open(sidecar, "w") as handle:
        json.dump(record, handle, indent=2)
    return path
