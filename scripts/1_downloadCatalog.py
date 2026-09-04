"""
Step 1 - acquire an earthquake catalog.

Two paths, one output. Either download fresh from an FDSN event service, or
normalise a catalog you already have; both write the same canonical CSV plus a
JSON provenance sidecar into output/1_catalog/, so every later step reads one
format and never needs to know which path was taken.

Examples
--------
# fresh download, named region
python scripts/1_downloadCatalog.py --region parkfield

# fresh download, explicit box
python scripts/1_downloadCatalog.py --datacenter NCEDC \
    --bbox 35.6 36.3 -120.8 -120.0 --start 1981-01-01 --end 2026-09-02 \
    --name parkfield

# use a catalog you already have (format auto-detected)
python scripts/1_downloadCatalog.py --file data/1977_20240611.txt --name socal_scedc

# list what is available
python scripts/1_downloadCatalog.py --list-regions
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src import catalog  # noqa: E402

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_DIR = os.path.join(PROJECT_DIR, "output", "1_catalog")

# Named regions are a convenience only - the pipeline itself is area-agnostic
# and any box can be given with --bbox. 'parkfield' reproduces the study area
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


def parse_args():
    parser = argparse.ArgumentParser(
        description="Step 1: download a catalog, or normalise one you already have.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    source = parser.add_argument_group("source (choose one)")
    source.add_argument("--region", choices=sorted(REGIONS),
                        help="named region preset to download")
    source.add_argument("--file", help="path to a catalog already on disk")
    source.add_argument("--list-regions", action="store_true",
                        help="show the presets and exit")

    download = parser.add_argument_group("download options")
    download.add_argument("--datacenter", default=None,
                          help=f"FDSN service: {', '.join(catalog.FDSN_SERVICES)}, "
                               f"or a full query URL")
    download.add_argument("--bbox", nargs=4, type=float, metavar=("MINLAT", "MAXLAT",
                          "MINLON", "MAXLON"), help="bounding box in degrees")
    download.add_argument("--start", help="start time, e.g. 1981-01-01")
    download.add_argument("--end", help="end time (default: now)")
    download.add_argument("--minmag", type=float, help="minimum magnitude")
    download.add_argument("--maxdepth", type=float, help="maximum depth, km")

    other = parser.add_argument_group("file options")
    other.add_argument("--format", default="auto",
                       help="override format detection: fdsn_text, scedc, "
                            "growclust, hys_mat, csv")

    parser.add_argument("--name", help="output basename (default: region or file stem)")
    parser.add_argument("--force", action="store_true",
                        help="overwrite an existing catalog of the same name")
    return parser.parse_args()


def main():
    args = parse_args()

    if args.list_regions:
        print("Named regions:\n")
        for name, cfg in sorted(REGIONS.items()):
            lat0, lat1, lon0, lon1 = cfg["bbox"]
            print(f"  {name}")
            print(f"    datacenter {cfg['datacenter']}")
            print(f"    box        {lat0} .. {lat1} N,  {lon0} .. {lon1} E")
            print(f"    from       {cfg['start']}")
            print(f"    {cfg['note']}\n")
        print("Any other area: --datacenter ... --bbox MINLAT MAXLAT MINLON MAXLON")
        return 0

    if bool(args.region) == bool(args.file):
        print("error: give exactly one of --region / --file "
              "(or --bbox with --datacenter). See --help.", file=sys.stderr)
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
              f"       Re-run with --force to overwrite, or pass a different --name.",
              file=sys.stderr)
        return 1

    # ---- fetch ------------------------------------------------------------
    if args.file:
        print(f"Reading {args.file}")
        events = catalog.read_file(args.file, fmt=args.format)
        provenance["detected_format"] = catalog.detect_format(args.file) \
            if args.format == "auto" else args.format
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
    print(f"\nWrote {os.path.relpath(destination, PROJECT_DIR)}")
    print(f"      {os.path.relpath(os.path.splitext(destination)[0] + '.json', PROJECT_DIR)}")
    return 0


def pd_now():
    import pandas as pd
    return pd.Timestamp.utcnow().tz_localize(None).strftime("%Y-%m-%d")


if __name__ == "__main__":
    sys.exit(main())
