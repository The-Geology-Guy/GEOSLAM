"""Command line: ``geoslam run <comcat-event-id> --out results/``."""
from __future__ import annotations

import argparse


def main(argv=None):
    ap = argparse.ArgumentParser(prog="geoslam", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run SLAM for a ComCat event id")
    r.add_argument("event_id")
    r.add_argument("--out", default="geoslam_out")
    r.add_argument("--method", choices=["envelope", "cronin"], default="envelope")
    r.add_argument("--dem-source", default="auto",
                   choices=["auto", "3dep", "copernicus", "opentopography", "opentopography-usgs"])
    r.add_argument("--cellsize", type=float, default=None)
    r.add_argument("--half-extent-km", type=float, default=None)
    r.add_argument("--plane-uncertainty", type=float, nargs=2, metavar=("STRIKE", "DIP"),
                   help="degrees, used when the catalog gives none")
    r.add_argument("--location-uncertainty", type=float, nargs=2, metavar=("EH_M", "EZ_M"),
                   help="metres, used when the catalog gives none")
    r.add_argument("--depth-datum", choices=["sea_level", "surface"], default="sea_level")
    a = ap.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")
    from .workflow import analyze_event

    res = analyze_event(a.event_id, dem_source=a.dem_source, cellsize=a.cellsize, method=a.method,
                        default_plane_uncertainty=a.plane_uncertainty,
                        unknown_location_uncertainty=a.location_uncertainty,
                        depth_datum=a.depth_datum,
                        half_extent_m=a.half_extent_km * 1000 if a.half_extent_km else None)
    out = res.save(a.out)
    for r_ in res.results:
        print(r_.summary(), "\n")
    for n in res.earthquake.notes:
        print("note:", n)
    print(f"wrote {out}")


if __name__ == "__main__":  # pragma: no cover
    main()
