"""Validation case: 2014-08-24 M6.0 South Napa earthquake (Worrell, 2016 MS thesis).

Uses recorded data in tests/fixtures (USGS 3DEP DEM at 150 m in UTM 10N, USGS
Quaternary faults, ComCat aftershocks) so it runs offline.

Expected (Worrell 2016): the NP1 seismo-lineament (dip direction 075, dip 85,
i.e. strike 345) contains the West Napa fault / Browns Valley strand and the
mapped 2014 ruptures north of the epicentre.

Run:  python examples/napa_2014_validation.py   (writes PNG + GeoJSON to examples/output)
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from geoslam import DEMGrid, Hypocenter, HypocenterUncertainty, NodalPlane, UTMZone, seismo_lineament
from geoslam.geodesy import lonlat_to_utm
from geoslam.plotting import plot_slam
from geoslam.sources.comcat import parse_detail

FIX = Path(__file__).resolve().parents[1] / "tests" / "fixtures"
OUT = Path(__file__).resolve().parent / "output"


def load_dem():
    d = np.load(FIX / "napa_3dep_150m_utm10n.npz")
    return DEMGrid(d["z"].astype(float), float(d["x_min"]), float(d["y_max"]), float(d["cellsize"]),
                   UTMZone(int(d["zone"])), "NAVD88", "USGS 3DEP (150 m)")


def load_faults():
    c = json.loads((FIX / "napa_qfaults_compact.json").read_text())
    out = []
    for idx, flat in c["lines"]:
        xy = np.asarray(flat, float).reshape(-1, 2)
        out.append({"properties": {"fault_name": c["names"][idx]},
                    "lines": [np.column_stack([xy[:, 0] / 1e4 - 122.5, xy[:, 1] / 1e4 + 38])]})
    return out


def load_aftershocks():
    r = np.asarray(json.loads((FIX / "napa_aftershocks_compact.json").read_text())["rows"], float)
    return np.column_stack([r[:, 0] / 1e4 - 122.5, r[:, 1] / 1e4 + 38, r[:, 3] / 10])


def fraction_inside(result, faults, name):
    """Fraction of a named fault's mapped length lying inside the swath."""
    dem = result.dem
    inside = total = 0.0
    for f in faults:
        if f["properties"]["fault_name"] != name:
            continue
        for line in f["lines"]:
            x, y, _ = lonlat_to_utm(line[:, 0], line[:, 1], dem.zone)
            # densify
            t = np.linspace(0, 1, 20)
            for k in range(len(x) - 1):
                px = x[k] + t * (x[k + 1] - x[k])
                py = y[k] + t * (y[k + 1] - y[k])
                seg = np.hypot(x[k + 1] - x[k], y[k + 1] - y[k]) / len(t)
                col = ((px - dem.x_min) / dem.cellsize).astype(int)
                row = ((dem.y_max - py) / dem.cellsize).astype(int)
                ok = (col >= 0) & (col < dem.shape[1]) & (row >= 0) & (row < dem.shape[0])
                total += seg * ok.sum()
                inside += seg * result.swath[row[ok], col[ok]].sum()
    return inside / total if total else float("nan"), total / 1000


def main():
    OUT.mkdir(exist_ok=True)
    dem, faults, aft = load_dem(), load_faults(), load_aftershocks()

    # --- Case A: Worrell (2016) Table 1-2 inputs (double-difference location, NCEDC FPFIT)
    hyp_w = Hypocenter(-122.29787, 38.20666, 8.919,
                       HypocenterUncertainty.from_cronin(28, 90, 28, 29, source="Worrell 2016"),
                       event_id="nc72282711 (Worrell)")
    np1_w = NodalPlane.from_dip_azimuth(75, 85, -170, 5, 10, 25, label="NP1")
    np2_w = NodalPlane(254, 80, -5, 5, 10, 25, label="NP2")

    # --- Case B: everything straight from ComCat (recorded detail + per-product QuakeML)
    detail = json.loads((FIX / "comcat_nc72282711_detail.json").read_text())
    eq = parse_detail(detail, fetch=lambda url: (FIX / "comcat_nc72282711_fm1.quakeml").read_bytes()
                      if "focal-mechanism" in url else (_ for _ in ()).throw(OSError("offline")))
    fm = eq.mechanisms[0]
    hyp_c = eq.hypocenter()

    rows = []
    results = {}
    for tag, hyp, plane in [("worrell_np1", hyp_w, np1_w), ("worrell_np2", hyp_w, np2_w),
                            ("comcat_np1", hyp_c, fm.np1), ("comcat_np2", hyp_c, fm.np2)]:
        for method in ("cronin", "envelope"):
            r = seismo_lineament(dem, hyp, plane, method=method)
            frac, km = fraction_inside(r, faults, "West Napa")
            rows.append((tag, method, r.area_km2, frac, km))
            results[(tag, method)] = r

    lines = ["case          method    swath_km2  WestNapa_inside  (mapped km in DEM)"]
    for tag, m, a, f, km in rows:
        lines.append(f"{tag:<13} {m:<9} {a:9.1f}  {f:15.2f}  ({km:.1f})")
    report = "\n".join(lines)
    print(report)
    print()
    print(results[("worrell_np1", "envelope")].summary())
    print()
    print(results[("comcat_np1", "envelope")].summary())
    (OUT / "napa_2014_report.txt").write_text(report + "\n\n" + results[("worrell_np1", "envelope")].summary()
                                             + "\n\n" + results[("comcat_np1", "envelope")].summary())

    fig, axes = plt.subplots(1, 3, figsize=(24, 10.5), layout="constrained")
    plot_slam(results[("worrell_np1", "cronin")], ax=axes[0], faults=faults, aftershocks=aft,
              title="NP1 strike 345 - Cronin 7-slab union (Worrell inputs)")
    plot_slam(results[("worrell_np1", "envelope")], ax=axes[1], faults=faults, aftershocks=aft,
              title="NP1 strike 345 - full envelope (Worrell inputs)")
    plot_slam(results[("worrell_np2", "envelope")], ax=axes[2], faults=faults, aftershocks=aft,
              title="NP2 strike 254 - full envelope (auxiliary plane)")
    fig.suptitle("2014 M6.0 South Napa: seismo-lineaments vs USGS Quaternary faults (red) "
                 "and 30-day aftershocks (orange)", fontsize=14)

    fig.savefig(OUT / "napa_2014_slam.png", dpi=110)
    (OUT / "napa_2014_np1.geojson").write_text(json.dumps(results[("worrell_np1", "envelope")].to_geojson()))
    return rows


if __name__ == "__main__":
    main()
