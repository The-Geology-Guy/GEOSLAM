"""The browser bridge (geoslam.web) with recorded South Napa data."""
import json
from pathlib import Path

import numpy as np

from geoslam import web

FIX = Path(__file__).parent / "fixtures"
QML_URL = "https://earthquake.usgs.gov/product/focal-mechanism/nc72282711_fm1/nc/1486207774840/quakeml.xml"


def _case():
    detail = (FIX / "comcat_nc72282711_detail.json").read_text()
    return web.case_from_comcat(detail, json.dumps({QML_URL: (FIX / "comcat_nc72282711_fm1.quakeml").read_text()}))


def test_case_and_plan():
    c = json.loads(_case())
    assert c["mechanisms"][0]["label"] == "NC FPFIT" and c["mechanisms"][0]["has_uncertainty"]
    plan = json.loads(web.dem_plan(_case()))
    assert plan["epsg"] == 32610 and plan["size"][0] <= 800
    assert 14_000 <= (plan["bbox"][2] - plan["bbox"][0]) / 2 <= 40_000


def test_run_returns_polygons_and_traces():
    d = np.load(FIX / "napa_3dep_150m_utm10n.npz")
    info = json.loads(web.set_dem_utm(d["z"].astype(np.float32).ravel(), 200, 160, float(d["x_min"]),
                                      float(d["y_max"]), 150.0, 32610))
    assert info["valid_fraction"] == 1.0
    c = json.loads(_case())
    m = c["mechanisms"][0]
    opts = {"method": "envelope", "planes": [{"label": "NP1", "plane": m["np1"]}, {"label": "NP2", "plane": m["np2"]}]}
    out = json.loads(web.run(_case(), json.dumps(opts)))
    p1, p2 = out["planes"]
    assert p1["polygons"] and p1["traces"]["Mean"] and abs(p1["area_km2"] - 178.8) < 2
    lon, lat = np.array(p1["polygons"][0][0]).T
    assert -122.45 < lon.mean() < -122.15 and 38.05 < lat.mean() < 38.35
    rgba = web.hillshade(p1["light_azimuth"], 15)[0]
    assert len(rgba) == 200 * 160 * 4
