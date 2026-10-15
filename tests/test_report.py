"""Report figures (Cronin-style maps) and terrain / DEM downloads."""
import json
from pathlib import Path

import numpy as np

from geoslam import web

FIX = Path(__file__).parent / "fixtures"
QML_URL = "https://earthquake.usgs.gov/product/focal-mechanism/nc72282711_fm1/nc/1486207774840/quakeml.xml"


def _case():
    detail = (FIX / "comcat_nc72282711_detail.json").read_text()
    return json.loads(web.case_from_comcat(
        detail, json.dumps({QML_URL: (FIX / "comcat_nc72282711_fm1.quakeml").read_text()})))


def _load_dem():
    d = np.load(FIX / "napa_3dep_150m_utm10n.npz")
    web.set_dem_utm(d["z"].astype(np.float32).ravel(), 200, 160, float(d["x_min"]),
                    float(d["y_max"]), 150.0, 32610)


def _items():
    a = _case()
    b = json.loads(json.dumps(a))
    b["event"].update(id="fake", lon=-122.22, lat=38.29, depth_km=6.0)
    m = a["mechanisms"][0]
    return [dict(label="A1", event_label="A", case=a, plane=m["np1"]),
            dict(label="A2", event_label="A", case=a, plane=m["np2"]),
            dict(label="B1", event_label="B", case=b, plane=m["np1"], method="cronin")]


def test_report_plan_covers_all_events():
    a, b = _items()[0]["case"], _items()[2]["case"]
    p = json.loads(web.report_plan(json.dumps([a, b])))
    lo0, la0, lo1, la1 = p["lonlat_bounds"]
    for c in (a, b):
        assert lo0 < c["event"]["lon"] < lo1 and la0 < c["event"]["lat"] < la1
    assert max(p["size"]) <= 1200 and p["epsg"] == 32610


def test_report_png_pdf_svg():
    _load_dem()
    faults = [[[-122.31, 38.25], [-122.33, 38.30]]]
    quakes = [[-122.30, 38.22, 3.1], [-122.32, 38.24, 2.0]]
    png = web.report_render(json.dumps(_items()), json.dumps(
        dict(dpi=80, legend=True, details=True, title="Test", faults=faults, quakes=quakes)))
    assert png[:8] == b"\x89PNG\r\n\x1a\n" and len(png) > 20_000
    pdf = web.report_render(json.dumps(_items()[:1]), json.dumps(dict(fmt="pdf", page="letter-portrait")))
    assert pdf[:5] == b"%PDF-"
    svg = web.report_render(json.dumps(_items()[:1]), json.dumps(dict(fmt="svg", terrain="slope")))
    assert b"<svg" in svg[:500]


def test_terrain_and_dem_downloads():
    _load_dem()
    for style in ("hillshade", "hillshade-color", "elevation", "slope"):
        assert web.terrain_png(style, 315, 30)[:4] == b"\x89PNG"
    asc = web.dem_ascii()
    assert asc.startswith("ncols 160\nnrows 200\n")
    assert len(web.world_file().split()) == 6
    assert 'UTM_Zone_10N' in web.projection_wkt() and "-123" in web.projection_wkt()
    assert web.save_dem("main") and json.loads(web.load_dem("main"))["shape"] == [200, 160]


def test_terrain_zip_and_light():
    import io
    import zipfile
    _load_dem()
    az = web.report_light_azimuth(json.dumps(_items()[0]))
    assert 0 <= az < 360
    z = zipfile.ZipFile(io.BytesIO(web.terrain_zip("t", "hillshade", az, 15)))
    assert sorted(z.namelist()) == ["README.txt", "t.pgw", "t.png", "t.prj"]
    z = zipfile.ZipFile(io.BytesIO(web.dem_zip("d")))
    assert "d.asc" in z.namelist()
