"""Parsers against recorded service responses, plus opt-in live tests.

Live tests run only when GEOSLAM_NETWORK_TESTS=1 (GitHub Actions sets it).
"""
import json
import os
from pathlib import Path

import numpy as np
import pytest

from geoslam.orientation import plane_normal
from geoslam.sources.comcat import origin_uncertainty_from_geojson, parse_detail
from geoslam.sources.quakeml import parse_quakeml

FIX = Path(__file__).parent / "fixtures"
NETWORK = os.environ.get("GEOSLAM_NETWORK_TESTS") == "1"


def _detail():
    return json.loads((FIX / "comcat_nc72282711_detail.json").read_text())


def _fake_fetch(url):
    assert "focal-mechanism/nc72282711_fm1" in url or "moment-tensor" in url
    if "focal-mechanism" in url:
        return (FIX / "comcat_nc72282711_fm1.quakeml").read_bytes()
    raise OSError("not recorded")


def test_comcat_origin_ellipsoid_units_are_metres():
    props = _detail()["properties"]["products"]["origin"][0]["properties"]
    u = origin_uncertainty_from_geojson(props)
    assert np.allclose(sorted(u.semi_axes), [168, 264, 384])
    # major axis plunges 86 deg: almost all of the 384 m is vertical
    assert u.axes[0][2] < -0.99
    assert u.half_width(plane_normal(90, 0)) < 300   # horizontal half-width << 384
    us = origin_uncertainty_from_geojson(_detail()["properties"]["products"]["origin"][1]["properties"])
    assert np.allclose(sorted(us.semi_axes), [1600, 1600, 3400])  # km -> m


def test_comcat_detail_parse_with_quakeml_uncertainties():
    eq = parse_detail(_detail(), fetch=_fake_fetch)
    assert eq.event_id == "nc72282711" and abs(eq.depth_km - 11.12) < 1e-9
    fm = eq.mechanisms[0]
    assert fm.product == "focal-mechanism" and fm.method == "FPFIT"
    assert (fm.np1.strike, fm.np1.dip, fm.np1.rake) == (345.0, 85.0, -170.0)
    assert (fm.np1.strike_uncertainty, fm.np1.dip_uncertainty, fm.np1.rake_uncertainty) == (5, 10, 25)
    assert fm.np2.strike == 254.0 and fm.np2.dip_uncertainty == 10
    # moment tensors without recorded QuakeML fall back to GeoJSON planes
    mww = [m for m in eq.mechanisms if m.extra.get("code") == "us_b000s5tp_mww"][0]
    assert (mww.np1.strike, mww.np1.dip, mww.np1.rake) == (64, 80, -13)
    # the nc TMTS product has no nodal planes in GeoJSON and none recorded -> skipped
    assert not any(m.extra.get("code") == "nc72282711_mt1" for m in eq.mechanisms)


def test_quakeml_origin_and_mechanism():
    evs = parse_quakeml((FIX / "comcat_nc72282711_fm1.quakeml").read_bytes())
    assert len(evs) == 1
    e = evs[0]
    assert abs(e.lat - 38.2151667) < 1e-9 and abs(e.depth_km - 11.12) < 1e-9
    assert e.uncertainty.confidence == 0.95
    assert e.mechanisms[0].np1.dip_uncertainty == 10.0


def test_quakeml_hypoinverse_ellipse():
    xml = b"""<q:quakeml xmlns="http://quakeml.org/xmlns/bed/1.2" xmlns:q="http://quakeml.org/xmlns/quakeml/1.2">
<eventParameters><event publicID="smi:nc/Event/1"><origin publicID="o1">
<originUncertainty><minHorizontalUncertainty>300</minHorizontalUncertainty>
<maxHorizontalUncertainty>540</maxHorizontalUncertainty><azimuthMaxHorizontalUncertainty>64</azimuthMaxHorizontalUncertainty>
</originUncertainty><time><value>1992-08-30T23:42:08Z</value></time>
<latitude><value>39.42</value></latitude><longitude><value>-120.19</value></longitude>
<depth><value>5344</value><uncertainty>540</uncertainty></depth></origin>
<focalMechanism publicID="f1"><nodalPlanes><nodalPlane1><strike><value>211</value><uncertainty>15</uncertainty></strike>
<dip><value>61</value><uncertainty>25</uncertainty></dip><rake><value>12</value></rake></nodalPlane1></nodalPlanes>
</focalMechanism></event></eventParameters></q:quakeml>"""
    e = parse_quakeml(xml)[0]
    assert np.allclose(e.uncertainty.semi_axes, [540, 300, 540])
    assert abs(np.degrees(np.arctan2(*e.uncertainty.axes[0][:2])) - 64) < 1e-9
    m = e.mechanisms[0]
    assert m.np2 is not None and abs(m.np2.strike - 115) < 2   # derived auxiliary plane


@pytest.mark.network
@pytest.mark.skipif(not NETWORK, reason="set GEOSLAM_NETWORK_TESTS=1")
def test_live_comcat_napa():
    from geoslam.sources import comcat
    eq = comcat.get_event("nc72282711")
    fm = [m for m in eq.mechanisms if m.product == "focal-mechanism"][0]
    assert fm.np1.dip_uncertainty > 0
    assert eq.uncertainty is not None


@pytest.mark.network
@pytest.mark.skipif(not NETWORK, reason="set GEOSLAM_NETWORK_TESTS=1")
def test_live_ncedc_mechanisms():
    from geoslam.sources import fdsn
    evs = fdsn.query("ncedc", start="2014-08-24T10:20:00", end="2014-08-24T10:21:30",
                     minmag=5.5, minlat=38, maxlat=38.5, minlon=-122.6, maxlon=-122.0)
    assert evs and evs[0].mechanisms


@pytest.mark.network
@pytest.mark.skipif(not NETWORK, reason="set GEOSLAM_NETWORK_TESTS=1")
def test_live_3dep_and_qfaults():
    from geoslam.sources.dems import get_dem
    from geoslam.sources.faults import usgs_qfaults
    dem = get_dem(-122.3123, 38.2152, 3000, source="3dep", cellsize=30)
    assert dem.zone.epsg == 32610 and np.isfinite(dem.z).mean() > 0.9
    assert 0 < np.nanmedian(dem.z) < 500
    f = usgs_qfaults(-122.5, 38.0, -122.1, 38.4)
    assert len(f) > 0
