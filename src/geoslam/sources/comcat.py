"""USGS ANSS Comprehensive Catalog (ComCat) via the FDSN event web service.

Endpoint: https://earthquake.usgs.gov/fdsnws/event/1/   (no key; 20,000-event cap)

Things that tripped up the 2020 code and that this module handles:

* **Units differ by field.**  In the GeoJSON ``origin`` product,
  ``horizontal-error`` and ``vertical-error`` are in km, but the
  ``error-ellipse-major/-intermediate/-minor`` semi-axes are in **metres**.
  QuakeML is all SI (metres).
* **The error ellipse is 3-D.**  ``error-ellipse-plunge`` is the plunge of the
  major axis (86 deg for the 2014 South Napa mainshock - a near-vertical
  major axis).  Treating it as a horizontal ellipse, as GEOSLAM 0.0.x did,
  overstates the horizontal error and the swath width.
* **Nodal-plane uncertainties are not in the GeoJSON**, nor in the aggregated
  ``format=quakeml`` response.  They are only in each product's own
  ``quakeml.xml``; :func:`get_event` follows those links.
* ``nodal-plane-1-rake`` is the live property name (PDL docs say ``-slip``);
  both are accepted.
* ``products['origin'][0]`` is the preferred origin (highest weight).
"""
from __future__ import annotations

import json

from ..orientation import NodalPlane
from ..uncertainty import HypocenterUncertainty
from .base import Earthquake, FocalMechanism, http_get
from .quakeml import parse_quakeml

BASE = "https://earthquake.usgs.gov/fdsnws/event/1/query"


def search(minlatitude, maxlatitude, minlongitude, maxlongitude, starttime=None, endtime=None,
           minmagnitude=None, maxmagnitude=None, mindepth=None, maxdepth=None,
           producttype: str | None = "focal-mechanism", catalog=None, limit=2000,
           orderby="time", **extra) -> list[dict]:
    """Search ComCat; returns GeoJSON summary features (cheap - no products).

    ``producttype`` restricts to events having that product (``"focal-mechanism"``,
    ``"moment-tensor"`` or ``None`` for all).  Note min/max longitude are the
    western/eastern bounds (the 2020 code passed them swapped).
    """
    params = dict(format="geojson", minlatitude=minlatitude, maxlatitude=maxlatitude,
                  minlongitude=minlongitude, maxlongitude=maxlongitude, starttime=starttime,
                  endtime=endtime, minmagnitude=minmagnitude, maxmagnitude=maxmagnitude,
                  mindepth=mindepth, maxdepth=maxdepth, producttype=producttype, catalog=catalog,
                  limit=limit, orderby=orderby, **extra)
    params = {k: v for k, v in params.items() if v is not None}
    r = http_get(BASE, params)
    if not r.content:
        return []
    return r.json().get("features", [])


def _f(props, key):
    v = props.get(key)
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def origin_uncertainty_from_geojson(props) -> HypocenterUncertainty | None:
    """Location uncertainty from an ``origin`` product's properties."""
    major = _f(props, "error-ellipse-major")
    inter = _f(props, "error-ellipse-intermediate")
    minor = _f(props, "error-ellipse-minor")
    az = _f(props, "error-ellipse-azimuth")
    pl = _f(props, "error-ellipse-plunge")
    rot = _f(props, "error-ellipse-rotation") or 0.0
    if None not in (major, inter, minor, az, pl) and major > 0:
        return HypocenterUncertainty.from_quakeml_ellipsoid(
            major, inter, minor, az, pl, rot, confidence=None,
            source="comcat:error-ellipse (m)")
    herr, verr = _f(props, "horizontal-error"), _f(props, "vertical-error")
    if herr is not None and verr is not None and (herr > 0 or verr > 0):
        return HypocenterUncertainty.isotropic_horizontal(herr * 1000.0, verr * 1000.0,
                                                          source="comcat:horizontal/vertical-error (km)")
    return None


def _planes_from_props(p):
    def g(i, k):
        v = _f(p, f"nodal-plane-{i}-{k}")
        if v is None and k == "rake":
            v = _f(p, f"nodal-plane-{i}-slip")
        return v

    if g(1, "strike") is None:
        return None
    np1 = NodalPlane(g(1, "strike"), g(1, "dip"), g(1, "rake"), label="NP1")
    if g(2, "strike") is not None:
        np2 = NodalPlane(g(2, "strike"), g(2, "dip"), g(2, "rake"), label="NP2")
    elif np1.rake is not None:
        np2 = np1.auxiliary()
        np2.label = "NP2"
    else:
        return None
    return np1, np2


def parse_detail(detail: dict, fetch_uncertainties=True, fetch=None) -> Earthquake:
    """Build an :class:`Earthquake` from a ComCat GeoJSON *detail* feature.

    ``fetch`` (url -> bytes) is used to read per-product QuakeML for nodal-plane
    uncertainties; it defaults to an HTTP GET with caching.
    """
    fetch = fetch or (lambda url: http_get(url, cache=True))
    props = detail["properties"]
    products = props.get("products", {})
    origins = products.get("origin", [])
    lon, lat, depth = detail["geometry"]["coordinates"][:3]
    unc, notes = None, []
    if origins:
        op = origins[0]["properties"]
        lon = _f(op, "longitude") or lon
        lat = _f(op, "latitude") or lat
        depth = _f(op, "depth") if _f(op, "depth") is not None else depth
        unc = origin_uncertainty_from_geojson(op)
        if op.get("depth-type"):
            notes.append(f"depth-type: {op['depth-type']}")
    eq = Earthquake(event_id=detail["id"],
                    time=(origins[0]["properties"].get("eventtime") if origins else str(props.get("time"))),
                    lon=float(lon), lat=float(lat), depth_km=float(depth),
                    magnitude=props.get("mag"), magnitude_type=props.get("magType", ""),
                    uncertainty=unc, catalog="comcat", url=props.get("url", ""), notes=notes)
    for ptype in ("focal-mechanism", "moment-tensor"):
        for prod in products.get(ptype, []):
            planes = _planes_from_props(prod["properties"])
            qurl = prod.get("contents", {}).get("quakeml.xml", {}).get("url", "")
            method = prod["properties"].get("beachball-type", "") or prod["properties"].get("derived-magnitude-type", "")
            fm = None
            if fetch_uncertainties and qurl:
                try:
                    evs = parse_quakeml(fetch(qurl), catalog=prod.get("source", ""))
                    if evs and evs[0].mechanisms:
                        fm = evs[0].mechanisms[0]
                except Exception as exc:  # noqa: BLE001
                    eq.notes.append(f"could not read {qurl}: {exc}")
            if fm is None and planes is None:
                continue
            if fm is None:
                fm = FocalMechanism(*planes)
            fm.source = prod.get("source", "")
            fm.product = ptype
            fm.method = (fm.method or method.rsplit("/", 1)[-1])
            fm.quakeml_url = qurl
            fm.extra.update(code=prod.get("code"), weight=prod.get("preferredWeight"),
                            review=prod["properties"].get("review-status"))
            eq.mechanisms.append(fm)
    return eq


def get_event(event_id: str, fetch_uncertainties=True) -> Earthquake:
    """Fetch one event (all products) and its mechanisms with uncertainties."""
    r = http_get(BASE, dict(eventid=event_id, format="geojson"))
    return parse_detail(r.json(), fetch_uncertainties)


def get_events(features_or_ids, fetch_uncertainties=True, progress=False) -> list[Earthquake]:
    out = []
    items = list(features_or_ids)
    for k, item in enumerate(items):
        eid = item if isinstance(item, str) else item["id"]
        if progress:
            print(f"[{k + 1}/{len(items)}] {eid}")
        out.append(get_event(eid, fetch_uncertainties))
    return out


def load_detail_file(path) -> dict:
    with open(path) as fh:
        return json.load(fh)
