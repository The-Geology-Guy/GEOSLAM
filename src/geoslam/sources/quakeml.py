"""Minimal QuakeML 1.2 reader for what SLAM needs (stdlib only).

Extracts, per event: preferred origin position, its uncertainty (confidence
ellipsoid, or the horizontal ellipse min/max/azimuth + depth uncertainty that
NCEDC/SCEDC report from Hypoinverse), magnitude, and every focal mechanism's
nodal planes *with their uncertainties*.

Units in QuakeML are SI: depth and uncertainties in metres.
"""
from __future__ import annotations

import xml.etree.ElementTree as ET

from ..orientation import NodalPlane
from ..uncertainty import HypocenterUncertainty
from .base import Earthquake, FocalMechanism


def _local(tag):
    return tag.rsplit("}", 1)[-1]


def _child(el, *names):
    for name in names:
        if el is None:
            return None
        el = next((c for c in el if _local(c.tag) == name), None)
    return el


def _children(el, name):
    return [c for c in el if _local(c.tag) == name] if el is not None else []


def _float(e):
    try:
        return float(e.text)
    except (TypeError, ValueError, AttributeError):
        return None


def _num(el, *path):
    """Number at ``path``; if that element is a RealQuantity, its <value>."""
    e = _child(el, *path)
    if e is None:
        return None
    v = _child(e, "value")
    return _float(v if v is not None else e)


def _quantity(el, name):
    """(value, uncertainty) of a RealQuantity child."""
    q = _child(el, name)
    if q is None:
        return None, None
    return _float(_child(q, "value")), _float(_child(q, "uncertainty"))


def _id(el):
    return el.get("publicID", "") if el is not None else ""


def parse_origin_uncertainty(origin) -> HypocenterUncertainty | None:
    """Best available location-uncertainty description of a QuakeML origin."""
    ou = _child(origin, "originUncertainty")
    _, dz = _quantity(origin, "depth")
    conf = None
    if ou is not None:
        cl = _child(ou, "confidenceLevel")
        conf = float(cl.text) / 100.0 if cl is not None and cl.text else None
        ell = _child(ou, "confidenceEllipsoid")
        if ell is not None:
            vals = {k: _num(ell, k) for k in ("semiMajorAxisLength", "semiIntermediateAxisLength",
                                              "semiMinorAxisLength", "majorAxisAzimuth",
                                              "majorAxisPlunge", "majorAxisRotation")}
            if all(vals[k] is not None for k in list(vals)[:5]):
                return HypocenterUncertainty.from_quakeml_ellipsoid(
                    vals["semiMajorAxisLength"], vals["semiIntermediateAxisLength"],
                    vals["semiMinorAxisLength"], vals["majorAxisAzimuth"],
                    vals["majorAxisPlunge"], vals["majorAxisRotation"] or 0.0,
                    confidence=conf, source="quakeml:confidenceEllipsoid")
        mx, mn = _num(ou, "maxHorizontalUncertainty"), _num(ou, "minHorizontalUncertainty")
        az = _num(ou, "azimuthMaxHorizontalUncertainty")
        if mx is not None and mn is not None and dz is not None:
            return HypocenterUncertainty.from_cronin(mx, az or 0.0, mn, dz, conf,
                                                     "quakeml:horizontalEllipse+depth")
        h = _num(ou, "horizontalUncertainty")
        if h is not None and dz is not None:
            return HypocenterUncertainty.isotropic_horizontal(h, dz, conf,
                                                              "quakeml:horizontalUncertainty+depth")
    lat_u = _quantity(origin, "latitude")[1]
    if dz is not None and lat_u is not None:
        h = lat_u * 111_195.0
        return HypocenterUncertainty.isotropic_horizontal(h, dz, conf, "quakeml:lat/depth uncertainty")
    return None


def _plane(np_el, label):
    s, su = _quantity(np_el, "strike")
    d, du = _quantity(np_el, "dip")
    r, ru = _quantity(np_el, "rake")
    if s is None or d is None:
        return None
    return NodalPlane(s, d, r, su or 0.0, du or 0.0, ru or 0.0, label)


def parse_focal_mechanism(fm, source="") -> FocalMechanism | None:
    nps = _child(fm, "nodalPlanes")
    if nps is None:
        return None
    p1 = _plane(_child(nps, "nodalPlane1"), "NP1")
    p2 = _plane(_child(nps, "nodalPlane2"), "NP2")
    if p1 is None:
        return None
    if p2 is None and p1.rake is not None:
        p2 = p1.auxiliary()
        p2.label = "NP2"
    pref = nps.get("preferredPlane")
    method = _child(fm, "methodID")
    extra = {}
    for k in ("azimuthalGap", "stationPolarityCount", "misfit", "stationDistributionRatio"):
        e = _child(fm, k)
        if e is not None and e.text:
            extra[k] = e.text
    q = _child(fm, "creationInfo", "agencyID")
    return FocalMechanism(p1, p2, source=source or (q.text if q is not None else ""),
                          product="quakeml", method=(method.text.rsplit("/", 1)[-1] if method is not None and method.text else ""),
                          preferred_plane=int(pref) if pref else None, extra=extra)


def parse_quakeml(xml_bytes, catalog="") -> list[Earthquake]:
    root = ET.fromstring(xml_bytes)
    ep = next((e for e in root.iter() if _local(e.tag) == "eventParameters"), root)
    out = []
    for ev in _children(ep, "event"):
        origins = {_id(o): o for o in _children(ev, "origin")}
        pref = _child(ev, "preferredOriginID")
        origin = origins.get(pref.text.strip() if pref is not None and pref.text else "",
                             next(iter(origins.values()), None))
        if origin is None:
            continue
        mags = {_id(m): m for m in _children(ev, "magnitude")}
        pm = _child(ev, "preferredMagnitudeID")
        mag = mags.get(pm.text.strip() if pm is not None and pm.text else "", next(iter(mags.values()), None))
        t = _child(origin, "time", "value")
        eid = _id(ev)
        for attr, v in ev.attrib.items():
            if attr.endswith("eventid"):
                eid = v
        eq = Earthquake(
            event_id=eid.rsplit("/", 1)[-1].split("=")[-1],
            time=t.text if t is not None else "",
            lon=_num(origin, "longitude"), lat=_num(origin, "latitude"),
            depth_km=(_num(origin, "depth") or 0.0) / 1000.0,
            magnitude=_num(mag, "mag") if mag is not None else None,
            magnitude_type=(_child(mag, "type").text if mag is not None and _child(mag, "type") is not None else ""),
            uncertainty=parse_origin_uncertainty(origin), catalog=catalog,
        )
        for fm in _children(ev, "focalMechanism"):
            m = parse_focal_mechanism(fm, catalog)
            if m is not None:
                eq.mechanisms.append(m)
        out.append(eq)
    return out
