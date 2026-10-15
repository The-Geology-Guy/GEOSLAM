"""Mapped fault traces for comparison with seismo-lineaments.

* USGS Quaternary Fault and Fold Database (US) - ArcGIS REST query by bounding
  box, layer 21 (national).  Replaces the 16 MB ``Qfaults_GIS.zip`` download the
  2020 code fetched on every run (the static zip is still published at
  https://earthquake.usgs.gov/static/lfs/nshm/qfaults/Qfaults_GIS.zip).
* GEM Global Active Faults Database (global; Styron & Pagani 2020, CC-BY-SA 4.0).

Faults are returned as plain dicts: ``{"properties": {...}, "lines": [ndarray(N,2) lon/lat, ...]}``.
"""
from __future__ import annotations

import json

import numpy as np

from .base import cache_dir, http_get

QFAULTS = "https://earthquake.usgs.gov/arcgis/rest/services/haz/Qfaults/MapServer/{layer}/query"
GEM_URL = ("https://raw.githubusercontent.com/GEMScienceTools/gem-global-active-faults/master/"
           "geojson/gem_active_faults_harmonized.geojson")


def _lines(geom):
    if geom is None:
        return []
    if geom["type"] == "LineString":
        return [np.asarray(geom["coordinates"], float)[:, :2]]
    if geom["type"] == "MultiLineString":
        return [np.asarray(c, float)[:, :2] for c in geom["coordinates"]]
    return []


def usgs_qfaults(west, south, east, north, layer=21, page=2000) -> list[dict]:
    out, offset = [], 0
    while True:
        js = http_get(QFAULTS.format(layer=layer), dict(
            geometry=f"{west},{south},{east},{north}", geometryType="esriGeometryEnvelope",
            inSR=4326, spatialRel="esriSpatialRelIntersects", outFields="*", outSR=4326,
            f="geojson", resultOffset=offset, resultRecordCount=page)).json()
        feats = js.get("features", [])
        out += [{"properties": f.get("properties", {}), "lines": _lines(f.get("geometry"))} for f in feats]
        if len(feats) < page and not js.get("exceededTransferLimit"):
            return out
        offset += len(feats)


def gem_faults(west, south, east, north) -> list[dict]:
    path = cache_dir() / "gem_active_faults_harmonized.geojson"
    if not path.exists():
        path.write_bytes(http_get(GEM_URL, timeout=600).content)
    data = json.loads(path.read_text())
    out = []
    for f in data["features"]:
        ls = _lines(f.get("geometry"))
        if any((l[:, 0].max() >= west) and (l[:, 0].min() <= east) and (l[:, 1].max() >= south)
               and (l[:, 1].min() <= north) for l in ls):
            out.append({"properties": f.get("properties", {}), "lines": ls})
    return out


def faults_in_box(west, south, east, north, source="auto"):
    if source in ("auto", "usgs"):
        try:
            res = usgs_qfaults(west, south, east, north)
            if res or source == "usgs":
                return res
        except Exception:
            if source == "usgs":
                raise
    return gem_faults(west, south, east, north)
