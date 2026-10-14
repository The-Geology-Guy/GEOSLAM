"""TexNet (Texas Seismological Network, BEG / UT Austin) earthquake catalog.

ArcGIS REST MapServer (the FeatureServer returns HTTP 500):
https://maps.texnet.beg.utexas.edu/arcgis/rest/services/catalog/TexNetCatalog/MapServer
Layer 0 = reviewed earthquakes.  Max 2000 records per page.

TexNet provides location errors (``LatitudeError``, ``LongitudeError``,
``DepthUncertainty``, km) and focal mechanisms (``Strike1/Dip1/Rake1``,
``Strike2/Dip2/Rake2``) with quality metrics (azimuthal gap, misfit) but **no
per-plane angular uncertainty**; supply one (e.g. Rasaka 2016 used 10 deg).
"""
from __future__ import annotations

from ..orientation import NodalPlane
from ..uncertainty import HypocenterUncertainty
from .base import Earthquake, FocalMechanism, http_get

URL = "https://maps.texnet.beg.utexas.edu/arcgis/rest/services/catalog/TexNetCatalog/MapServer/0/query"


def _get(a, *names):
    for n in names:
        for k, v in a.items():
            if k.lower() == n.lower() and v not in (None, ""):
                return v
    return None


def query(minlon, minlat, maxlon, maxlat, where="1=1", with_mechanism=True,
          default_plane_uncertainty=(10.0, 10.0, 10.0), page=2000) -> list[Earthquake]:
    """Events in a lon/lat box.  ``where`` is an ArcGIS SQL filter, e.g.
    ``"Magnitude >= 3"``.  Plane uncertainties default to 10 deg (strike, dip, rake)."""
    if with_mechanism:
        where = f"({where}) AND Strike1 IS NOT NULL"
    out, offset = [], 0
    while True:
        r = http_get(URL, dict(where=where, outFields="*", f="geojson", resultOffset=offset,
                               resultRecordCount=page, geometry=f"{minlon},{minlat},{maxlon},{maxlat}",
                               geometryType="esriGeometryEnvelope", inSR=4326, outSR=4326,
                               spatialRel="esriSpatialRelIntersects"))
        js = r.json()
        feats = js.get("features", [])
        for f in feats:
            out.append(_parse(f, default_plane_uncertainty))
        if len(feats) < page and not js.get("exceededTransferLimit"):
            break
        offset += len(feats)
    return out


def _parse(f, dpu):
    a = f["properties"]
    lon, lat = f["geometry"]["coordinates"][:2]
    depth = float(_get(a, "Depth", "DepthKm", "Depth_km") or 0.0)
    elat, elon = _get(a, "LatitudeError"), _get(a, "LongitudeError")
    edep = _get(a, "DepthUncertainty", "DepthError")
    unc = None
    if elat is not None and elon is not None and edep is not None:
        elat, elon = float(elat) * 1000, float(elon) * 1000
        # independent N-S / E-W errors -> axis-aligned horizontal ellipse
        if elat >= elon:
            unc = HypocenterUncertainty.from_cronin(elat, 0.0, elon, float(edep) * 1000, source="texnet")
        else:
            unc = HypocenterUncertainty.from_cronin(elon, 90.0, elat, float(edep) * 1000, source="texnet")
    eq = Earthquake(event_id=str(_get(a, "EventId", "EventID", "Event_ID") or ""),
                    time=str(_get(a, "Event_Date", "OriginTime", "Origin_Time") or ""),
                    lon=float(lon), lat=float(lat), depth_km=depth,
                    magnitude=_get(a, "Magnitude", "LocalMagnitude"),
                    magnitude_type=str(_get(a, "MagnitudeType", "Magnitude_Type") or ""),
                    uncertainty=unc, catalog="texnet")
    s1 = _get(a, "Strike1")
    if s1 is not None:
        su, du, ru = dpu
        p1 = NodalPlane(float(s1), float(_get(a, "Dip1")), float(_get(a, "Rake1")), su, du, ru, "NP1")
        if _get(a, "Strike2") is not None:
            p2 = NodalPlane(float(_get(a, "Strike2")), float(_get(a, "Dip2")), float(_get(a, "Rake2")),
                            su, du, ru, "NP2")
        else:
            p2 = p1.auxiliary()
            p2.strike_uncertainty, p2.dip_uncertainty, p2.rake_uncertainty = su, du, ru
        eq.mechanisms.append(FocalMechanism(p1, p2, source="tx", product="texnet",
                                            method=str(_get(a, "FocalMechanismMethodID") or ""),
                                            extra=dict(gap=_get(a, "FocalMechanismAzimuthalGap"),
                                                       misfit=_get(a, "FocalMechanismMisfit"))))
    return eq
