"""Regional FDSN event services that return mechanisms inside QuakeML.

* NCEDC (northern California):  https://service.ncedc.org/fdsnws/event/1/query
* SCEDC (southern California):  https://service.scedc.caltech.edu/fdsnws/event/1/query

Both accept ``includemechanisms=true`` and report Hypoinverse location errors
as a horizontal ellipse (``min/maxHorizontalUncertainty``,
``azimuthMaxHorizontalUncertainty``) plus depth uncertainty - exactly the
``eh1, eh2, eh1Az, ez`` inputs of Cronin's code.  NCEDC uses the short FDSN
parameter names (``start``, ``minlat`` ...), which are valid everywhere.

For the best hypocentres use the relocated catalogs Cronin's students used:
Waldhauser's NCAeqDD (https://nocaldd.ldeo.columbia.edu) and Hauksson-Yang-
Shearer (https://scedc.caltech.edu/data/alt-2011-dd-hauksson-yang-shearer.html);
see :func:`geoslam.io.read_catalog_table` to bring those in.
"""
from __future__ import annotations

from .base import Earthquake, http_get
from .quakeml import parse_quakeml

SERVICES = {
    "ncedc": "https://service.ncedc.org/fdsnws/event/1/query",
    "scedc": "https://service.scedc.caltech.edu/fdsnws/event/1/query",
    "usgs": "https://earthquake.usgs.gov/fdsnws/event/1/query",
    "isc": "https://www.isc.ac.uk/fdsnws/event/1/query",
    "emsc": "https://www.seismicportal.eu/fdsnws/event/1/query",
    "geofon": "https://geofon.gfz.de/fdsnws/event/1/query",
}

_MECH_PARAM = {"ncedc": "includemechanisms", "scedc": "includemechanisms",
               "geofon": "includefocalmechanism"}


def query(service: str, eventid=None, start=None, end=None, minlat=None, maxlat=None,
          minlon=None, maxlon=None, minmag=None, limit=500, **extra) -> list[Earthquake]:
    """Query an FDSN event service and parse events with their mechanisms."""
    url = SERVICES.get(service, service)
    params = dict(eventid=eventid, start=start, end=end, minlat=minlat, maxlat=maxlat,
                  minlon=minlon, maxlon=maxlon, minmag=minmag, limit=limit, format="xml", **extra)
    mp = _MECH_PARAM.get(service)
    if mp:
        params[mp] = "true"
    if service == "usgs":
        params["format"] = "quakeml"
    params = {k: v for k, v in params.items() if v is not None}
    r = http_get(url, params)
    if not r.content:
        return []
    return parse_quakeml(r.content, catalog=service)
