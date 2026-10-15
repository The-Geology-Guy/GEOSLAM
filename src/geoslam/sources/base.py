"""Common data model for earthquakes from any catalog, plus small HTTP helpers."""
from __future__ import annotations

import hashlib
import os
import time
from dataclasses import dataclass, field
from pathlib import Path

from ..orientation import NodalPlane
from ..slam import Hypocenter
from ..uncertainty import HypocenterUncertainty

USER_AGENT = "geoslam/1.0 (+https://github.com/The-Geology-Guy/GEOSLAM)"


def cache_dir() -> Path:
    p = Path(os.environ.get("GEOSLAM_CACHE", Path.home() / ".cache" / "geoslam"))
    p.mkdir(parents=True, exist_ok=True)
    return p


def http_get(url, params=None, timeout=90, retries=3, stream=False, cache=False):
    """GET with a descriptive User-Agent, retries and optional on-disk caching.

    Returns ``requests.Response`` (or bytes when ``cache=True``).
    """
    import requests

    key = None
    if cache:
        key = hashlib.sha1((url + repr(sorted((params or {}).items()))).encode()).hexdigest()
        f = cache_dir() / "http" / key
        if f.exists():
            return f.read_bytes()
    last = None
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=timeout, stream=stream,
                             headers={"User-Agent": USER_AGENT})
            if r.status_code == 204:
                r._content = b""
                return r.content if cache else r
            if r.status_code in (429, 500, 502, 503, 504):
                raise requests.HTTPError(f"{r.status_code} from {url}", response=r)
            r.raise_for_status()
            if cache:
                f = cache_dir() / "http" / key
                f.parent.mkdir(exist_ok=True)
                f.write_bytes(r.content)
                return r.content
            return r
        except Exception as exc:  # noqa: BLE001 - network errors of all kinds
            last = exc
            if getattr(getattr(exc, "response", None), "status_code", 0) in (400, 401, 403, 404):
                raise
            time.sleep(2 ** attempt)
    raise last


@dataclass
class FocalMechanism:
    """A double-couple solution: two nodal planes (true-north azimuths)."""

    np1: NodalPlane
    np2: NodalPlane
    source: str = ""          # e.g. "nc", "us", "ci", "tx"
    product: str = ""         # "focal-mechanism" | "moment-tensor" | "fdsn"
    method: str = ""          # FPFIT, HASH, Mww, TMTS ...
    quality: str = ""
    preferred_plane: int | None = None
    quakeml_url: str = ""
    extra: dict = field(default_factory=dict)

    @property
    def planes(self):
        return [self.np1, self.np2]


@dataclass
class Earthquake:
    event_id: str
    time: str
    lon: float
    lat: float
    depth_km: float
    magnitude: float | None = None
    magnitude_type: str = ""
    uncertainty: HypocenterUncertainty | None = None
    mechanisms: list[FocalMechanism] = field(default_factory=list)
    catalog: str = ""
    url: str = ""
    notes: list[str] = field(default_factory=list)

    def hypocenter(self, depth_datum="sea_level", uncertainty=None) -> Hypocenter:
        return Hypocenter(self.lon, self.lat, self.depth_km, uncertainty or self.uncertainty,
                          depth_datum, self.event_id, self.time, self.magnitude)

    def as_row(self) -> dict:
        """Flat record in the spirit of Cronin's input spreadsheet."""
        row = dict(event_id=self.event_id, time=self.time, catalog=self.catalog,
                   lat=self.lat, lon=self.lon, depth_km=self.depth_km,
                   mag=self.magnitude, mag_type=self.magnitude_type)
        u = self.uncertainty
        if u is not None:
            row.update(err_semi_axes_m=tuple(round(float(a), 1) for a in u.semi_axes),
                       err_source=u.source)
        for k, m in enumerate(self.mechanisms):
            for j, p in enumerate(m.planes, start=1):
                row[f"fm{k}_np{j}"] = (p.strike, p.dip, p.rake)
                row[f"fm{k}_np{j}_unc"] = (p.strike_uncertainty, p.dip_uncertainty, p.rake_uncertainty)
            row[f"fm{k}_src"] = f"{m.source}:{m.product}:{m.method}"
        return row
