"""One-call SLAM workflow: event -> DEM -> both nodal planes -> maps and GeoJSON.

Mirrors Cronin's lab workflow (SLAMWorkflow.html):

1. get the hypocentre, its uncertainty and the focal mechanism;
2. size the map from depth and the gentlest plausible dip;
3. get a DEM and project it to UTM;
4. run SLAM for *both* nodal planes (the radiation pattern cannot tell the
   fault plane from the auxiliary plane);
5. overlay mapped Quaternary faults and seismicity for interpretation.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .dem import DEMGrid
from .orientation import NodalPlane
from .slam import SLAMResult, seismo_lineament, suggest_half_extent
from .sources.base import Earthquake, FocalMechanism
from .uncertainty import HypocenterUncertainty


@dataclass
class EventAnalysis:
    earthquake: Earthquake
    mechanism: FocalMechanism
    dem: DEMGrid
    results: list[SLAMResult]
    faults: list = field(default_factory=list)

    def plot(self, **kw):
        import matplotlib.pyplot as plt
        from .plotting import plot_slam

        fig, axes = plt.subplots(1, len(self.results), figsize=(8 * len(self.results), 8))
        for ax, r in zip(np.atleast_1d(axes), self.results):
            plot_slam(r, ax=ax, faults=self.faults, **kw)
        fig.tight_layout()
        return fig

    def save(self, outdir):
        out = Path(outdir)
        out.mkdir(parents=True, exist_ok=True)
        feats = []
        for r in self.results:
            feats += r.to_geojson()["features"]
        (out / f"{self.earthquake.event_id}_slam.geojson").write_text(
            json.dumps({"type": "FeatureCollection", "features": feats}))
        (out / f"{self.earthquake.event_id}_summary.txt").write_text(
            "\n\n".join(r.summary() for r in self.results))
        fig = self.plot()
        fig.savefig(out / f"{self.earthquake.event_id}_slam.png", dpi=150)
        return out


def pick_mechanism(eq: Earthquake, prefer=("focal-mechanism", "moment-tensor")) -> FocalMechanism:
    """Prefer first-motion mechanisms (they carry plane uncertainties), then
    the highest-weight product."""
    if not eq.mechanisms:
        raise ValueError(f"{eq.event_id} has no focal mechanism")
    def key(m):
        has_unc = m.np1.strike_uncertainty > 0 or m.np1.dip_uncertainty > 0
        rank = prefer.index(m.product) if m.product in prefer else len(prefer)
        return (not has_unc, rank, -(m.extra.get("weight") or 0))
    return sorted(eq.mechanisms, key=key)[0]


def analyze_event(event: str | Earthquake, dem: DEMGrid | None = None, dem_source="auto",
                  cellsize=None, method="envelope", planes=("NP1", "NP2"),
                  default_plane_uncertainty=None, unknown_location_uncertainty=None,
                  depth_datum="sea_level", half_extent_m=None, faults=True,
                  mechanism: FocalMechanism | None = None, multiplier=1.0) -> EventAnalysis:
    """Run SLAM for an event.

    ``event`` is a ComCat event id (e.g. ``"nc72282711"``) or an
    :class:`Earthquake`.  ``default_plane_uncertainty`` = (strike, dip) degrees
    used when the catalog gives none; ``unknown_location_uncertainty`` = (eh, ez)
    metres used when the hypocentre has none.
    """
    if isinstance(event, str):
        from .sources import comcat
        eq = comcat.get_event(event)
    else:
        eq = event
    fm = mechanism or pick_mechanism(eq)
    chosen: list[NodalPlane] = [p for p in fm.planes if p.label in planes] or fm.planes
    if default_plane_uncertainty is not None:
        for p in chosen:
            if p.strike_uncertainty == 0 and p.dip_uncertainty == 0:
                p.strike_uncertainty, p.dip_uncertainty = default_plane_uncertainty
    unk = None
    if unknown_location_uncertainty is not None:
        unk = HypocenterUncertainty.isotropic_horizontal(*unknown_location_uncertainty,
                                                         source="assumed")
    hyp = eq.hypocenter(depth_datum)
    if dem is None:
        from .sources.dems import get_dem
        half = half_extent_m or max(suggest_half_extent(hyp, p) for p in chosen)
        dem = get_dem(eq.lon, eq.lat, half, source=dem_source, cellsize=cellsize)
    results = [seismo_lineament(dem, hyp, p, method=method, unknown_uncertainty=unk,
                                multiplier=multiplier) for p in chosen]
    flist = []
    if faults:
        from .sources.dems import lonlat_bounds_of_utm_box
        from .sources.faults import faults_in_box
        try:
            flist = faults_in_box(*lonlat_bounds_of_utm_box(dem.zone, (dem.x_min, dem.x_max,
                                                                       dem.y_min, dem.y_max), 0))
        except Exception as exc:  # noqa: BLE001
            eq.notes.append(f"fault download failed: {exc}")
    return EventAnalysis(eq, fm, dem, results, flist)
