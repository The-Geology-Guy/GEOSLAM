"""Core of the Seismo-Lineament Analysis Method (SLAM).

Method (Cronin et al., 2008, Environ. Eng. Geosci. 14(3):199-219; Cronin,
SLAMcode20210103.nb):

1.  Put the hypocentre and a DEM in one metric frame: UTM metres, z up,
    hypocentre depth relative to sea level (``z_h = -depth``).
2.  Rotate the nodal plane's dip azimuth from true to grid north
    (``grid = true - gamma``).
3.  A candidate fault plane passes through the hypocentre with upward unit
    normal ``n``.  A ground point ``P`` (relative to the hypocentre) lies on the
    plane's *uncertainty slab* when ``|n . P| <= h(n)``, where ``h`` is the
    distance from the hypocentre to the planes parallel to the fault plane and
    tangent to the hypocentre uncertainty ellipsoid.
4.  The seismo-lineament is the set of ground points that lie on the slab of
    at least one plausible fault plane, i.e. any plane whose dip and dip
    azimuth lie within their published uncertainties.

Two implementations of step 4 are provided:

``method="cronin"``
    Exactly Cronin's notebook: the union of seven slabs - the mean plane with
    a thin width ``widthFactor * cellsize * sin(dip)`` (his "MThin" trace),
    the four corner planes (dip -/+ sigma_dip, dip azimuth -/+ sigma_az: GS,
    GL, SS, SL), the steep/mean-azimuth plane (SM) and the mean plane with
    the ellipsoid width.  Use this to reproduce published SLAM maps.

``method="envelope"`` (default)
    The full envelope over the continuous rectangle of plausible planes.  A
    point is inside if some plane in the rectangle passes through it exactly
    (the signed distance changes sign across the rectangle) or if it lies on
    the slab of a plane on the rectangle's boundary.  Because ``n . P`` has no
    interior minimum of ``|n . P|`` other than its zeros, this is the true
    envelope up to the (slow) variation of ``h``.  This replaces both the gaps
    between Cronin's seven discrete slabs far from the epicentre and the
    ad-hoc ``fill_swaths`` / ``fill_corners`` heuristics in GEOSLAM 0.0.x.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from .dem import DEMGrid
from .geodesy import grid_convergence, lonlat_to_utm, utm_to_lonlat
from .orientation import NodalPlane, plane_normal
from .uncertainty import HypocenterUncertainty


@dataclass
class Hypocenter:
    """Earthquake focus.  ``depth_km`` is positive down.

    ``uncertainty`` is the location confidence ellipsoid with azimuths relative
    to *true* north (as catalogs report them); it is rotated to grid north
    internally.  ``depth_datum`` is ``"sea_level"`` (Cronin's convention) or
    ``"surface"`` (depth below the ground at the epicentre).
    """

    lon: float
    lat: float
    depth_km: float
    uncertainty: HypocenterUncertainty | None = None
    depth_datum: Literal["sea_level", "surface"] = "sea_level"
    event_id: str = ""
    time: str = ""
    magnitude: float | None = None


@dataclass
class PlaneSpec:
    name: str
    dip: float           # may be <0 or >90 (equivalent folded plane)
    dip_azimuth: float   # grid azimuth
    half_width: float    # metres


@dataclass
class SLAMResult:
    dem: DEMGrid
    hypocenter: Hypocenter
    plane: NodalPlane
    method: str
    gamma: float                      # grid convergence at the epicentre (deg)
    focus_xyz: np.ndarray             # hypocentre in UTM metres (z relative to sea level)
    planes: list[PlaneSpec]
    swath: np.ndarray                 # bool grid: the seismo-lineament
    slabs: dict[str, np.ndarray]      # bool grid per discrete plane
    mean_signed_distance: np.ndarray  # n0 . P for the mean plane (m); its 0-contour is the trace
    params: dict = field(default_factory=dict)
    signed_distances: dict[str, np.ndarray] = field(default_factory=dict, repr=False)
    margin: np.ndarray | None = field(default=None, repr=False)
    """Continuous field (m): negative inside the swath, positive outside.  Its 0-contour
    is a sub-cell, smooth version of the swath boundary (used for plotting)."""

    # ------------------------------------------------------------ geometry
    @property
    def area_km2(self) -> float:
        return float(self.swath.sum() * self.dem.cellsize**2 / 1e6)

    def _contours(self, field_, level):
        from contourpy import contour_generator

        x, y = self.dem.xy()
        gen = contour_generator(x, y, np.nan_to_num(field_, nan=-1e9))
        return gen.lines(level)

    def boundary_lines(self, which="swath", lonlat=True):
        """Boundary polylines of the swath (or a named slab) as (N,2) arrays."""
        arr = self.swath if which == "swath" else self.slabs[which]
        lines = self._contours(arr.astype(float), 0.5)
        if lonlat:
            lines = [np.column_stack(utm_to_lonlat(l[:, 0], l[:, 1], self.dem.zone)) for l in lines]
        return lines

    def mean_trace(self, lonlat=True):
        """Ground trace of the mean nodal plane (zero contour of n0 . P)."""
        lines = self._contours(self.mean_signed_distance, 0.0)
        if lonlat:
            lines = [np.column_stack(utm_to_lonlat(l[:, 0], l[:, 1], self.dem.zone)) for l in lines]
        return lines

    def to_geojson(self) -> dict:
        """Swath boundary and mean trace as a GeoJSON FeatureCollection (WGS84)."""
        feats = []
        props = {"event_id": self.hypocenter.event_id, "method": self.method,
                 "strike": self.plane.strike, "dip": self.plane.dip, "rake": self.plane.rake,
                 "label": self.plane.label}
        for line in self.boundary_lines():
            feats.append({"type": "Feature", "properties": dict(props, kind="swath_boundary"),
                          "geometry": {"type": "LineString", "coordinates": line.round(7).tolist()}})
        for line in self.mean_trace():
            feats.append({"type": "Feature", "properties": dict(props, kind="mean_trace"),
                          "geometry": {"type": "LineString", "coordinates": line.round(7).tolist()}})
        feats.append({"type": "Feature", "properties": dict(props, kind="epicenter"),
                      "geometry": {"type": "Point",
                                   "coordinates": [self.hypocenter.lon, self.hypocenter.lat]}})
        return {"type": "FeatureCollection", "features": feats}

    def summary(self) -> str:
        p = self.plane
        rows = [f"SLAM ({self.method}) event {self.hypocenter.event_id or '?'}  "
                f"plane {p.label or ''} strike/dip/rake {p.strike:.1f}/{p.dip:.1f}/{p.rake}",
                f"  zone {self.dem.zone}, grid convergence {self.gamma:+.3f} deg, "
                f"DEM {self.dem.shape[1]}x{self.dem.shape[0]} @ {self.dem.cellsize:g} m",
                f"  focus (E, N, z) = ({self.focus_xyz[0]:.1f}, {self.focus_xyz[1]:.1f}, {self.focus_xyz[2]:.1f}) m",
                f"  swath area {self.area_km2:.2f} km^2"]
        for s in self.planes:
            rows.append(f"  {s.name:>5}: dip {s.dip:6.2f}  dip-az(grid) {s.dip_azimuth % 360:7.2f}  half-width {s.half_width:8.1f} m")
        return "\n".join(rows)


# ----------------------------------------------------------------------------
def _ground_points(dem: DEMGrid, focus):
    x, y = dem.xy()
    xx, yy = np.meshgrid(x - focus[0], y - focus[1])
    return np.stack([xx, yy, dem.z - focus[2]], axis=-1)


def _focus(dem: DEMGrid, hypo: Hypocenter):
    e, n, _ = lonlat_to_utm(hypo.lon, hypo.lat, dem.zone)
    depth = hypo.depth_km * 1000.0
    if hypo.depth_datum == "surface":
        z = float(dem.elevation_at(e, n)[0]) - depth
    else:
        z = -depth
    return np.array([float(e), float(n), z])


def cronin_planes(dip, dip_az, sd, sa):
    """Cronin's six orientation variants (MThin and Mean share the mean orientation).

    Returns (name, dip, dip_azimuth) using continuous angles; Cronin's fold-over
    rule (|dip - sd| with azimuth +180, or 180 - (dip + sd) with azimuth +180)
    describes the same planes.
    """
    return [
        ("GS", dip - sd, dip_az - sa),
        ("GL", dip - sd, dip_az + sa),
        ("SS", dip + sd, dip_az - sa),
        ("SL", dip + sd, dip_az + sa),
        ("SM", dip + sd, dip_az),
        ("Mean", dip, dip_az),
    ]


def seismo_lineament(dem: DEMGrid, hypocenter: Hypocenter, plane: NodalPlane,
                     method: Literal["envelope", "cronin"] = "envelope",
                     width_factor: float = 1.5, multiplier: float = 1.0,
                     angular_step: float = 0.5, chunk: int = 250_000,
                     unknown_uncertainty: HypocenterUncertainty | None = None) -> SLAMResult:
    """Compute the seismo-lineament of one nodal plane on a UTM DEM.

    Parameters
    ----------
    dem : DEMGrid in UTM metres; its zone defines grid north.
    hypocenter : Hypocenter (lon/lat WGS84, depth km).
    plane : NodalPlane with strike/dip (true north) and their uncertainties.
    method : ``"envelope"`` (continuous) or ``"cronin"`` (7 discrete slabs).
    width_factor : Cronin's ``widthFactor`` for the thin mean-plane trace.
    multiplier : scale applied to every ellipsoid half-width (1 = as published).
    angular_step : sampling step (deg) along the uncertainty rectangle edges.
    unknown_uncertainty : used when the hypocentre has none.  If both are
        missing the swath degenerates to the thin mean-plane trace.
    """
    gamma = float(grid_convergence(hypocenter.lon, hypocenter.lat, dem.zone))
    focus = _focus(dem, hypocenter)
    unc = hypocenter.uncertainty or unknown_uncertainty
    unc_grid = unc.rotated_azimuths(gamma).scaled(multiplier) if unc is not None else None

    dip0 = float(plane.dip)
    az0 = (plane.dip_azimuth - gamma) % 360.0
    sd, sa = float(plane.dip_uncertainty), float(plane.strike_uncertainty)

    P = _ground_points(dem, focus)
    shape = P.shape[:2]
    Pf = P.reshape(-1, 3)
    valid = ~np.isnan(Pf[:, 2])

    def half_width(normals):
        if unc_grid is None:
            return np.zeros(normals.shape[:-1])
        return unc_grid.half_width(normals)

    n0 = plane_normal(dip0, az0)
    mean_sd = (Pf @ n0).reshape(shape)
    thin_w = width_factor * dem.cellsize * abs(np.sin(np.radians(dip0)))

    planes: list[PlaneSpec] = [PlaneSpec("MThin", dip0, az0, thin_w)]
    slabs: dict[str, np.ndarray] = {"MThin": np.abs(mean_sd) <= thin_w}
    signed: dict[str, np.ndarray] = {}
    for name, d, a in cronin_planes(dip0, az0, sd, sa):
        n = plane_normal(d, a)
        h = float(half_width(n))
        planes.append(PlaneSpec(name, d, a, h))
        sdist = (Pf @ n).reshape(shape)
        signed[name] = sdist
        slabs[name] = np.abs(sdist) <= h

    if method == "cronin":
        swath = np.zeros(shape, bool)
        for m in slabs.values():
            swath |= m
        margin = np.abs(mean_sd) - thin_w
        for ps in planes[1:]:
            margin = np.minimum(margin, np.abs(signed[ps.name]) - ps.half_width)
    elif method == "envelope":
        env, margin = _envelope(Pf, dip0, az0, sd, sa, half_width, angular_step, chunk)
        swath = env.reshape(shape) | slabs["MThin"]
        margin = np.minimum(margin.reshape(shape), np.abs(mean_sd) - thin_w)
    else:
        raise ValueError("method must be 'envelope' or 'cronin'")

    swath &= valid.reshape(shape)
    return SLAMResult(dem, hypocenter, plane, method, gamma, focus, planes, swath, slabs,
                      mean_sd, dict(width_factor=width_factor, multiplier=multiplier,
                                    angular_step=angular_step, grid_dip_azimuth=az0), signed,
                      np.where(valid.reshape(shape), margin, np.nan))


def _rect_boundary(d0, a0, sd, sa, step):
    """(dip, az) samples along the boundary of the uncertainty rectangle, plus a
    coarse interior lattice for robustness."""
    if sd == 0 and sa == 0:
        return np.array([[d0, a0]])
    nd = max(int(np.ceil(2 * sd / step)), 1)
    na = max(int(np.ceil(2 * sa / step)), 1)
    ds = np.linspace(d0 - sd, d0 + sd, nd + 1)
    az = np.linspace(a0 - sa, a0 + sa, na + 1)
    edges = [np.column_stack([ds, np.full_like(ds, a0 - sa)]),
             np.column_stack([ds, np.full_like(ds, a0 + sa)]),
             np.column_stack([np.full_like(az, d0 - sd), az]),
             np.column_stack([np.full_like(az, d0 + sd), az])]
    gi, ga = np.meshgrid(np.linspace(d0 - sd, d0 + sd, 5), np.linspace(a0 - sa, a0 + sa, 5))
    edges.append(np.column_stack([gi.ravel(), ga.ravel()]))
    return np.unique(np.vstack(edges).round(9), axis=0)


def _in_rect(d, a, d0, a0, sd, sa):
    """Is (d, a) - or an equivalent parametrisation - inside the rectangle?"""
    out = np.zeros(np.shape(d), bool)
    # n(d, a) == n(-d, a + 180): the only other parametrisation of one normal
    for dd, aa in ((d, a), (-d, a + 180.0)):
        da = ((aa - a0 + 180.0) % 360.0) - 180.0
        out |= (np.abs(dd - d0) <= sd + 1e-9) & (np.abs(da) <= sa + 1e-9)
    return out


def _envelope(Pf, d0, a0, sd, sa, half_width, step, chunk):
    samples = _rect_boundary(d0, a0, sd, sa, step)
    normals = plane_normal(samples[:, 0], samples[:, 1])   # (K, 3)
    widths = np.atleast_1d(half_width(normals))            # (K,)
    # float32 halves memory traffic; positions are relative to the focus (< ~1e5 m),
    # so single precision still resolves distances to ~1 cm.
    normals_t = normals.T.astype(np.float32)
    widths32 = widths.astype(np.float32)
    out = np.zeros(len(Pf), bool)
    margin = np.zeros(len(Pf))
    hmin = float(widths.min())
    d_lo, d_hi = d0 - sd, d0 + sd
    for s in range(0, len(Pf), chunk):
        p = Pf[s:s + chunk]
        f = p.astype(np.float32) @ normals_t                 # (m, K)
        slab_margin = np.min(np.abs(f) - widths32, axis=1).astype(float)
        inside = slab_margin <= 0.0
        fmin, fmax = f.min(axis=1), f.max(axis=1)
        # interior critical points of f: n parallel to +/-P (f = +/-|P|)
        r = np.linalg.norm(p, axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            dip_p = np.degrees(np.arccos(np.clip(p[:, 2] / r, -1, 1)))
        az_p = np.degrees(np.arctan2(p[:, 0], p[:, 1]))
        plus = _in_rect(dip_p, az_p, d0, a0, sd, sa)           # n = +P/|P|  -> f = +|P|
        minus = _in_rect(dip_p - 180.0, az_p, d0, a0, sd, sa)  # n = -P/|P|  -> f = -|P|
        fmax = np.where(plus, r, fmax)
        fmin = np.where(minus, -r, fmin)
        # saddle on a horizontal plane (dip 0 or 180 inside the dip range): f = +/-Pz
        if d_lo <= 0.0 <= d_hi:
            fmin, fmax = np.minimum(fmin, p[:, 2]), np.maximum(fmax, p[:, 2])
        if d_lo <= 180.0 <= d_hi or d_lo <= -180.0 <= d_hi:
            fmin, fmax = np.minimum(fmin, -p[:, 2]), np.maximum(fmax, -p[:, 2])
        crosses = (fmin <= 0.0) & (fmax >= 0.0)
        out[s:s + chunk] = inside | crosses
        # smooth margin: outside the crossing region it is the distance beyond the
        # nearest bounding slab; inside, depth into the crossing region (offset by -h)
        margin[s:s + chunk] = np.where(crosses, -np.minimum(fmax, -fmin) - hmin, slab_margin)
    return out, margin


def suggest_half_extent(hypocenter: Hypocenter, plane: NodalPlane, mean_elevation_m=0.0,
                        min_m=5_000.0, max_m=60_000.0, margin=1.5) -> float:
    """Half-width (m) of a map that contains the seismo-lineament.

    Cronin's workflow: d = (mean epicentral elevation + |depth|) / tan(gentlest
    plausible dip), times a safety margin.
    """
    gentle = max(abs(plane.dip) - abs(plane.dip_uncertainty), 5.0)
    d = (mean_elevation_m + hypocenter.depth_km * 1000.0) / np.tan(np.radians(gentle))
    if hypocenter.uncertainty is not None:
        d += float(hypocenter.uncertainty.semi_axes.max())
    return float(np.clip(margin * d, min_m, max_m))
