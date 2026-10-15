"""WGS84 Universal Transverse Mercator projection and grid convergence.

SLAM is computed in a metric, right-handed Cartesian frame (x = grid east,
y = grid north, z = up).  Cronin's Mathematica code works in UTM metres and
converts azimuths measured from geographic (true) north to grid north with the
grid (meridian) convergence angle::

    gridConvergence = ArcTan[Tan[(midLong - zoneMeridian) Degree] * Sin[midLat Degree]]
    npDipTr         = npDipTrend - gridConvergence

This module implements the projection exactly enough for that purpose with no
GIS dependencies, using the Krüger series (Karney 2011, "Transverse Mercator
with an accuracy of a few nanometers", J. Geodesy 85:475-485; third-order terms,
sub-millimetre inside a UTM zone).  It also returns the *exact* ellipsoidal
convergence angle rather than the spherical approximation, which differ by a
few hundredths of a degree near zone edges.

Cronin's ``convertToUTM`` uses the Clarke 1866 ellipsoid (NAD27).  Modern
catalogs (ComCat, NCEDC, SCEDC, TexNet) and DEMs (SRTM, Copernicus, 3DEP) are
referenced to WGS84 / NAD83, so WGS84 is used here.  Mixing NAD27 epicentres
with a NAD83 DEM produces horizontal errors of tens to ~100 m in CONUS.

Sign convention: ``convergence`` (gamma) is the angle from true north to grid
north, positive clockwise.  It is positive east of the central meridian in the
northern hemisphere, so ``grid_azimuth = true_azimuth - gamma``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# WGS84
A_WGS84 = 6378137.0
F_WGS84 = 1.0 / 298.257223563
K0 = 0.9996
E0 = 500000.0
N0_SOUTH = 10000000.0

_n = F_WGS84 / (2.0 - F_WGS84)
_A = A_WGS84 / (1.0 + _n) * (1.0 + _n**2 / 4.0 + _n**4 / 64.0)
_ALPHA = (
    _n / 2.0 - 2.0 * _n**2 / 3.0 + 5.0 * _n**3 / 16.0,
    13.0 * _n**2 / 48.0 - 3.0 * _n**3 / 5.0,
    61.0 * _n**3 / 240.0,
)
_BETA = (
    _n / 2.0 - 2.0 * _n**2 / 3.0 + 37.0 * _n**3 / 96.0,
    _n**2 / 48.0 + _n**3 / 15.0,
    17.0 * _n**3 / 480.0,
)
_DELTA = (
    2.0 * _n - 2.0 * _n**2 / 3.0 - 2.0 * _n**3,
    7.0 * _n**2 / 3.0 - 8.0 * _n**3 / 5.0,
    56.0 * _n**3 / 15.0,
)
_C = 2.0 * np.sqrt(_n) / (1.0 + _n)


@dataclass(frozen=True)
class UTMZone:
    """A UTM zone: number (1-60) and hemisphere."""

    number: int
    north: bool = True

    @property
    def central_meridian(self) -> float:
        return -183.0 + 6.0 * self.number

    @property
    def epsg(self) -> int:
        """EPSG code of the WGS84 UTM zone (326xx north, 327xx south)."""
        return (32600 if self.north else 32700) + self.number

    @property
    def false_northing(self) -> float:
        return 0.0 if self.north else N0_SOUTH

    def __str__(self) -> str:  # pragma: no cover - cosmetic
        return f"UTM {self.number}{'N' if self.north else 'S'} (EPSG:{self.epsg})"


def utm_zone(lon: float, lat: float) -> UTMZone:
    """UTM zone containing a point, including the Norway/Svalbard exceptions."""
    lon = ((float(lon) + 180.0) % 360.0) - 180.0
    number = int((lon + 180.0) // 6.0) + 1
    if 56.0 <= lat < 64.0 and 3.0 <= lon < 12.0:
        number = 32
    elif 72.0 <= lat < 84.0 and lon >= 0.0:
        if lon < 9.0:
            number = 31
        elif lon < 21.0:
            number = 33
        elif lon < 33.0:
            number = 35
        elif lon < 42.0:
            number = 37
    return UTMZone(min(number, 60), lat >= 0.0)


def _forward(lon, lat, lon0):
    phi = np.radians(np.asarray(lat, dtype=float))
    dlam = np.radians(np.asarray(lon, dtype=float) - lon0)
    sphi = np.sin(phi)
    t = np.sinh(np.arctanh(sphi) - _C * np.arctanh(_C * sphi))
    xi_p = np.arctan2(t, np.cos(dlam))
    eta_p = np.arctanh(np.sin(dlam) / np.sqrt(1.0 + t * t))
    xi, eta = xi_p.copy(), eta_p.copy()
    sigma = np.ones_like(xi_p)
    tau = np.zeros_like(xi_p)
    for j, a in enumerate(_ALPHA, start=1):
        c2, s2 = np.cos(2 * j * xi_p), np.sin(2 * j * xi_p)
        ch2, sh2 = np.cosh(2 * j * eta_p), np.sinh(2 * j * eta_p)
        xi = xi + a * s2 * ch2
        eta = eta + a * c2 * sh2
        sigma = sigma + 2 * j * a * c2 * ch2
        tau = tau + 2 * j * a * s2 * sh2
    gamma = np.arctan2(
        tau * np.sqrt(1 + t * t) + sigma * t * np.tan(dlam),
        sigma * np.sqrt(1 + t * t) - tau * t * np.tan(dlam),
    )
    k = (K0 * _A / A_WGS84) * np.sqrt(
        (1 + ((1 - _n) / (1 + _n) * np.tan(phi)) ** 2)
        * (sigma**2 + tau**2)
        / (t * t + np.cos(dlam) ** 2)
    )
    return K0 * _A * eta, K0 * _A * xi, np.degrees(gamma), k


def lonlat_to_utm(lon, lat, zone: UTMZone | None = None):
    """Project WGS84 lon/lat (degrees) to UTM easting/northing (metres).

    Returns ``(easting, northing, zone)``.  Array inputs are supported; the zone
    is taken from the first point when not given.
    """
    lon_a, lat_a = np.asarray(lon, float), np.asarray(lat, float)
    if zone is None:
        zone = utm_zone(float(lon_a.flat[0]), float(lat_a.flat[0]))
    x, y, _, _ = _forward(lon_a, lat_a, zone.central_meridian)
    return x + E0, y + zone.false_northing, zone


def utm_to_lonlat(easting, northing, zone: UTMZone):
    """Inverse UTM: easting/northing (metres) to WGS84 lon/lat (degrees)."""
    xi = (np.asarray(northing, float) - zone.false_northing) / (K0 * _A)
    eta = (np.asarray(easting, float) - E0) / (K0 * _A)
    xi_p, eta_p = xi.copy(), eta.copy()
    for j, b in enumerate(_BETA, start=1):
        xi_p = xi_p - b * np.sin(2 * j * xi) * np.cosh(2 * j * eta)
        eta_p = eta_p - b * np.cos(2 * j * xi) * np.sinh(2 * j * eta)
    chi = np.arcsin(np.sin(xi_p) / np.cosh(eta_p))
    phi = chi.copy()
    for j, d in enumerate(_DELTA, start=1):
        phi = phi + d * np.sin(2 * j * chi)
    lam = np.arctan2(np.sinh(eta_p), np.cos(xi_p))
    return zone.central_meridian + np.degrees(lam), np.degrees(phi)


def grid_convergence(lon, lat, zone: UTMZone | None = None):
    """Exact ellipsoidal grid convergence gamma (degrees) at a point.

    gamma is positive when grid north lies clockwise (east) of true north.
    Convert an azimuth with ``grid = true - gamma``.
    """
    if zone is None:
        zone = utm_zone(float(np.asarray(lon).flat[0]), float(np.asarray(lat).flat[0]))
    _, _, gamma, _ = _forward(lon, lat, zone.central_meridian)
    return gamma


def grid_convergence_cronin(lon, lat, central_meridian):
    """Cronin's spherical formula, ``atan(tan(lon - lon0) * sin(lat))`` (degrees)."""
    dl = np.radians(np.asarray(lon, float) - central_meridian)
    return np.degrees(np.arctan(np.tan(dl) * np.sin(np.radians(np.asarray(lat, float)))))


def scale_factor(lon, lat, zone: UTMZone | None = None):
    """Point scale factor k of the UTM projection."""
    if zone is None:
        zone = utm_zone(float(np.asarray(lon).flat[0]), float(np.asarray(lat).flat[0]))
    return _forward(lon, lat, zone.central_meridian)[3]


def true_to_grid_azimuth(azimuth_true, gamma):
    """Azimuth relative to grid north, wrapped to [0, 360)."""
    return (np.asarray(azimuth_true, float) - gamma) % 360.0


def grid_to_true_azimuth(azimuth_grid, gamma):
    return (np.asarray(azimuth_grid, float) + gamma) % 360.0
