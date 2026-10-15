"""Unit tests of the SLAM mathematics against analytic results and Cronin's code."""
import numpy as np
import pytest

from geoslam.dem import DEMGrid, GeographicDEM
from geoslam.geodesy import (UTMZone, grid_convergence, grid_convergence_cronin,
                             lonlat_to_utm, scale_factor, utm_to_lonlat, utm_zone)
from geoslam.orientation import (NodalPlane, auxiliary_plane, fold_plane, plane_normal,
                                 plane_normal_cronin)
from geoslam.slam import Hypocenter, seismo_lineament
from geoslam.uncertainty import (HypocenterUncertainty, box_half_width_bayliss,
                                 find_half_width_cronin)

RNG = np.random.default_rng(42)


# ----------------------------------------------------------------- geodesy
def test_utm_round_trip():
    lon = RNG.uniform(-125, -65, 200)
    lat = RNG.uniform(25, 49, 200)
    for lo, la in zip(lon, lat):
        z = utm_zone(lo, la)
        e, n, _ = lonlat_to_utm(lo, la, z)
        lo2, la2 = utm_to_lonlat(e, n, z)
        assert abs(lo2 - lo) < 1e-8 and abs(la2 - la) < 1e-8


def test_utm_central_meridian_and_scale():
    z = UTMZone(10)
    e, n, _ = lonlat_to_utm(-123.0, 0.0, z)
    assert abs(e - 500000.0) < 1e-6 and abs(n) < 1e-6
    assert abs(scale_factor(-123.0, 40.0, z) - 0.9996) < 1e-9
    # northing on the CM equals k0 x meridian arc length (integrated numerically)
    from scipy.integrate import quad
    a, f = 6378137.0, 1 / 298.257223563
    e2 = f * (2 - f)
    arc = quad(lambda p: a * (1 - e2) / (1 - e2 * np.sin(p) ** 2) ** 1.5,
               np.radians(39), np.radians(40))[0]
    _, n1, _ = lonlat_to_utm(-123.0, 39.0, z)
    _, n2, _ = lonlat_to_utm(-123.0, 40.0, z)
    assert abs((n2 - n1) / 0.9996 - arc) < 0.001


def test_zone_and_epsg():
    assert utm_zone(-120.16, 39.438) == UTMZone(10, True)
    assert utm_zone(-97.5, 35.8).epsg == 32614       # Oklahoma (Rasaka: UTM 14N)
    assert utm_zone(-70, -33).epsg == 32719           # southern hemisphere uses 327xx
    assert utm_zone(5.0, 60.0).number == 32           # Norway exception


def test_grid_convergence_matches_cronin_formula():
    # Truckee 1966 epicentre, zone 10 (CM -123)
    g = grid_convergence(-120.16, 39.438, UTMZone(10))
    gc = grid_convergence_cronin(-120.16, 39.438, -123.0)
    assert abs(g - 1.8050) < 2e-3
    assert abs(g - gc) < 1e-3
    # sign: east of CM in the north is positive, west negative
    assert grid_convergence(-124.0, 40.0, UTMZone(10)) < 0


# ------------------------------------------------------------- orientation
def test_normals_match_cronin_construction():
    for _ in range(200):
        d, a = RNG.uniform(1, 89.9), RNG.uniform(0, 359.9)
        assert np.allclose(plane_normal(d, a), plane_normal_cronin(d, a), atol=1e-12)


def test_fold_over_is_same_plane():
    for d, a in [(95, 30), (-10, 200), (100, 359)]:
        d2, a2 = fold_plane(d, a)
        n1, n2 = plane_normal(d, a), plane_normal(d2, a2)
        assert abs(abs(n1 @ n2) - 1) < 1e-12
        assert 0 <= d2 <= 90


@pytest.mark.parametrize("np1,np2", [
    ((345, 85, -170), (254, 80, -5)),      # 2014 South Napa, ComCat nc FPFIT
    ((210, 80, 20), (116, 70, 169)),       # Truckee 1983-07-03
    ((115, 80, 150), (211, 61, 12)),       # Truckee 1992-08-30
    ((63, 44, -90), (243, 46, -90)),       # TexNet texnet2025iqwk
])
def test_auxiliary_plane_matches_catalogs(np1, np2):
    s, d, r = auxiliary_plane(*np1)
    assert abs(((s - np2[0] + 180) % 360) - 180) < 1.5
    assert abs(d - np2[1]) < 1.5
    assert abs(((r - np2[2] + 180) % 360) - 180) < 1.5


# --------------------------------------------------------- uncertainty
def test_half_width_equals_cronin_findHalfWidth():
    for _ in range(300):
        dipaz, dip = RNG.uniform(0, 360), RNG.uniform(1, 89)
        eh1, eh2, ez = RNG.uniform(10, 3000, 3)
        az1 = RNG.uniform(0, 360)
        u = HypocenterUncertainty.from_cronin(eh1, az1, eh2, ez)
        h = u.half_width(plane_normal(dip, dipaz))
        assert abs(h - find_half_width_cronin(dipaz, dip, eh1, az1, eh2, ez)) < 1e-6 * max(1, h)


def test_half_width_is_support_function():
    u = HypocenterUncertainty.from_quakeml_ellipsoid(900, 400, 150, 37, 20, 15)
    # sample the ellipsoid surface densely; max projection onto n is h(n)
    t = RNG.normal(size=(200_000, 3))
    t /= np.linalg.norm(t, axis=1)[:, None]
    pts = (t * u.semi_axes) @ u.axes
    for _ in range(20):
        n = RNG.normal(size=3); n /= np.linalg.norm(n)
        assert abs(np.max(pts @ n) - u.half_width(n)) < 0.01 * u.half_width(n) + 1.0


def test_half_width_limits():
    u = HypocenterUncertainty.isotropic_horizontal(500.0, 2000.0)
    assert abs(u.half_width(plane_normal(90, 10)) - 500.0) < 1e-9  # vertical plane: horizontal error
    assert abs(u.half_width(plane_normal(0, 10)) - 2000.0) < 1e-9  # horizontal plane: vertical error
    # ellipsoid support never exceeds Bayliss' box model
    for d in range(0, 91, 5):
        assert u.half_width(plane_normal(d, 0)) <= box_half_width_bayliss(d, 500, 2000) + 1e-9


def test_rotation_true_to_grid():
    u = HypocenterUncertainty.from_cronin(1000, 50.0, 100, 100)
    g = u.rotated_azimuths(2.0)
    v = g.axes[0]
    assert abs(np.degrees(np.arctan2(v[0], v[1])) - 48.0) < 1e-9


# ---------------------------------------------------------------- SLAM
def flat_dem(half=20_000.0, cell=50.0, elev=0.0, zone=UTMZone(10)):
    n = int(2 * half / cell)
    return DEMGrid(np.full((n, n), elev), 500_000.0 - half, 4_400_000.0 + half, cell, zone)


def hypo_at_center(dem, depth_km, unc=None):
    lon, lat = utm_to_lonlat(500_000.0, 4_400_000.0, dem.zone)   # on the CM: gamma = 0
    return Hypocenter(float(lon), float(lat), depth_km, unc)


def test_mean_trace_position_on_flat_ground():
    dem = flat_dem()
    hyp = hypo_at_center(dem, 8.0)
    for dip, dipaz in [(60, 90), (45, 0), (75, 225)]:
        plane = NodalPlane.from_dip_azimuth(dipaz, dip, 0)
        res = seismo_lineament(dem, hyp, plane)
        assert abs(res.gamma) < 1e-6
        # trace is offset up-dip (towards dipaz + 180) by depth / tan(dip)
        expected = 8000.0 / np.tan(np.radians(dip))
        a = np.radians(dipaz + 180)
        line = max(res.mean_trace(lonlat=False), key=len)
        d = (line[:, 0] - res.focus_xyz[0]) * np.sin(a) + (line[:, 1] - res.focus_xyz[1]) * np.cos(a)
        assert np.allclose(d, expected, atol=1.0)


def test_slab_width_on_flat_ground():
    dem = flat_dem(cell=20.0)
    unc = HypocenterUncertainty.isotropic_horizontal(300.0, 800.0)
    hyp = hypo_at_center(dem, 6.0, unc)
    dip = 70.0
    plane = NodalPlane.from_dip_azimuth(90, dip, 0)   # dips east; strike N-S
    res = seismo_lineament(dem, hyp, plane, method="cronin")
    h = unc.half_width(plane_normal(dip, 90))
    width_cells = res.swath[res.swath.shape[0] // 2].sum()
    assert abs(width_cells * 20.0 - 2 * h / np.sin(np.radians(dip))) <= 2 * 20.0


def test_envelope_contains_cronin_union_and_brute_force():
    dem = flat_dem(half=15_000.0, cell=100.0)
    unc = HypocenterUncertainty.from_cronin(400, 30, 200, 600)
    hyp = hypo_at_center(dem, 7.0, unc)
    plane = NodalPlane.from_dip_azimuth(120, 55, -90, dip_azimuth_uncertainty=15, dip_uncertainty=12)
    cro = seismo_lineament(dem, hyp, plane, method="cronin")
    env = seismo_lineament(dem, hyp, plane, method="envelope")
    assert np.all(env.swath[cro.swath])
    assert env.swath.sum() > cro.swath.sum()
    # brute force: dense lattice over the full rectangle
    P = np.stack(np.meshgrid(*dem.xy()), -1)
    P = np.concatenate([P - env.focus_xyz[:2], (dem.z - env.focus_xyz[2])[..., None]], -1).reshape(-1, 3)
    az0 = env.params["grid_dip_azimuth"]
    ds, az = np.meshgrid(np.linspace(55 - 12, 55 + 12, 97), np.linspace(az0 - 15, az0 + 15, 121))
    nrm = plane_normal(ds.ravel(), az.ravel())
    hw = unc.half_width(nrm)
    brute = np.zeros(len(P), bool)
    for s in range(0, len(P), 20_000):
        brute[s:s + 20_000] = np.any(np.abs(P[s:s + 20_000] @ nrm.T) <= hw, axis=1)
    brute = brute.reshape(dem.shape) | env.slabs["MThin"]
    disagree = (brute ^ env.swath).sum()
    assert disagree <= 0.002 * env.swath.sum()


def test_envelope_steep_fold_over():
    dem = flat_dem(half=10_000.0, cell=100.0)
    unc = HypocenterUncertainty.isotropic_horizontal(200, 300)
    hyp = hypo_at_center(dem, 5.0, unc)
    plane = NodalPlane.from_dip_azimuth(80, 85, 180, 10, 15)   # dip range 70..100 crosses vertical
    res = seismo_lineament(dem, hyp, plane)
    row = res.swath[res.swath.shape[0] // 2]
    xs = dem.xy()[0][row] - res.focus_xyz[0]
    # planes on both sides of vertical -> swath straddles the epicentre symmetrically-ish
    assert xs.min() < -5000 * np.tan(np.radians(10)) + 200
    assert xs.max() > 5000 * np.tan(np.radians(15)) - 200


def test_zero_uncertainty_reduces_to_thin_trace():
    dem = flat_dem(half=8_000.0, cell=25.0)
    hyp = hypo_at_center(dem, 4.0)
    res = seismo_lineament(dem, hyp, NodalPlane.from_dip_azimuth(270, 60, 90))
    assert np.array_equal(res.swath, res.slabs["MThin"])


def test_surface_depth_datum():
    dem = flat_dem(elev=1500.0)
    h1 = hypo_at_center(dem, 5.0)
    h2 = hypo_at_center(dem, 5.0)
    h2.depth_datum = "surface"
    p = NodalPlane.from_dip_azimuth(0, 45, 90)
    r1, r2 = seismo_lineament(dem, h1, p), seismo_lineament(dem, h2, p)
    assert abs(r1.focus_xyz[2] + 5000) < 1e-9 and abs(r2.focus_xyz[2] + 3500) < 1e-6


# ----------------------------------------------------------------- DEM
def test_geographic_to_utm_resampling():
    lon = np.linspace(-120.5, -119.9, 400)
    lat = np.linspace(39.7, 39.2, 300)
    LO, LA = np.meshgrid(lon, lat)
    z = 1000 + 2000 * (LO + 120.5) + 800 * (LA - 39.2)
    dlon, dlat = lon[1] - lon[0], lat[0] - lat[1]
    g = GeographicDEM(z, lon[0] - dlon / 2, lat[0] + dlat / 2, dlon, dlat)
    d = g.to_utm(cellsize=60.0)
    lo2, la2 = d.lonlat_grid()
    expect = 1000 + 2000 * (lo2 + 120.5) + 800 * (la2 - 39.2)
    assert np.nanmax(np.abs(d.z - expect)) < 0.05
    assert d.zone == UTMZone(10)


def test_arc_ascii_round_trip(tmp_path):
    dem = flat_dem(half=1000, cell=100, elev=12.5)
    dem.z[3, 4] = np.nan
    p = tmp_path / "dem.asc"
    dem.to_arc_ascii(p)
    back = DEMGrid.from_arc_ascii(p, dem.zone)
    assert back.shape == dem.shape and np.isnan(back.z[3, 4])
    assert abs(back.x_min - dem.x_min) < 1e-3 and abs(back.y_max - dem.y_max) < 1e-3
