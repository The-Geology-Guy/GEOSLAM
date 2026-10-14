"""Maps of seismo-lineaments (matplotlib only; no basemap tiles required).

The base map is a hillshade of the same UTM DEM used for the computation, lit
as Cronin does: from the direction opposite the nodal plane's dip (so scarps
facing the light along the fault are emphasised) at a low sun angle (15 deg,
``Pi/12`` in his notebook).  The 2020 port passed a math-convention angle to
matplotlib's compass-convention ``azdeg``; here the compass azimuth is used.
"""
from __future__ import annotations

import numpy as np

from .geodesy import lonlat_to_utm, utm_to_lonlat
from .orientation import plane_normal, slip_vector

SLAB_STYLE = {"GS": ":", "GL": ":", "SS": "--", "SL": "--", "SM": "-.", "Mean": "-", "MThin": "-"}
"""Optional dash patterns per plane, used only with ``trace_style="dashed"``."""


def light_azimuth(result) -> float:
    """Cronin: light from opposite the (grid) dip direction."""
    return (result.params["grid_dip_azimuth"] + 180.0) % 360.0


def hillshade(dem, azimuth=315.0, altitude=15.0, vert_exag=1.0):
    from matplotlib.colors import LightSource

    z = dem.z.copy()
    nan = np.isnan(z)
    if nan.any():
        z[nan] = np.nanmean(z)
    hs = LightSource(azdeg=azimuth, altdeg=altitude).hillshade(z, vert_exag=vert_exag,
                                                               dx=dem.cellsize, dy=dem.cellsize)
    hs[nan] = np.nan
    return hs


def beachball_image(strike, dip, rake, n=201):
    """Lower-hemisphere equal-area beachball as an RGBA image (compressional = dark)."""
    u = np.linspace(-1, 1, n)
    X, Y = np.meshgrid(u, -u)
    R = np.hypot(X, Y)
    inside = R <= 1.0
    # equal-area (Schmidt) inverse: r = sqrt(2) sin(i/2), i = angle from vertical
    inc = 2 * np.arcsin(np.clip(R / np.sqrt(2), 0, 1))
    az = np.arctan2(X, Y)
    ray = np.stack([np.sin(inc) * np.sin(az), np.sin(inc) * np.cos(az), -np.cos(inc)], -1)
    nrm = plane_normal(dip, (strike + 90) % 360)
    slp = slip_vector(strike, dip, rake)
    pol = (ray @ nrm) * (ray @ slp)
    img = np.ones((n, n, 4))
    img[..., :3] = np.where((pol > 0)[..., None], 0.1, 1.0)
    img[..., 3] = inside.astype(float)
    edge = np.abs(R - 1.0) < 1.5 / n
    img[edge] = (0, 0, 0, 1)
    return img


def plot_slam(result, ax=None, faults=None, show_slabs=True, show_hillshade=True,
              aftershocks=None, beachball=True, title=None, km=True, fault_color="#d62728",
              swath_color="#1f77b4", light_altitude=15.0, trace_style="solid"):
    """Map of one SLAM result.  ``faults`` is the output of
    :func:`geoslam.sources.faults.faults_in_box`; ``aftershocks`` an (N, 2+)
    array of lon, lat[, magnitude].

    All lines are drawn solid by default.  ``trace_style="dashed"`` restores
    per-plane dash patterns (dotted = gentle, dashed = steep, dash-dot = SM).
    """
    import matplotlib.pyplot as plt

    dem = result.dem
    if ax is None:
        _, ax = plt.subplots(figsize=(9, 9))
    s = 1e-3 if km else 1.0
    ext = [v * s for v in dem.extent]
    if show_hillshade:
        ax.imshow(hillshade(dem, light_azimuth(result), light_altitude), cmap="gray",
                  extent=ext, origin="upper", interpolation="bilinear")
    x, y = dem.xy()
    X, Y = np.meshgrid(x * s, y * s)
    widths = {p_.name: p_.half_width for p_ in result.planes}
    if result.method == "cronin" and result.signed_distances:
        # Each slab is the band -h <= n.P <= h of a smooth field, so contouring the field
        # itself gives continuous, sub-cell lines even for slabs thinner than one DEM cell.
        bands = [("MThin", result.mean_signed_distance)] + list(result.signed_distances.items())
        for name, sd in bands:
            h = widths[name]
            sdf = np.nan_to_num(sd)
            ax.contourf(X, Y, sdf, levels=[-h, h], colors=[swath_color], alpha=0.25)
            ax.contour(X, Y, sdf, levels=[-h, h], colors=[swath_color], linewidths=1.6,
                       linestyles="solid")
    elif result.margin is not None:
        m = np.nan_to_num(result.margin, nan=1e9)
        ax.contourf(X, Y, m, levels=[np.nanmin(m) - 1.0, 0.0], colors=[swath_color], alpha=0.25)
        ax.contour(X, Y, m, levels=[0.0], colors=[swath_color], linewidths=2.0, linestyles="solid")
    else:
        ax.contourf(X, Y, result.swath.astype(float), levels=[0.5, 1.5], colors=[swath_color], alpha=0.25)
        ax.contour(X, Y, result.swath.astype(float), levels=[0.5], colors=[swath_color], linewidths=2.0,
                   linestyles="solid")
    if show_slabs:
        # ground traces of Cronin's bounding planes (gentle/steep x smaller/larger azimuth)
        for name, sd in result.signed_distances.items():
            if name == "Mean":
                continue
            ax.contour(X, Y, np.nan_to_num(sd), levels=[0.0], colors=["k"],
                       linewidths=0.7, alpha=0.55,
                       linestyles=SLAB_STYLE.get(name, "-") if trace_style == "dashed" else "solid")
    ax.contour(X, Y, np.nan_to_num(result.mean_signed_distance), levels=[0.0], colors=["k"],
               linewidths=1.8, linestyles="solid")
    if faults:
        for f in faults:
            for line in f["lines"]:
                fx, fy, _ = lonlat_to_utm(line[:, 0], line[:, 1], dem.zone)
                ax.plot(fx * s, fy * s, color=fault_color, lw=1.4, alpha=0.9)
    if aftershocks is not None and len(aftershocks):
        a = np.asarray(aftershocks, float)
        axx, ayy, _ = lonlat_to_utm(a[:, 0], a[:, 1], dem.zone)
        size = 4 * np.clip(a[:, 2], 0.5, 8) ** 2 if a.shape[1] > 2 else 12
        ax.scatter(axx * s, ayy * s, s=size, facecolors="none", edgecolors="#ff7f0e", lw=0.8)
    fx, fy = result.focus_xyz[:2] * s
    ax.plot(fx, fy, marker="*", ms=16, mfc="yellow", mec="k", zorder=6)
    if beachball and result.plane.rake is not None:
        p = result.plane
        img = beachball_image(p.strike, p.dip, p.rake)
        w = 0.07 * (ext[1] - ext[0])
        ax.imshow(img, extent=[fx + 0.4 * w, fx + 1.4 * w, fy + 0.4 * w, fy + 1.4 * w], zorder=7)
    ax.set_xlim(ext[0], ext[1])
    ax.set_ylim(ext[2], ext[3])
    ax.set_aspect("equal")
    unit = "km" if km else "m"
    ax.set_xlabel(f"{dem.zone} easting ({unit})")
    ax.set_ylabel(f"northing ({unit})")
    _graticule(ax, dem, s)
    p = result.plane
    ax.set_title(title or (f"{result.hypocenter.event_id}  {p.label}  "
                           f"{p.strike:.0f}/{p.dip:.0f}/{p.rake if p.rake is not None else '?'}"
                           f"  (±{p.strike_uncertainty:.0f}/±{p.dip_uncertainty:.0f})  [{result.method}]"),
                 fontsize=10)
    return ax


def _graticule(ax, dem, s, n=4):
    lo, la = dem.lonlat_grid()
    lons = np.linspace(np.nanmin(lo), np.nanmax(lo), n + 2)[1:-1]
    lats = np.linspace(np.nanmin(la), np.nanmax(la), n + 2)[1:-1]
    lat_line = np.linspace(np.nanmin(la), np.nanmax(la), 50)
    lon_line = np.linspace(np.nanmin(lo), np.nanmax(lo), 50)
    for L in lons:
        x, y, _ = lonlat_to_utm(np.full(50, L), lat_line, dem.zone)
        ax.plot(x * s, y * s, color="w", lw=0.4, alpha=0.5)
        ax.annotate(f"{L:.2f}°", (x[0] * s, dem.y_min * s), color="0.3", fontsize=7,
                    ha="center", va="bottom")
    for L in lats:
        x, y, _ = lonlat_to_utm(lon_line, np.full(50, L), dem.zone)
        ax.plot(x * s, y * s, color="w", lw=0.4, alpha=0.5)
        ax.annotate(f"{L:.2f}°", (dem.x_min * s, y[0] * s), color="0.3", fontsize=7,
                    ha="left", va="center")
