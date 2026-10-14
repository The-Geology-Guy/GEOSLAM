"""DEM providers, all returning a :class:`~geoslam.dem.DEMGrid` in UTM metres.

=================  ===========================================================
``"3dep"``         USGS 3DEP ImageServer (US only, no key). Requested directly
                   in the UTM zone of the epicentre at the chosen cell size
                   (1/3 arc-second ~10 m data; 1 m lidar where available).
``"opentopography"`` OpenTopography Global DEM API.  **Requires a free API key**
                   (``OPENTOPOGRAPHY_API_KEY`` env var or ``api_key=``).
                   ``demtype`` e.g. COP30 (default), SRTMGL1, NASADEM, AW3D30.
``"opentopography-usgs"`` OpenTopography USGS 3DEP (USGS10m, USGS30m, USGS1m).
``"copernicus"``   Copernicus GLO-30 Cloud-Optimised GeoTIFFs on AWS Open Data
                   (no key, global).  Whole 1x1 degree tiles are cached.
=================  ===========================================================

Vertical datums differ (3DEP NAVD88; Copernicus EGM2008; SRTM EGM96).  The
differences (< ~30 m) shift the hypocentre-to-ground distance by the same
amount; that is negligible next to typical depth uncertainty but is reported in
``DEMGrid.vertical_datum``.

Copernicus and SRTM are *surface* models (canopy, buildings) - for geomorphic
lineament mapping prefer bare-earth lidar or 3DEP.
"""
from __future__ import annotations

import math
import os
import tempfile

import numpy as np

from ..dem import DEMGrid, GeographicDEM, read_geotiff
from ..geodesy import UTMZone, lonlat_to_utm, utm_to_lonlat, utm_zone
from .base import cache_dir, http_get

THREE_DEP = "https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer/exportImage"
OT_GLOBAL = "https://portal.opentopography.org/API/globaldem"
OT_USGS = "https://portal.opentopography.org/API/usgsdem"
COP_URL = "https://copernicus-dem-{res}m.s3.amazonaws.com/{name}/{name}.tif"


def utm_box(lon, lat, half_extent_m, zone: UTMZone | None = None):
    zone = zone or utm_zone(lon, lat)
    e, n, _ = lonlat_to_utm(lon, lat, zone)
    e, n = float(e), float(n)
    return zone, (e - half_extent_m, e + half_extent_m, n - half_extent_m, n + half_extent_m)


def lonlat_bounds_of_utm_box(zone, box, pad_deg=0.01):
    x0, x1, y0, y1 = box
    xs = np.r_[np.linspace(x0, x1, 20), np.full(20, x1), np.linspace(x1, x0, 20), np.full(20, x0)]
    ys = np.r_[np.full(20, y0), np.linspace(y0, y1, 20), np.full(20, y1), np.linspace(y1, y0, 20)]
    lo, la = utm_to_lonlat(xs, ys, zone)
    return lo.min() - pad_deg, la.min() - pad_deg, lo.max() + pad_deg, la.max() + pad_deg


def _in_us(lon, lat):
    return (-125.5 <= lon <= -66 and 24 <= lat <= 50) or (-170 <= lon <= -129 and 51 <= lat <= 72) \
        or (-161 <= lon <= -154 and 18 <= lat <= 23) or (-68 <= lon <= -65 and 17.5 <= lat <= 18.7)


def fetch_3dep(zone, box, cellsize=10.0) -> DEMGrid:
    x0, x1, y0, y1 = box
    w, h = int(round((x1 - x0) / cellsize)), int(round((y1 - y0) / cellsize))
    if w > 8000 or h > 8000:
        raise ValueError("3DEP export limited to 8000x8000 px; increase cellsize or reduce extent")
    meta = http_get(THREE_DEP, dict(bbox=f"{x0},{y0},{x1},{y1}", bboxSR=zone.epsg, imageSR=zone.epsg,
                                    size=f"{w},{h}", format="tiff", pixelType="F32",
                                    noDataInterpretation="esriNoDataMatchAny",
                                    interpolation="RSP_BilinearInterpolation", f="json")).json()
    if "href" not in meta:
        raise RuntimeError(f"3DEP export failed: {meta}")
    tif = http_get(meta["href"]).content
    ext = meta.get("extent", {})
    with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as fh:
        fh.write(tif)
        path = fh.name
    try:
        import tifffile
        z = tifffile.imread(path).astype(float)
    finally:
        os.unlink(path)
    z[(z < -1e4) | (z > 1e5)] = np.nan
    xmin, ymax = ext.get("xmin", x0), ext.get("ymax", y1)
    cs = (ext.get("xmax", x1) - xmin) / z.shape[1]
    return DEMGrid(z, float(xmin), float(ymax), float(cs), zone, "NAVD88", "USGS 3DEP")


def _download_geotiff(url, params=None, source=""):
    data = http_get(url, params, cache=True, timeout=300)
    with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as fh:
        fh.write(data)
        path = fh.name
    try:
        return read_geotiff(path, source=source)
    finally:
        os.unlink(path)


def fetch_opentopography(bounds, demtype="COP30", api_key=None, usgs=False):
    api_key = api_key or os.environ.get("OPENTOPOGRAPHY_API_KEY")
    if not api_key:
        raise RuntimeError("OpenTopography requires an API key: set OPENTOPOGRAPHY_API_KEY "
                           "(free at https://portal.opentopography.org/requestService)")
    w, s, e, n = bounds
    params = dict(south=s, north=n, west=w, east=e, outputFormat="GTiff", API_Key=api_key)
    if usgs:
        params["datasetName"] = demtype
        return _download_geotiff(OT_USGS, params, f"OpenTopography {demtype}")
    params["demtype"] = demtype
    return _download_geotiff(OT_GLOBAL, params, f"OpenTopography {demtype}")


def _cop_tile_name(lat0, lon0, res=30):
    ns = f"{'N' if lat0 >= 0 else 'S'}{abs(lat0):02d}_00"
    ew = f"{'E' if lon0 >= 0 else 'W'}{abs(lon0):03d}_00"
    return f"Copernicus_DSM_COG_{10 if res == 30 else 30}_{ns}_{ew}_DEM"


def fetch_copernicus(bounds, res=30) -> GeographicDEM:
    """Mosaic Copernicus GLO-30/90 tiles covering ``bounds`` (w, s, e, n)."""
    w, s, e, n = bounds
    tiles = []
    for la in range(math.floor(s), math.floor(n) + 1):
        for lo in range(math.floor(w), math.floor(e) + 1):
            name = _cop_tile_name(la, lo, res)
            tiles.append(_download_geotiff(COP_URL.format(res=res, name=name), source="Copernicus GLO-30"))
    if len(tiles) == 1:
        g = tiles[0]
    else:
        g = _mosaic(tiles)
    g.vertical_datum = "EGM2008"
    return g


def _mosaic(tiles):
    dlat = min(t.dlat for t in tiles)
    dlon = min(t.dlon for t in tiles)
    lo0 = min(t.lon_min for t in tiles)
    la1 = max(t.lat_max for t in tiles)
    lo1 = max(t.bounds[2] for t in tiles)
    la0 = min(t.bounds[1] for t in tiles)
    nc, nr = int(round((lo1 - lo0) / dlon)), int(round((la1 - la0) / dlat))
    z = np.full((nr, nc), np.nan)
    for t in tiles:
        if abs(t.dlon - dlon) > 1e-12 or abs(t.dlat - dlat) > 1e-12:
            from scipy.ndimage import zoom
            t = GeographicDEM(zoom(t.z, (t.dlat / dlat, t.dlon / dlon), order=1), t.lon_min,
                              t.lat_max, dlon, dlat)
        r0 = int(round((la1 - t.lat_max) / dlat))
        c0 = int(round((t.lon_min - lo0) / dlon))
        z[r0:r0 + t.z.shape[0], c0:c0 + t.z.shape[1]] = t.z[:nr - r0, :nc - c0]
    return GeographicDEM(z, lo0, la1, dlon, dlat, source=tiles[0].source)


def get_dem(lon, lat, half_extent_m, source="auto", cellsize=None, zone=None, **kw) -> DEMGrid:
    """DEM in UTM metres centred on (lon, lat), ``2*half_extent_m`` on a side."""
    zone, box = utm_box(lon, lat, half_extent_m, zone)
    if source == "auto":
        if _in_us(lon, lat):
            source = "3dep"
        elif os.environ.get("OPENTOPOGRAPHY_API_KEY") or kw.get("api_key"):
            source = "opentopography"
        else:
            source = "copernicus"
    if source == "3dep":
        return fetch_3dep(zone, box, cellsize or 10.0)
    bounds = lonlat_bounds_of_utm_box(zone, box)
    if source == "opentopography":
        g = fetch_opentopography(bounds, kw.get("demtype", "COP30"), kw.get("api_key"))
    elif source == "opentopography-usgs":
        g = fetch_opentopography(bounds, kw.get("demtype", "USGS10m"), kw.get("api_key"), usgs=True)
        if isinstance(g, DEMGrid):
            return g.crop(*box)
    elif source == "copernicus":
        g = fetch_copernicus(bounds, kw.get("res", 30))
    else:
        raise ValueError(f"unknown DEM source {source!r}")
    if isinstance(g, DEMGrid):
        return g.crop(*box)
    return g.to_utm(cellsize=cellsize, zone=zone).crop(*box)
