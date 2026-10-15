"""Digital elevation model handling: a UTM-metre grid for SLAM.

Cronin's SLAM code requires the DEM to be a regular grid in UTM metres (an
Esri ASCII grid with a six-line header).  The 2020 GEOSLAM port kept the DEM
in geographic coordinates and approximated the cell size, which distorts the
geometry by several percent at mid-latitudes.  Here every DEM is resampled
onto a true UTM grid before any SLAM computation.

Readers:

* :func:`DEMGrid.from_arc_ascii` - Cronin's native format (UTM metres).
* :func:`read_geotiff` - GeoTIFF via ``rasterio`` if installed, otherwise via
  ``tifffile`` (enough for the OpenTopography / 3DEP / Copernicus products).

Projection: :meth:`GeographicDEM.to_utm` resamples a lon/lat raster onto a UTM
grid with bilinear interpolation (``scipy.ndimage.map_coordinates``).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .geodesy import UTMZone, lonlat_to_utm, utm_to_lonlat, utm_zone


@dataclass
class DEMGrid:
    """Elevation grid in UTM metres.

    ``z[i, j]`` is the elevation (metres, NaN = no data) of the cell whose
    centre is at ``x = x_min + (j + 0.5) * cellsize`` and
    ``y = y_max - (i + 0.5) * cellsize``.  Row 0 is the northern edge.
    """

    z: np.ndarray
    x_min: float
    y_max: float
    cellsize: float
    zone: UTMZone
    vertical_datum: str = "unknown"
    source: str = ""

    @property
    def shape(self):
        return self.z.shape

    @property
    def x_max(self):
        return self.x_min + self.z.shape[1] * self.cellsize

    @property
    def y_min(self):
        return self.y_max - self.z.shape[0] * self.cellsize

    @property
    def extent(self):
        """(x_min, x_max, y_min, y_max) for matplotlib ``imshow``."""
        return (self.x_min, self.x_max, self.y_min, self.y_max)

    def xy(self):
        """Cell-centre coordinate vectors (x along columns, y along rows)."""
        nr, nc = self.z.shape
        x = self.x_min + (np.arange(nc) + 0.5) * self.cellsize
        y = self.y_max - (np.arange(nr) + 0.5) * self.cellsize
        return x, y

    def lonlat_grid(self):
        x, y = self.xy()
        xx, yy = np.meshgrid(x, y)
        return utm_to_lonlat(xx, yy, self.zone)

    def elevation_at(self, x, y):
        """Bilinear elevation at UTM point(s)."""
        from scipy.ndimage import map_coordinates

        col = (np.asarray(x, float) - self.x_min) / self.cellsize - 0.5
        row = (self.y_max - np.asarray(y, float)) / self.cellsize - 0.5
        return map_coordinates(self.z, [np.atleast_1d(row), np.atleast_1d(col)], order=1,
                               mode="nearest", cval=np.nan)

    def crop(self, x_min, x_max, y_min, y_max) -> "DEMGrid":
        x, y = self.xy()
        cols = np.where((x >= x_min) & (x <= x_max))[0]
        rows = np.where((y >= y_min) & (y <= y_max))[0]
        if cols.size == 0 or rows.size == 0:
            raise ValueError("crop window does not overlap the DEM")
        z = self.z[rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1]
        return DEMGrid(z.copy(), self.x_min + cols[0] * self.cellsize,
                       self.y_max - rows[0] * self.cellsize, self.cellsize, self.zone,
                       self.vertical_datum, self.source)

    def thin(self, factor: int) -> "DEMGrid":
        """Decimate by an integer factor (Cronin's 'SLAM data diet') using block means."""
        f = int(factor)
        if f <= 1:
            return self
        nr, nc = (self.z.shape[0] // f) * f, (self.z.shape[1] // f) * f
        blocks = self.z[:nr, :nc].reshape(nr // f, f, nc // f, f)
        with np.errstate(invalid="ignore"):
            z = np.nanmean(blocks, axis=(1, 3))
        return DEMGrid(z, self.x_min, self.y_max, self.cellsize * f, self.zone,
                       self.vertical_datum, self.source)

    # ------------------------------------------------------------------ I/O
    @classmethod
    def from_arc_ascii(cls, path, zone: UTMZone, registration: str = "corner", **kw) -> "DEMGrid":
        """Read an Esri ASCII grid in UTM metres (Cronin's DEM format).

        ``registration='cronin'`` reproduces Cronin's node positions, which
        place node (i, j) at ``xllcorner + cellsize*(j-1)`` (a half-cell shift
        relative to the Esri cell-centre convention).
        """
        header = {}
        with open(path) as fh:
            for _ in range(6):
                pos = fh.tell()
                line = fh.readline().split()
                if not line or not line[0][0].isalpha():
                    fh.seek(pos)
                    break
                header[line[0].lower()] = float(line[1])
            data = np.loadtxt(fh)
        nodata = header.get("nodata_value")
        z = data.astype(float)
        if nodata is not None:
            z[z == nodata] = np.nan
        cs = header["cellsize"]
        nr = int(header["nrows"])
        if "xllcenter" in header:
            x_min, y_min = header["xllcenter"] - cs / 2, header["yllcenter"] - cs / 2
        else:
            x_min, y_min = header["xllcorner"], header["yllcorner"]
        if registration == "cronin":
            x_min, y_min = x_min - cs / 2, y_min - cs / 2
        return cls(z, x_min, y_min + nr * cs, cs, zone, **kw)

    def to_arc_ascii(self, path, nodata=-9999.0):
        z = np.where(np.isnan(self.z), nodata, self.z)
        with open(path, "w") as fh:
            fh.write(f"ncols {self.z.shape[1]}\nnrows {self.z.shape[0]}\n")
            fh.write(f"xllcorner {self.x_min:.3f}\nyllcorner {self.y_min:.3f}\n")
            fh.write(f"cellsize {self.cellsize:.6f}\nNODATA_value {nodata}\n")
            np.savetxt(fh, z, fmt="%.2f")


@dataclass
class GeographicDEM:
    """Elevation raster on a regular lon/lat grid (as served by most DEM APIs)."""

    z: np.ndarray
    lon_min: float
    lat_max: float
    dlon: float
    dlat: float
    vertical_datum: str = "unknown"
    source: str = ""

    @property
    def bounds(self):
        nr, nc = self.z.shape
        return (self.lon_min, self.lat_max - nr * self.dlat, self.lon_min + nc * self.dlon, self.lat_max)

    def to_utm(self, cellsize: float | None = None, zone: UTMZone | None = None,
               center: tuple[float, float] | None = None) -> DEMGrid:
        """Resample onto a UTM grid (bilinear).

        ``cellsize`` defaults to the native north-south spacing in metres.
        ``zone`` defaults to the zone of ``center`` (or of the raster centre),
        which should normally be the epicentre.
        """
        from scipy.ndimage import map_coordinates

        lo0, la0, lo1, la1 = self.bounds
        if center is None:
            center = ((lo0 + lo1) / 2, (la0 + la1) / 2)
        if zone is None:
            zone = utm_zone(*center)
        if cellsize is None:
            cellsize = float(np.round(self.dlat * 111_132.954, 1))
        # bounding box of the raster in UTM; keep only the fully covered core
        blon = np.r_[np.linspace(lo0, lo1, 50), np.full(50, lo1), np.linspace(lo1, lo0, 50), np.full(50, lo0)]
        blat = np.r_[np.full(50, la0), np.linspace(la0, la1, 50), np.full(50, la1), np.linspace(la1, la0, 50)]
        bx, by, _ = lonlat_to_utm(blon, blat, zone)
        # inner rectangle: max of left edge / min of right edge etc.
        x_min = np.max(bx[150:]) ; x_max = np.min(bx[50:100])
        y_min = np.max(by[:50]) ; y_max = np.min(by[100:150])
        x_min = np.ceil(x_min / cellsize) * cellsize
        y_max = np.floor(y_max / cellsize) * cellsize
        nc = int((x_max - x_min) // cellsize)
        nr = int((y_max - y_min) // cellsize)
        x = x_min + (np.arange(nc) + 0.5) * cellsize
        y = y_max - (np.arange(nr) + 0.5) * cellsize
        xx, yy = np.meshgrid(x, y)
        lon, lat = utm_to_lonlat(xx, yy, zone)
        col = (lon - self.lon_min) / self.dlon - 0.5
        row = (self.lat_max - lat) / self.dlat - 0.5
        zsrc = self.z.astype(float)
        mask = np.isnan(zsrc)
        if mask.any():
            fill = np.where(mask, np.nanmean(zsrc), zsrc)
            z = map_coordinates(fill, [row, col], order=1, mode="nearest")
            bad = map_coordinates(mask.astype(float), [row, col], order=1, mode="nearest") > 0.01
            z[bad] = np.nan
        else:
            z = map_coordinates(zsrc, [row, col], order=1, mode="nearest")
        return DEMGrid(z, float(x_min), float(y_max), float(cellsize), zone,
                       self.vertical_datum, self.source)


def _read_geotiff_tifffile(path):
    import tifffile

    with tifffile.TiffFile(path) as tif:
        page = tif.pages[0]
        z = page.asarray().astype(float)
        tags = page.tags
        scale = tags["ModelPixelScaleTag"].value
        tie = tags["ModelTiepointTag"].value
        geokeys = tags["GeoKeyDirectoryTag"].value if "GeoKeyDirectoryTag" in tags else ()
        nodata = tags["GDAL_NODATA"].value if "GDAL_NODATA" in tags else None
        area_or_point = 1
        epsg = None
        keys = list(geokeys)
        for k in range(4, len(keys), 4):
            key_id, loc, _cnt, val = keys[k:k + 4]
            if loc != 0:
                continue
            if key_id == 1025:
                area_or_point = val  # 1 = PixelIsArea, 2 = PixelIsPoint
            if key_id in (2048, 3072):
                epsg = val
    if nodata is not None:
        z[z == float(str(nodata).strip("\x00"))] = np.nan
    i, j, _, x, y, _ = tie[:6]
    dx, dy = scale[0], scale[1]
    x0 = x - i * dx
    y0 = y + j * dy
    if area_or_point == 2:  # tiepoint at pixel centre
        x0 -= dx / 2
        y0 += dy / 2
    return z, x0, y0, dx, dy, epsg


def read_geotiff(path, source=""):
    """Read a GeoTIFF; returns :class:`GeographicDEM` (EPSG:4326/4269) or :class:`DEMGrid` (UTM)."""
    try:
        import rasterio

        with rasterio.open(path) as ds:
            z = ds.read(1, masked=True).filled(np.nan).astype(float)
            t = ds.transform
            epsg = ds.crs.to_epsg() if ds.crs else None
            x0, y0, dx, dy = t.c, t.f, t.a, -t.e
    except ImportError:
        z, x0, y0, dx, dy, epsg = _read_geotiff_tifffile(path)
    z[z < -1e4] = np.nan
    if epsg in (4326, 4269, 4979, 4258) or (epsg is None and abs(x0) <= 360 and dx < 1):
        return GeographicDEM(z, x0, y0, dx, dy, source=source)
    if epsg and (32601 <= epsg <= 32660 or 32701 <= epsg <= 32760):
        zone = UTMZone(epsg % 100, epsg < 32700)
        if abs(dx - dy) > 1e-6:
            raise ValueError("non-square UTM cells are not supported; resample first")
        return DEMGrid(z, x0, y0, dx, zone, source=source)
    if epsg and 26901 <= epsg <= 26923:  # NAD83 / UTM north - treat as WGS84 (sub-metre)
        if abs(dx - dy) > 1e-6:
            raise ValueError("non-square UTM cells are not supported; resample first")
        return DEMGrid(z, x0, y0, dx, UTMZone(epsg - 26900, True), source=source)
    raise ValueError(f"unsupported CRS EPSG:{epsg}; reproject to WGS84 geographic or UTM")
