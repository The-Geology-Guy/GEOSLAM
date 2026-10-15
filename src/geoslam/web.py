"""JSON-in / JSON-out bridge used by the browser app (runs under Pyodide).

The web app fetches all data itself (ComCat, USGS 3DEP, faults) and hands
plain Python objects to these functions, so nothing here does network I/O.

``case`` objects exchanged with JavaScript::

    {"event": {"id", "lon", "lat", "depth_km", "mag", "mag_type", "time", "title", "url"},
     "uncertainty": null | {"semi_axes": [a1, a2, a3] (m), "axes": [[e1], [e2], [e3]],
                            "confidence": 0.95 | null, "source": str},
     "mechanisms": [{"label", "source", "product", "method",
                     "np1": {"strike", "dip", "rake", "su", "du", "ru"},
                     "np2": {...}}]}
"""
from __future__ import annotations

import json

import numpy as np

from .dem import DEMGrid, GeographicDEM
from .geodesy import UTMZone, lonlat_to_utm, utm_to_lonlat, utm_zone
from .orientation import NodalPlane
from .slam import Hypocenter, seismo_lineament, suggest_half_extent
from .sources.comcat import parse_detail
from .uncertainty import HypocenterUncertainty

_DEM: DEMGrid | None = None


# --------------------------------------------------------------------- helpers
def _plane_to_dict(p: NodalPlane) -> dict:
    return dict(strike=round(float(p.strike), 2), dip=round(float(p.dip), 2),
                rake=None if p.rake is None else round(float(p.rake), 2),
                su=float(p.strike_uncertainty), du=float(p.dip_uncertainty),
                ru=float(p.rake_uncertainty))


def _plane_from_dict(d: dict, label="") -> NodalPlane:
    return NodalPlane(float(d["strike"]), float(d["dip"]),
                      None if d.get("rake") is None else float(d["rake"]),
                      float(d.get("su") or 0), float(d.get("du") or 0), float(d.get("ru") or 0), label)


def _unc_to_dict(u: HypocenterUncertainty | None):
    if u is None:
        return None
    return dict(semi_axes=[float(a) for a in u.semi_axes], axes=u.axes.round(6).tolist(),
                confidence=u.confidence, source=u.source)


def _unc_from_dict(d):
    if not d:
        return None
    if "eh1" in d:  # Cronin-style manual entry (metres, true-north azimuth)
        return HypocenterUncertainty.from_cronin(float(d["eh1"]), float(d.get("eh1_az") or 0),
                                                 float(d["eh2"]), float(d["ez"]), source="entered")
    return HypocenterUncertainty(np.asarray(d["semi_axes"], float), np.asarray(d["axes"], float),
                                 d.get("confidence"), d.get("source", ""))


def describe_uncertainty(d) -> dict:
    """Horizontal/vertical half-widths in plain terms, for display."""
    if isinstance(d, str):
        d = json.loads(d)
    u = _unc_from_dict(d)
    if u is None:
        return json.dumps({"text": "No location uncertainty reported"})
    from .orientation import plane_normal
    h = float(u.half_width(plane_normal(90, 0)))
    h2 = float(u.half_width(plane_normal(90, 90)))
    v = float(u.half_width(plane_normal(0, 0)))
    return json.dumps({"horizontal_m": [round(min(h, h2)), round(max(h, h2))], "vertical_m": round(v),
                       "confidence": u.confidence, "source": u.source})


# ---------------------------------------------------------------- event intake
def case_from_comcat(detail_json: str, quakeml: dict) -> str:
    """Parse a ComCat GeoJSON detail plus a {url: xml_text} map of product QuakeML."""
    detail = json.loads(detail_json)
    if detail.get("type") == "FeatureCollection":
        detail = detail["features"][0]
    if isinstance(quakeml, str):
        quakeml = json.loads(quakeml)

    def fetch(url):
        if url in quakeml and quakeml[url]:
            return quakeml[url].encode() if isinstance(quakeml[url], str) else quakeml[url]
        raise OSError("not fetched")

    eq = parse_detail(detail, fetch=fetch)
    mechs = []
    for m in eq.mechanisms:
        has_unc = m.np1.strike_uncertainty > 0 or m.np1.dip_uncertainty > 0
        mechs.append(dict(
            label=f"{m.source.upper()} {m.method or m.product}".strip(),
            source=m.source, product=m.product, method=m.method,
            has_uncertainty=bool(has_unc), weight=m.extra.get("weight"),
            np1=_plane_to_dict(m.np1), np2=_plane_to_dict(m.np2)))
    # first-motion mechanisms with uncertainties first, then by product weight
    mechs.sort(key=lambda m: (not m["has_uncertainty"], m["product"] != "focal-mechanism",
                              -(m["weight"] or 0)))
    props = detail["properties"]
    case = dict(
        event=dict(id=eq.event_id, lon=eq.lon, lat=eq.lat, depth_km=eq.depth_km,
                   mag=eq.magnitude, mag_type=eq.magnitude_type, time=eq.time,
                   title=props.get("title", ""), url=eq.url),
        uncertainty=_unc_to_dict(eq.uncertainty), mechanisms=mechs, notes=eq.notes)
    return json.dumps(case)


def auxiliary(plane_json: str) -> str:
    p = _plane_from_dict(json.loads(plane_json))
    a = p.auxiliary()
    a.strike_uncertainty, a.dip_uncertainty, a.rake_uncertainty = (
        p.strike_uncertainty, p.dip_uncertainty, p.rake_uncertainty)
    return json.dumps(_plane_to_dict(a))


# ------------------------------------------------------------------- DEM setup
def dem_plan(case_json: str, half_extent_km=None, cellsize=None, max_cells=800,
             min_half_km=15.0, max_half_km=40.0) -> str:
    """Choose map extent and cell size; return the 3DEP request in UTM."""
    case = json.loads(case_json)
    ev = case["event"]
    hyp = Hypocenter(ev["lon"], ev["lat"], ev["depth_km"], _unc_from_dict(case.get("uncertainty")))
    if half_extent_km:
        half = float(half_extent_km) * 1000
    else:
        half = 0.0
        for m in case["mechanisms"][:1] or []:
            for k in ("np1", "np2"):
                half = max(half, suggest_half_extent(hyp, _plane_from_dict(m[k])))
        half = min(max(half, min_half_km * 1000), max_half_km * 1000)
    if not cellsize:
        cellsize = max(10.0, np.ceil(2 * half / max_cells / 10.0) * 10.0)
    cellsize = float(cellsize)
    n = int(np.ceil(2 * half / cellsize))
    half = n * cellsize / 2
    zone = utm_zone(ev["lon"], ev["lat"])
    e, nn, _ = lonlat_to_utm(ev["lon"], ev["lat"], zone)
    e = round(float(e) / cellsize) * cellsize
    nn = round(float(nn) / cellsize) * cellsize
    box = [e - half, nn - half, e + half, nn + half]
    corners = _corners(zone, box)
    lo, la = utm_to_lonlat(np.array([box[0], box[2], box[2], box[0]]),
                           np.array([box[1], box[1], box[3], box[3]]), zone)
    return json.dumps(dict(epsg=zone.epsg, zone=zone.number, north=zone.north, bbox=box,
                           cellsize=cellsize, size=[n, n], corners=corners,
                           lonlat_bounds=[float(lo.min()), float(la.min()), float(lo.max()), float(la.max())]))


def _corners(zone, box):
    """Image corners (lon, lat) in MapLibre order: TL, TR, BR, BL."""
    x0, y0, x1, y1 = box
    lo, la = utm_to_lonlat(np.array([x0, x1, x1, x0]), np.array([y1, y1, y0, y0]), zone)
    return [[round(float(a), 7), round(float(b), 7)] for a, b in zip(lo, la)]


def set_dem_utm(values, nrows, ncols, x_min, y_max, cellsize, epsg, source="USGS 3DEP"):
    """Install a UTM DEM from a flat float32 buffer (row 0 = north)."""
    global _DEM
    z = np.asarray(values, dtype=np.float32).astype(float)[: nrows * ncols].reshape(nrows, ncols)
    z[(z < -1e4) | (z > 1e5)] = np.nan
    epsg = int(epsg)
    if 26901 <= epsg <= 26923:
        zone = UTMZone(epsg - 26900, True)
    else:
        zone = UTMZone(epsg % 100, epsg < 32700)
    _DEM = DEMGrid(z, float(x_min), float(y_max), float(cellsize), zone, source=source)
    return _dem_info()


def set_dem_geographic(values, nrows, ncols, lon_min, lat_max, dlon, dlat, center_lon, center_lat,
                       cellsize=None, box=None, source="uploaded GeoTIFF"):
    """Install a lon/lat DEM (e.g. an uploaded GeoTIFF) by resampling to UTM."""
    global _DEM
    z = np.asarray(values, dtype=np.float32).astype(float)[: nrows * ncols].reshape(nrows, ncols)
    z[(z < -1e4) | (z > 1e5)] = np.nan
    g = GeographicDEM(z, float(lon_min), float(lat_max), float(dlon), float(dlat), source=source)
    zone = utm_zone(center_lon, center_lat)
    d = g.to_utm(cellsize=cellsize, zone=zone)
    if box:
        try:
            d = d.crop(box[0], box[2], box[1], box[3])
        except ValueError:
            pass
    _DEM = d
    return _dem_info()


def set_dem_ascii(text, zone_number, north=True, source="uploaded ASCII grid"):
    global _DEM
    import io
    import tempfile
    import os
    with tempfile.NamedTemporaryFile("w", suffix=".asc", delete=False) as fh:
        fh.write(text)
        path = fh.name
    try:
        _DEM = DEMGrid.from_arc_ascii(path, UTMZone(int(zone_number), bool(north)), source=source)
    finally:
        os.unlink(path)
    return _dem_info()


def crop_dem(box_json):
    """Crop the current DEM to [x_min, y_min, x_max, y_max] when they overlap."""
    global _DEM
    x0, y0, x1, y1 = json.loads(box_json)
    try:
        _DEM = _DEM.crop(x0, x1, y0, y1)
    except ValueError:
        pass
    return _dem_info()


def thin_dem(max_cells=800):
    global _DEM
    f = int(np.ceil(max(_DEM.shape) / max_cells))
    if f > 1:
        _DEM = _DEM.thin(f)
    return _dem_info()


def _dem_info():
    d = _DEM
    finite = np.isfinite(d.z)
    return json.dumps(dict(shape=list(d.shape), cellsize=d.cellsize, epsg=d.zone.epsg,
                           box=[d.x_min, d.y_min, d.x_max, d.y_max],
                           corners=_corners(d.zone, [d.x_min, d.y_min, d.x_max, d.y_max]),
                           valid_fraction=float(finite.mean()),
                           elev_range=[float(np.nanmin(d.z)), float(np.nanmax(d.z))] if finite.any() else None,
                           source=d.source))


def hillshade(azimuth=315.0, altitude=15.0, exaggeration=2.0):
    """Grey hillshade (RGBA bytes) of the current DEM, lit Cronin's way."""
    d = _DEM
    z = d.z.copy()
    nan = np.isnan(z)
    if nan.any():
        z[nan] = np.nanmean(z) if (~nan).any() else 0.0
    g_row, g_col = np.gradient(z * exaggeration, d.cellsize)
    dzdx, dzdy = g_col, -g_row                            # rows run north -> south
    az, alt = np.radians(azimuth), np.radians(altitude)
    light = np.array([np.sin(az) * np.cos(alt), np.cos(az) * np.cos(alt), np.sin(alt)])
    shade = (-dzdx * light[0] - dzdy * light[1] + light[2]) / np.sqrt(1 + dzdx**2 + dzdy**2)
    shade = np.clip(shade, 0, 1)
    lo, hi = np.percentile(shade[~nan], [1, 99]) if (~nan).any() else (0, 1)
    g = np.clip((shade - lo) / max(hi - lo, 1e-6), 0, 1)
    rgba = np.empty(z.shape + (4,), np.uint8)
    rgba[..., 0] = rgba[..., 1] = rgba[..., 2] = (40 + 205 * g).astype(np.uint8)
    rgba[..., 3] = np.where(nan, 0, 255)
    return rgba.tobytes(), d.shape[1], d.shape[0]


def beachball(strike, dip, rake, n=161):
    from .plotting import beachball_image
    img = beachball_image(float(strike), float(dip), float(rake), n)
    return (img * 255).astype(np.uint8).tobytes(), n


# ------------------------------------------------------------------------ run
def _lines_lonlat(lines, zone, step=1):
    out = []
    for l in lines:
        if len(l) < 2:
            continue
        lo, la = utm_to_lonlat(l[::step, 0], l[::step, 1], zone)
        out.append(np.column_stack([lo, la]).round(6).tolist())
    return out


def _polygons(result):
    """Filled seismo-lineament as MultiPolygon coordinates (lon/lat)."""
    from contourpy import contour_generator

    x, y = result.dem.xy()
    m = result.margin if result.margin is not None else np.where(result.swath, -1.0, 1.0)
    m = np.where(np.isnan(m), 1e9, m)
    gen = contour_generator(x, y, m, fill_type="OuterOffset")
    pts_list, offs_list = gen.filled(-1e15, 0.0)
    polys = []
    for pts, offs in zip(pts_list, offs_list):
        rings = []
        for a, b in zip(offs[:-1], offs[1:]):
            ring = pts[a:b]
            lo, la = utm_to_lonlat(ring[:, 0], ring[:, 1], result.dem.zone)
            rings.append(np.column_stack([lo, la]).round(6).tolist())
        polys.append(rings)
    return polys


def run(case_json: str, options_json: str) -> str:
    """Run SLAM for the selected planes; returns GeoJSON layers and summaries."""
    if _DEM is None:
        raise RuntimeError("no DEM loaded")
    case, opt = json.loads(case_json), json.loads(options_json)
    ev = case["event"]
    unc = _unc_from_dict(case.get("uncertainty"))
    unknown = None
    if unc is None and opt.get("assumed_location_m"):
        eh, ez = opt["assumed_location_m"]
        unknown = HypocenterUncertainty.isotropic_horizontal(float(eh), float(ez), source="assumed")
    hyp = Hypocenter(ev["lon"], ev["lat"], ev["depth_km"], unc,
                     opt.get("depth_datum", "sea_level"), ev.get("id", ""))
    out = dict(planes=[])
    from contourpy import contour_generator

    for item in opt["planes"]:
        plane = _plane_from_dict(item["plane"], item.get("label", ""))
        r = seismo_lineament(_DEM, hyp, plane, method=opt.get("method", "envelope"),
                             multiplier=float(opt.get("multiplier", 1.0)),
                             unknown_uncertainty=unknown, chunk=int(opt.get("chunk", 20000)),
                             angular_step=float(opt.get("angular_step", 1.0)))
        x, y = _DEM.xy()
        traces = {}
        for name, sd in [("Mean", r.mean_signed_distance)] + [
                (k, v) for k, v in r.signed_distances.items() if k != "Mean"]:
            lines = contour_generator(x, y, np.nan_to_num(sd)).lines(0.0)
            traces[name] = _lines_lonlat(lines, _DEM.zone)
        out["planes"].append(dict(
            label=item.get("label", ""), strike=plane.strike, dip=plane.dip, rake=plane.rake,
            su=plane.strike_uncertainty, du=plane.dip_uncertainty,
            area_km2=round(r.area_km2, 2), gamma=round(r.gamma, 4),
            grid_dip_azimuth=round(r.params["grid_dip_azimuth"], 3),
            light_azimuth=round((r.params["grid_dip_azimuth"] + 180) % 360, 2),
            focus=[round(float(v), 1) for v in r.focus_xyz],
            planes=[dict(name=p.name, dip=round(p.dip, 2), dip_azimuth=round(p.dip_azimuth % 360, 2),
                         half_width_m=round(p.half_width, 1)) for p in r.planes],
            polygons=_polygons(r), traces=traces, summary=r.summary(),
            location_uncertainty=(unc or unknown).source if (unc or unknown) else "none"))
    return json.dumps(out)


# --------------------------------------------------------------------- report
_DEM_SLOTS: dict[str, DEMGrid] = {}


def save_dem(slot: str) -> bool:
    """Keep the current DEM under a name (so a report DEM can be fetched and the
    analysis DEM restored afterwards without re-downloading it)."""
    if _DEM is None:
        return False
    _DEM_SLOTS[slot] = _DEM
    return True


def load_dem(slot: str) -> str:
    global _DEM
    if slot not in _DEM_SLOTS:
        return ""
    _DEM = _DEM_SLOTS[slot]
    return _dem_info()


def report_plan(cases_json: str, half_extent_km=None, cellsize=None, max_cells=1200,
                margin_km=None) -> str:
    """One map frame covering every event in a report.

    Each event gets the square ``dem_plan`` would give it (or ``half_extent_km``);
    the frame is the union of those squares in the UTM zone of their centroid.
    """
    cases = json.loads(cases_json)
    if not cases:
        raise ValueError("no events in the report")
    lons = [c["event"]["lon"] for c in cases]
    lats = [c["event"]["lat"] for c in cases]
    zone = utm_zone(float(np.mean(lons)), float(np.mean(lats)))
    boxes = []
    for c in cases:
        p = json.loads(dem_plan(json.dumps(c), half_extent_km, 100.0))
        half = (p["bbox"][2] - p["bbox"][0]) / 2
        if margin_km is not None and len(cases) > 1:
            half = float(margin_km) * 1000
        e, n, _ = lonlat_to_utm(c["event"]["lon"], c["event"]["lat"], zone)
        boxes.append([float(e) - half, float(n) - half, float(e) + half, float(n) + half])
    b = np.array(boxes)
    x0, y0, x1, y1 = b[:, 0].min(), b[:, 1].min(), b[:, 2].max(), b[:, 3].max()
    if not cellsize:
        cellsize = max(10.0, np.ceil(max(x1 - x0, y1 - y0) / max_cells / 10.0) * 10.0)
    cellsize = float(cellsize)
    x0, y0 = np.floor(x0 / cellsize) * cellsize, np.floor(y0 / cellsize) * cellsize
    x1, y1 = np.ceil(x1 / cellsize) * cellsize, np.ceil(y1 / cellsize) * cellsize
    box = [float(x0), float(y0), float(x1), float(y1)]
    lo, la = utm_to_lonlat(np.array([x0, x1, x1, x0]), np.array([y0, y0, y1, y1]), zone)
    return json.dumps(dict(epsg=zone.epsg, zone=zone.number, north=zone.north, bbox=box,
                           cellsize=cellsize,
                           size=[int(round((x1 - x0) / cellsize)), int(round((y1 - y0) / cellsize))],
                           corners=_corners(zone, box),
                           lonlat_bounds=[float(lo.min()), float(la.min()), float(lo.max()), float(la.max())]))


def _report_result(item, dem):
    case = item["case"]
    ev = case["event"]
    unc = _unc_from_dict(case.get("uncertainty"))
    unknown = None
    if unc is None and item.get("assumed_location_m"):
        eh, ez = item["assumed_location_m"]
        unknown = HypocenterUncertainty.isotropic_horizontal(float(eh), float(ez), source="assumed")
    hyp = Hypocenter(ev["lon"], ev["lat"], ev["depth_km"], unc, item.get("depth_datum", "sea_level"),
                     ev.get("id", ""), ev.get("time", "") or "", ev.get("mag"))
    plane = _plane_from_dict(item["plane"], item.get("label", ""))
    r = seismo_lineament(dem, hyp, plane, method=item.get("method", "envelope"),
                         multiplier=float(item.get("multiplier", 1.0)), unknown_uncertainty=unknown,
                         chunk=20000, angular_step=float(item.get("angular_step", 1.0)))
    u = unc or unknown
    loc = ""
    if u is not None:
        a = sorted(float(v) for v in u.semi_axes)
        loc = f"{a[-1]:.0f}/{a[1]:.0f}/{a[0]:.0f} m" + (" (assumed)" if unc is None else "")
    details = dict(date=(ev.get("time") or "")[:19].replace("T", " "),
                   mag=f"{ev['mag']:.1f}" if ev.get("mag") is not None else "",
                   location_unc=loc, source=item.get("source", ""))
    return r, details


def report_render(items_json: str, options_json: str) -> bytes:
    """Render a Cronin-style report figure of the current DEM; returns file bytes.

    ``items``: [{label, event_label, case, plane, method, multiplier, depth_datum,
    assumed_location_m, source}].  ``options``: ReportOptions fields plus
    ``faults`` (list of lon/lat polylines) and ``quakes`` ([[lon, lat, mag], ...]).
    """
    from .report import ReportItem, ReportOptions, render

    if _DEM is None:
        raise RuntimeError("no DEM loaded")
    items, opt = json.loads(items_json), json.loads(options_json)
    faults = [{"lines": [np.asarray(l, float)], "properties": {}} for l in (opt.pop("faults", None) or [])
              if len(l) > 1]
    quakes = opt.pop("quakes", None)
    quakes = np.asarray(quakes, float) if quakes else None
    known = ReportOptions.__dataclass_fields__
    ro = ReportOptions(**{k: v for k, v in opt.items() if k in known})
    ritems = []
    for it in items:
        r, details = _cached_result(it, _DEM)
        ritems.append(ReportItem(it["label"], it["event_label"], r, it.get("color"), details))
    return render(_DEM, ritems, ro, faults=faults, quakes=quakes)


def terrain_png(style="hillshade", azimuth=315.0, altitude=15.0, exaggeration=2.0) -> bytes:
    """The current DEM as a terrain image (one pixel per cell, north up)."""
    from .report import terrain_image
    return terrain_image(_DEM, style, float(azimuth), float(altitude), float(exaggeration))


def world_file() -> str:
    """ESRI world file (.pgw) for :func:`terrain_png` (pixel centres, UTM metres)."""
    d = _DEM
    c = d.cellsize
    return f"{c:.6f}\n0.0\n0.0\n{-c:.6f}\n{d.x_min + c / 2:.3f}\n{d.y_max - c / 2:.3f}\n"


def projection_wkt() -> str:
    """ESRI .prj text for the current DEM's WGS84 / UTM zone."""
    z = _DEM.zone
    name = f"WGS_1984_UTM_Zone_{z.number}{'N' if z.north else 'S'}"
    return (f'PROJCS["{name}",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,'
            f'298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],'
            f'PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],'
            f'PARAMETER["False_Northing",{0.0 if z.north else 10000000.0}],'
            f'PARAMETER["Central_Meridian",{z.central_meridian}],PARAMETER["Scale_Factor",0.9996],'
            f'PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]')


def dem_ascii() -> str:
    """The current DEM as an Esri ASCII grid (Cronin's SLAM input format)."""
    import os
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".asc", delete=False) as fh:
        path = fh.name
    try:
        _DEM.to_arc_ascii(path)
        with open(path) as fh:
            return fh.read()
    finally:
        os.unlink(path)


_RESULT_CACHE: dict = {}


def _cached_result(item, dem):
    key = (id(dem), json.dumps({k: v for k, v in item.items() if k not in ("label", "event_label", "color")},
                               sort_keys=True))
    if key not in _RESULT_CACHE:
        if len(_RESULT_CACHE) > 24:
            _RESULT_CACHE.clear()
        _RESULT_CACHE[key] = _report_result(item, dem)
    return _RESULT_CACHE[key]


def _zip(files: dict) -> bytes:
    import io
    import zipfile
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name, data in files.items():
            z.writestr(name, data)
    return buf.getvalue()


def terrain_zip(name="terrain", style="hillshade", azimuth=315.0, altitude=15.0, exaggeration=2.0) -> bytes:
    """Terrain PNG with its world file (.pgw) and projection (.prj), ready for GIS."""
    return _zip({f"{name}.png": terrain_png(style, azimuth, altitude, exaggeration),
                 f"{name}.pgw": world_file(), f"{name}.prj": projection_wkt(),
                 "README.txt": _readme("terrain image", style)})


def dem_zip(name="dem") -> bytes:
    """The DEM as an Esri ASCII grid (.asc + .prj), the format Cronin's SLAMcode reads."""
    return _zip({f"{name}.asc": dem_ascii(), f"{name}.prj": projection_wkt(),
                 "README.txt": _readme("elevation grid")})


def _readme(what, style=""):
    d = _DEM
    z = d.zone
    return (f"GEOSLAM {what}{' (' + style + ')' if style else ''}\n"
            f"Source: {d.source}\nProjection: WGS 84 / UTM zone {z.number}{'N' if z.north else 'S'} "
            f"(EPSG:{z.epsg})\nCell size: {d.cellsize:g} m; {d.shape[1]} columns x {d.shape[0]} rows\n"
            f"Extent (m): x {d.x_min:.1f} to {d.x_max:.1f}, y {d.y_min:.1f} to {d.y_max:.1f}\n"
            "Elevations in metres; NODATA -9999.\nMade with GEOSLAM, https://geoslam.org\n")


def report_light_azimuth(item_json: str) -> float:
    """Cronin's light direction (opposite the grid dip direction) for a report item."""
    r, _ = _cached_result(json.loads(item_json), _DEM)
    return float((r.params["grid_dip_azimuth"] + 180.0) % 360.0)
