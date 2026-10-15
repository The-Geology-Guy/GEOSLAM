"""Publication-style SLAM figures in the manner of Cronin's SLAM maps.

A report figure shows one terrain frame (UTM, so distances and angles are
true) with longitude/latitude labels on all four edges, and for every
selected nodal plane:

* the seismo-lineament (light band with dashed boundaries),
* the ground trace of the mean nodal plane (solid black),
* optionally the traces of the bounding planes,

plus mapped Quaternary faults, seismicity, and a labelled beachball for each
event (``A1`` = event A, nodal plane 1).  Output is PNG, PDF or SVG; PDF and
SVG keep the linework and beachballs as vectors.

Only matplotlib is needed (it is available in Pyodide).
"""
from __future__ import annotations

import io
from dataclasses import dataclass, field

import numpy as np

from .dem import DEMGrid
from .geodesy import lonlat_to_utm, utm_to_lonlat
from .orientation import plane_normal, slip_vector

# Cronin-style defaults; further colours are used when several events are mapped.
SWATH_COLORS = ["#9fd8ec", "#f2b26b", "#b9a3e3", "#8fd19e", "#f08f8f", "#e3d06b", "#7fc4c4"]
FAULT_COLOR = "#e6d43a"
QUAKE_COLOR = "#ff9f1c"

PAGES = {  # inches (width, height) of the whole figure
    "map": None,              # sized to the map's aspect ratio
    "letter-portrait": (8.5, 11.0),
    "letter-landscape": (11.0, 8.5),
    "a4-portrait": (8.27, 11.69),
    "a4-landscape": (11.69, 8.27),
    "square": (9.0, 9.0),
}


@dataclass
class ReportItem:
    """One nodal plane to draw.  ``result`` is a :class:`~geoslam.slam.SLAMResult`."""

    label: str                      # e.g. "A1"
    event_label: str                # e.g. "A"
    result: object
    color: str | None = None
    details: dict = field(default_factory=dict)


@dataclass
class ReportOptions:
    terrain: str = "hillshade"      # hillshade | hillshade-color | elevation | slope | none
    light: str = "plane"            # plane (Cronin) | nw | custom
    light_azimuth: float = 315.0
    light_altitude: float = 15.0
    exaggeration: float = 2.0
    swath_fill: bool = True
    swath_edges: bool = True
    bounding_traces: bool = False
    mean_trace: bool = True
    faults: bool = True
    seismicity: bool = True
    beachballs: bool = True
    labels: bool = True
    graticule: bool = True
    scale_bar: bool = True
    north_arrow: bool = True
    legend: bool = False
    details: bool = False
    title: str = ""
    color_by_event: bool | None = None   # None = automatic (only when >1 event)
    page: str = "map"
    dpi: int = 300
    fmt: str = "png"                # png | pdf | svg
    beachball_size: float = 0.75    # inches
    beachball_offset: bool = True   # place beside the epicentre with a leader line
    width_in: float = 8.0           # map width when page == "map"


# ------------------------------------------------------------------ terrain
def shade(dem: DEMGrid, azimuth: float, altitude: float, exaggeration: float = 2.0):
    """Lambertian hillshade in [0, 1] (NaN where the DEM has no data)."""
    z = dem.z.astype(float)
    nan = np.isnan(z)
    if nan.all():
        return np.full(z.shape, np.nan)
    zf = np.where(nan, np.nanmean(z), z)
    g_row, g_col = np.gradient(zf * exaggeration, dem.cellsize)
    dzdx, dzdy = g_col, -g_row
    az, alt = np.radians(azimuth), np.radians(altitude)
    lx, ly, lz = np.sin(az) * np.cos(alt), np.cos(az) * np.cos(alt), np.sin(alt)
    s = (-dzdx * lx - dzdy * ly + lz) / np.sqrt(1 + dzdx**2 + dzdy**2)
    s = np.clip(s, 0, 1)
    lo, hi = np.nanpercentile(s[~nan], [1, 99.5])
    s = np.clip((s - lo) / max(hi - lo, 1e-9), 0, 1)
    s[nan] = np.nan
    return s


def terrain_rgb(dem: DEMGrid, style: str, azimuth: float, altitude: float, exaggeration=2.0):
    """RGB(A) float image (rows north->south) for a terrain style."""
    import matplotlib
    from matplotlib.colors import LinearSegmentedColormap

    nan = np.isnan(dem.z)
    if style == "none":
        img = np.ones(dem.shape + (4,))
        return img
    if style == "slope":
        z = np.where(nan, np.nanmean(dem.z), dem.z)
        gy, gx = np.gradient(z, dem.cellsize)
        sl = np.degrees(np.arctan(np.hypot(gx, gy)))
        v = np.clip(sl / max(np.nanpercentile(sl, 99), 1e-6), 0, 1)
        img = matplotlib.colormaps["Greys"](v)
    else:
        hs = shade(dem, azimuth, altitude, exaggeration)
        if style == "hillshade":
            g = 0.18 + 0.72 * np.nan_to_num(hs, nan=0.5)      # dark-grey Cronin look
            img = np.dstack([g, g, g, np.ones_like(g)])
        else:
            terrain = LinearSegmentedColormap.from_list(
                "slam_terrain", ["#4f7a5a", "#93ad78", "#d8cf9a", "#c7a274", "#9c7a5c", "#f2efe9"])
            z = dem.z
            lo, hi = np.nanpercentile(z, [1, 99])
            c = terrain(np.clip((np.nan_to_num(z, nan=lo) - lo) / max(hi - lo, 1e-6), 0, 1))
            if style == "elevation":
                img = c
            else:  # hillshade-color: soft-light blend
                h = np.nan_to_num(hs, nan=0.5)[..., None]
                rgb = c[..., :3]
                img = np.concatenate([np.clip(rgb * (0.45 + 0.75 * h), 0, 1), np.ones(h.shape)], -1)
    img = np.array(img, float)
    img[nan, 3] = 0.0
    return img


# ---------------------------------------------------------------- beachball
def _beachball(ax, cx, cy, radius, strike, dip, rake, color="k", zorder=30, lw=1.2):
    """Vector lower-hemisphere beachball centred at (cx, cy) in data units."""
    from matplotlib.patches import Circle

    n = 181
    u = np.linspace(-1, 1, n)
    X, Y = np.meshgrid(u, u)          # Y up = north
    R = np.hypot(X, Y)
    inc = 2 * np.arcsin(np.clip(R / np.sqrt(2), 0, 1))
    az = np.arctan2(X, Y)
    ray = np.stack([np.sin(inc) * np.sin(az), np.sin(inc) * np.cos(az), -np.cos(inc)], -1)
    pol = (ray @ plane_normal(dip, (strike + 90) % 360)) * (ray @ slip_vector(strike, dip, rake))
    pol = np.where(R <= 1.0, pol, np.nan)
    ax.add_patch(Circle((cx, cy), radius, facecolor="white", edgecolor="none", zorder=zorder))
    ax.contourf(cx + X * radius, cy + Y * radius, pol, levels=[0, np.inf], colors=[color], zorder=zorder + 0.1,
                antialiased=True)
    ax.add_patch(Circle((cx, cy), radius, facecolor="none", edgecolor=color, lw=lw, zorder=zorder + 0.2))


# ---------------------------------------------------------------- graticule
def _nice_step(span_deg):
    for s in (0.005, 0.01, 0.02, 0.025, 0.05, 0.1, 0.2, 0.25, 0.5, 1.0, 2.0):
        if span_deg / s <= 7:
            return s
    return 5.0


def _fmt(v, pos_hemi, neg_hemi, step):
    dec = max(0, int(np.ceil(-np.log10(step) - 1e-9)))
    return f"{abs(v):.{dec}f}°{pos_hemi if v >= 0 else neg_hemi}"


def _edge_ticks(dem, kind, edge, n=400):
    """Positions along a frame edge where lon (kind='lon') or lat is a round value."""
    x0, x1, y0, y1 = dem.x_min, dem.x_max, dem.y_min, dem.y_max
    if edge in ("bottom", "top"):
        xs = np.linspace(x0, x1, n)
        ys = np.full(n, y0 if edge == "bottom" else y1)
    else:
        ys = np.linspace(y0, y1, n)
        xs = np.full(n, x0 if edge == "left" else x1)
    lon, lat = utm_to_lonlat(xs, ys, dem.zone)
    val = lon if kind == "lon" else lat
    step = _nice_step(abs(val.max() - val.min()) or 1e-3)
    ticks = []
    for t in np.arange(np.ceil(val.min() / step) * step, val.max() + 1e-12, step):
        k = np.searchsorted(val if val[-1] > val[0] else val[::-1], t)
        if val[-1] < val[0]:
            k = n - k
        if 0 < k < n:
            a, b = k - 1, k
            f = (t - val[a]) / (val[b] - val[a])
            pos = (xs if edge in ("bottom", "top") else ys)[a] + f * (
                (xs if edge in ("bottom", "top") else ys)[b] - (xs if edge in ("bottom", "top") else ys)[a])
            ticks.append((pos, t))
    return ticks, step


def _graticule(ax, dem, fontsize):
    for edge in ("bottom", "top", "left", "right"):
        kind = "lon" if edge in ("bottom", "top") else "lat"
        ticks, step = _edge_ticks(dem, kind, edge)
        pos = [p for p, _ in ticks]
        labs = [_fmt(t, "E", "W", step) if kind == "lon" else _fmt(t, "N", "S", step) for _, t in ticks]
        if edge == "bottom":
            ax.set_xticks(pos, labs)
        elif edge == "left":
            ax.set_yticks(pos, labs)
        elif edge == "top":
            sec = ax.secondary_xaxis("top")
            sec.set_xticks(pos, labs)
            sec.tick_params(labelsize=fontsize, length=4, width=0.8, direction="out")
        else:
            sec = ax.secondary_yaxis("right")
            sec.set_yticks(pos, labs)
            sec.tick_params(labelsize=fontsize, length=4, width=0.8, direction="out")
    ax.tick_params(labelsize=fontsize, length=4, width=0.8, direction="out")


def _scale_bar(ax, dem, fontsize):
    import matplotlib.patheffects as pe

    width = dem.x_max - dem.x_min
    target = width / 5
    nice = [100, 200, 250, 500, 1000, 2000, 2500, 5000, 10000, 20000, 25000, 50000]
    L = max([v for v in nice if v <= target] or [nice[0]])
    x0 = dem.x_min + 0.04 * width
    y0 = dem.y_min + 0.04 * (dem.y_max - dem.y_min)
    h = 0.008 * (dem.y_max - dem.y_min)
    halo = [pe.withStroke(linewidth=3, foreground="white")]
    for i in range(4):
        ax.add_patch(__import__("matplotlib.patches", fromlist=["Rectangle"]).Rectangle(
            (x0 + i * L / 4, y0), L / 4, h, facecolor="k" if i % 2 == 0 else "white",
            edgecolor="k", lw=0.8, zorder=40))
    lab = f"{L / 1000:g} km" if L >= 1000 else f"{L:g} m"
    ax.text(x0 + L / 2, y0 + 2.2 * h, lab, ha="center", va="bottom", fontsize=fontsize,
            zorder=41, path_effects=halo)


def _north_arrow(ax, dem, gamma, fontsize):
    """Arrow pointing to true north (grid north rotated by -gamma)."""
    import matplotlib.patheffects as pe

    W, H = dem.x_max - dem.x_min, dem.y_max - dem.y_min
    cx, cy = dem.x_max - 0.06 * W, dem.y_min + 0.08 * H
    L = 0.06 * H
    g = np.radians(-gamma)
    dx, dy = L * np.sin(g), L * np.cos(g)
    ax.annotate("", xy=(cx + dx, cy + dy), xytext=(cx, cy),
                arrowprops=dict(arrowstyle="-|>", lw=1.6, color="k", mutation_scale=14,
                                path_effects=[pe.withStroke(linewidth=4, foreground="white")]), zorder=42)
    ax.text(cx + 1.25 * dx, cy + 1.25 * dy, "N", ha="center", va="bottom", fontsize=fontsize + 1,
            fontweight="bold", zorder=42, path_effects=[pe.withStroke(linewidth=3, foreground="white")])


def _beachball_spots(items, dem, radius, offset=True):
    """One beachball per event, beside its epicentre (off the mean traces), nudged
    apart if they overlap and kept inside the frame."""
    events = {}
    for it in items:
        ev = events.setdefault(it.event_label, dict(items=[], x=None, y=None))
        ev["items"].append(it)
        ev["x"], ev["y"] = it.result.focus_xyz[:2]
    pts = [(k, np.array([v["x"], v["y"]], float)) for k, v in events.items()]
    offs = {k: np.array([0.0, 0.0]) for k, _ in pts}
    if offset:
        for k, p in pts:
            # directions of the mean traces through the epicentre (grid strike)
            strikes = [np.radians(it.result.params["grid_dip_azimuth"] - 90.0) for it in events[k]["items"]]
            best, best_score = None, -np.inf
            for ang in np.radians(np.arange(0, 360, 15)):
                d = np.array([np.sin(ang), np.cos(ang)])
                # angular clearance from every trace (both directions)
                clear = min(abs(np.sin(ang - s)) for s in strikes)
                q = p + 2.6 * radius * d
                inside = (dem.x_min + 1.3 * radius < q[0] < dem.x_max - 1.3 * radius and
                          dem.y_min + 1.3 * radius < q[1] < dem.y_max - 2.2 * radius)
                score = clear + (2.0 if inside else 0.0) + 0.05 * np.cos(ang - np.radians(45))
                if score > best_score:
                    best, best_score = d, score
            offs[k] = 2.6 * radius * best
    for _ in range(60):
        moved = False
        for i in range(len(pts)):
            for j in range(i + 1, len(pts)):
                ki, pi = pts[i]
                kj, pj = pts[j]
                d = (pi + offs[ki]) - (pj + offs[kj])
                dist = np.hypot(*d)
                if dist < 2.4 * radius:
                    push = (2.4 * radius - dist) / 2 * (d / dist if dist > 1e-9 else np.array([1.0, 0.0]))
                    offs[ki] += push
                    offs[kj] -= push
                    moved = True
        if not moved:
            break
    return events, offs


def render(dem: DEMGrid, items: list[ReportItem], opt: ReportOptions | None = None,
           faults=None, quakes=None) -> bytes:
    """Draw the report figure and return its bytes in ``opt.fmt``.

    ``faults``: list of {"lines": [array(N,2) lon/lat], "properties": {...}}.
    ``quakes``: array (N, 3) lon, lat, magnitude.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.patheffects as pe
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    opt = opt or ReportOptions()
    color_by_event = opt.color_by_event
    events = sorted({it.event_label for it in items})
    if color_by_event is None:
        color_by_event = len(events) > 1
    for it in items:
        if it.color is None:
            idx = events.index(it.event_label) if color_by_event else 0
            it.color = SWATH_COLORS[idx % len(SWATH_COLORS)]

    W, H = dem.x_max - dem.x_min, dem.y_max - dem.y_min
    fs = 9
    # ---- figure geometry (inches): ticks/legend/table stack below the map
    below = 0.55 + (0.55 if opt.legend else 0) + (0.3 + 0.22 * (len(items) + 1) if opt.details else 0)
    if PAGES.get(opt.page) is None:
        map_w = opt.width_in
        map_h = map_w * H / W
        top = 0.5 + (0.45 if opt.title else 0)
        fig = plt.figure(figsize=(map_w + 1.6, map_h + below + top))
        fw, fh = fig.get_size_inches()
        ax = fig.add_axes([0.85 / fw, below / fh, map_w / fw, map_h / fh])
    else:
        fw, fh = PAGES[opt.page]
        fig = plt.figure(figsize=(fw, fh))
        top_pad = 0.9 if opt.title else 0.55
        bot_pad = below
        avail_w, avail_h = fw - 1.7, fh - top_pad - bot_pad
        scale = min(avail_w / W, avail_h / H)
        map_w, map_h = W * scale, H * scale
        ax = fig.add_axes([(fw - map_w) / 2 / fw, (bot_pad + (avail_h - map_h)) / fh, map_w / fw, map_h / fh])

    ext = (dem.x_min, dem.x_max, dem.y_min, dem.y_max)
    first = items[0].result if items else None

    # ---- terrain
    if opt.light == "plane" and first is not None:
        light_az = (first.params["grid_dip_azimuth"] + 180.0) % 360.0
    elif opt.light == "nw":
        light_az = 315.0
    else:
        light_az = float(opt.light_azimuth)
    img = terrain_rgb(dem, opt.terrain, light_az, opt.light_altitude, opt.exaggeration)
    ax.imshow(img, extent=ext, origin="upper", interpolation="bilinear", zorder=0)
    ax.set_facecolor("#f4f4f2")

    x, y = dem.xy()
    X, Y = np.meshgrid(x, y)
    halo = [pe.withStroke(linewidth=3.5, foreground="white")]

    # ---- seismo-lineaments
    for it in items:
        r = it.result
        m = r.margin if r.margin is not None else np.where(r.swath, -1.0, 1.0)
        m = np.where(np.isnan(m), 1e9, m)
        if opt.swath_fill:
            ax.contourf(X, Y, m, levels=[-1e15, 0], colors=["white"], alpha=0.30, zorder=5)
        if opt.swath_edges:
            ax.contour(X, Y, m, levels=[0], colors=[it.color], linewidths=2.4, linestyles="--", zorder=8)
        if opt.bounding_traces:
            for ps in r.planes[1:]:
                sd = r.signed_distances.get(ps.name)
                if sd is None or ps.half_width <= 0:
                    continue
                ax.contour(X, Y, np.nan_to_num(sd), levels=[-ps.half_width, ps.half_width],
                           colors=[it.color], linewidths=0.6, alpha=0.8, zorder=7,
                           linestyles="solid")
        if opt.mean_trace:
            ax.contour(X, Y, np.nan_to_num(r.mean_signed_distance), levels=[0], colors=["#111111"],
                       linewidths=2.4, alpha=0.85, zorder=9)

    # ---- faults and seismicity
    if opt.faults and faults:
        for f in faults:
            for line in f["lines"]:
                fx, fy, _ = lonlat_to_utm(line[:, 0], line[:, 1], dem.zone)
                ax.plot(fx, fy, color=FAULT_COLOR, lw=1.5, alpha=0.9, zorder=10, solid_capstyle="round")
    if opt.seismicity and quakes is not None and len(quakes):
        q = np.asarray(quakes, float)
        qx, qy, _ = lonlat_to_utm(q[:, 0], q[:, 1], dem.zone)
        size = 6 * np.clip(q[:, 2], 0.5, 8) ** 2 if q.shape[1] > 2 else 20
        ax.scatter(qx, qy, s=size, facecolors="none", edgecolors=QUAKE_COLOR, linewidths=0.9, zorder=11)

    ax.set_xlim(dem.x_min, dem.x_max)
    ax.set_ylim(dem.y_min, dem.y_max)
    ax.set_aspect("equal")

    # ---- beachballs and labels
    data_per_inch = W / map_w
    radius = opt.beachball_size / 2 * data_per_inch
    if (opt.beachballs or opt.labels) and items:
        evs, offs = _beachball_spots(items, dem, radius, opt.beachball_offset)
        for key, ev in evs.items():
            it0 = ev["items"][0]
            p = it0.result.plane
            ex, ey = ev["x"], ev["y"]
            bx, by = np.array([ex, ey]) + offs[key]
            if np.hypot(*(offs[key])) > 0:
                ax.plot([ex, bx], [ey, by], color="k", lw=1.0, zorder=29,
                        path_effects=[pe.withStroke(linewidth=2.6, foreground="white")])
            ax.plot(ex, ey, marker="*", ms=9, mfc="#ffd400", mec="k", mew=0.8, zorder=31)
            if opt.beachballs and p.rake is not None:
                _beachball(ax, bx, by, radius, p.strike, p.dip, p.rake)
            if opt.labels:
                nums = [it.label[len(key):] or it.label for it in ev["items"]]
                lab = key + ",".join(nums) if all(n.isdigit() for n in nums) else ", ".join(it.label for it in ev["items"])
                ax.text(bx - radius * 0.95, by + radius * 1.05, lab, ha="left", va="bottom",
                        fontsize=fs + 9, fontweight="bold", color="black", zorder=33,
                        path_effects=[pe.withStroke(linewidth=4, foreground="white")])
        # plane labels where each mean trace leaves the map (useful with several swaths)
        if opt.labels and len(items) > 1 and opt.mean_trace:
            for it in items:
                lines = it.result.mean_trace(lonlat=False)
                if not lines:
                    continue
                l = max(lines, key=len)
                seg = np.hypot(*np.diff(l, axis=0).T)
                cum = np.concatenate([[0], np.cumsum(seg)])
                if cum[-1] <= 0:
                    continue
                inset = min(0.07 * min(W, H), 0.3 * cum[-1])
                # label the end farther from the epicentre
                ex0, ey0 = it.result.focus_xyz[:2]
                far_end = np.hypot(l[-1, 0] - ex0, l[-1, 1] - ey0) >= np.hypot(l[0, 0] - ex0, l[0, 1] - ey0)
                tpos = cum[-1] - inset if far_end else inset
                end = (np.interp(tpos, cum, l[:, 0]), np.interp(tpos, cum, l[:, 1]))
                ax.text(end[0], end[1], it.label, fontsize=fs + 2, fontweight="bold", color="k",
                        ha="center", va="center", zorder=34, clip_on=True,
                        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=it.color, lw=1.2))

    # ---- frame furniture
    if opt.graticule:
        _graticule(ax, dem, fs)
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    if opt.scale_bar:
        _scale_bar(ax, dem, fs)
    if opt.north_arrow and first is not None:
        _north_arrow(ax, dem, first.gamma, fs)
    for s in ax.spines.values():
        s.set_linewidth(0.9)

    if opt.title:
        ax.set_title(opt.title, fontsize=fs + 5, fontweight="bold", pad=22 if opt.graticule else 8)

    if opt.legend:
        handles = []
        seen = set()
        for it in items:
            if it.event_label in seen and color_by_event:
                continue
            seen.add(it.event_label)
            handles.append(Line2D([], [], color=it.color, lw=2.4, ls="--",
                                  label=f"Seismo-lineament {it.event_label if color_by_event else ''}".strip()))
            if not color_by_event:
                break
        handles.append(Line2D([], [], color="#111111", lw=2.4, label="Mean nodal-plane trace"))
        if opt.faults and faults:
            handles.append(Line2D([], [], color=FAULT_COLOR, lw=1.8, label="Quaternary fault (USGS)"))
        if opt.seismicity and quakes is not None and len(quakes):
            handles.append(Line2D([], [], marker="o", ls="", mfc="none", mec=QUAKE_COLOR, label="Seismicity"))
        ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.045 if opt.graticule else -0.01),
                  ncol=min(len(handles), 3), fontsize=fs - 0.5, frameon=False, handlelength=2.6)

    if opt.details and items:
        rows = []
        for it in items:
            d = it.details or {}
            p = it.result.plane
            rows.append([it.label, d.get("date", ""), d.get("mag", ""),
                         f"{it.result.hypocenter.lat:.4f}, {it.result.hypocenter.lon:.4f}",
                         f"{it.result.hypocenter.depth_km:.2f}",
                         f"{p.strike:.0f}/{p.dip:.0f}/{'' if p.rake is None else f'{p.rake:.0f}'}",
                         f"±{p.strike_uncertainty:.0f}/±{p.dip_uncertainty:.0f}",
                         d.get("location_unc", ""), d.get("source", "")])
        cols = ["Label", "Date (UTC)", "M", "Lat, Lon", "Depth km", "Strike/dip/rake",
                "σ strike/dip", "Location ±", "Source"]
        bbox = ax.get_position()
        h_in = 0.22 * (len(rows) + 1)
        top_of_table = bbox.y0 * fh - 0.55 - (0.55 if opt.legend else 0) - 0.1
        tax = fig.add_axes([bbox.x0, (top_of_table - h_in) / fh, bbox.width, h_in / fh])
        tax.axis("off")
        t = tax.table(cellText=rows, colLabels=cols, loc="center", cellLoc="left", colLoc="left",
                      edges="horizontal")
        t.auto_set_font_size(False)
        t.set_fontsize(fs - 1)
        t.auto_set_column_width(list(range(len(cols))))
        for (r_, c_), cell in t.get_celld().items():
            cell.set_linewidth(0.5)
            if r_ == 0:
                cell.set_text_props(fontweight="bold")

    buf = io.BytesIO()
    fig.savefig(buf, format=opt.fmt, dpi=opt.dpi if opt.fmt == "png" else None,
                bbox_inches="tight", pad_inches=0.15, facecolor="white")
    plt.close(fig)
    return buf.getvalue()


def terrain_image(dem: DEMGrid, style="hillshade", azimuth=315.0, altitude=15.0, exaggeration=2.0) -> bytes:
    """Terrain alone as a PNG at the DEM's native resolution (one pixel per cell)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    img = terrain_rgb(dem, style, azimuth, altitude, exaggeration)
    buf = io.BytesIO()
    plt.imsave(buf, np.clip(img, 0, 1), format="png")
    return buf.getvalue()
