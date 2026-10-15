"""Tabular input: Cronin's SLAM spreadsheet layout and a named-column CSV.

Cronin's ``Sample-SLAM-FM-data.xlsx`` is read by column position
(SLAMcode20210103.nb, 1-based)::

    7 latitude   8 longitude   9 depth (km)
    10 eh1 (km)  11 eh2 (km)   12 eh1 azimuth (deg)   13 ez (km)
    15 dip trend (deg, true N)  16 dip (deg)
    18 dip-trend uncertainty    19 dip uncertainty
    21 UTM easting  22 UTM northing  (optional; recomputed here from lat/lon)

Columns 17 and 20 are taken as rake and rake uncertainty when numeric.

The named-column CSV (``read_events_csv``) needs::

    event_id, lat, lon, depth_km,
    eh1_km, eh2_km, eh1_az, ez_km,                 (optional)
    dip_azimuth | strike, dip, rake,
    strike_unc | dip_azimuth_unc, dip_unc, rake_unc  (optional)
"""
from __future__ import annotations

import numpy as np

from .orientation import NodalPlane
from .slam import Hypocenter
from .uncertainty import HypocenterUncertainty


def _row_to_case(lat, lon, depth, eh1, eh2, eh1az, ez, dipaz, dip, rake, sa, sd, sr, eid=""):
    unc = None
    if all(np.isfinite([eh1, eh2, ez])) and (eh1 > 0 or ez > 0):
        unc = HypocenterUncertainty.from_cronin(eh1 * 1000, eh1az if np.isfinite(eh1az) else 0.0,
                                                eh2 * 1000, ez * 1000, source="table (km)")
    hyp = Hypocenter(float(lon), float(lat), float(depth), unc, event_id=str(eid))
    plane = NodalPlane.from_dip_azimuth(float(dipaz), float(dip),
                                        float(rake) if np.isfinite(rake) else None,
                                        float(sa) if np.isfinite(sa) else 0.0,
                                        float(sd) if np.isfinite(sd) else 0.0,
                                        float(sr) if np.isfinite(sr) else 0.0)
    return hyp, plane


def read_cronin_table(path, sheet=0):
    """Rows of Cronin's input sheet -> list of (Hypocenter, NodalPlane)."""
    import pandas as pd

    df = pd.read_excel(path, sheet_name=sheet, header=None) if str(path).endswith(("xls", "xlsx")) \
        else pd.read_csv(path, header=None)
    num = df.apply(pd.to_numeric, errors="coerce")
    out = []
    for i, r in num.iterrows():
        v = r.values
        def c(k):
            return v[k - 1] if len(v) >= k else np.nan
        if not np.isfinite(c(7)) or not np.isfinite(c(16)):
            continue  # header or blank row
        out.append(_row_to_case(c(7), c(8), c(9), c(10), c(11), c(12), c(13), c(15), c(16),
                                c(17), c(18), c(19), c(20), eid=str(df.iloc[i, 0])))
    return out


def read_events_csv(path):
    import pandas as pd

    df = pd.read_csv(path)
    g = lambda r, *ks: next((float(r[k]) for k in ks if k in r and pd.notna(r[k])), np.nan)
    out = []
    for _, r in df.iterrows():
        dipaz = g(r, "dip_azimuth")
        if not np.isfinite(dipaz):
            dipaz = (g(r, "strike") + 90.0) % 360.0
        out.append(_row_to_case(g(r, "lat"), g(r, "lon"), g(r, "depth_km"), g(r, "eh1_km"),
                                g(r, "eh2_km"), g(r, "eh1_az"), g(r, "ez_km"), dipaz, g(r, "dip"),
                                g(r, "rake"), g(r, "strike_unc", "dip_azimuth_unc"), g(r, "dip_unc"),
                                g(r, "rake_unc"), eid=r.get("event_id", "")))
    return out
