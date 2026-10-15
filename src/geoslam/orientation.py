"""Nodal-plane orientation vectors and focal-mechanism conventions.

Frame: x = (grid) east, y = (grid) north, z = up.  Angles in degrees.
Azimuths are measured clockwise from north; plunges are positive downward.

Cronin's definitions (SLAMcode20210103.nb)::

    makeVector[plunge_, trend_] := {Cos[plunge] Sin[trend], Cos[plunge] Cos[trend], -Sin[plunge]}
    npStrike     = If[npDipTr < 90, npDipTr + 270, npDipTr - 90]      (right-hand rule)
    npStrikeVect = {Sin[npStrike], Cos[npStrike], 0}
    npNormalVect = unitVector[Cross[npDipVect, npStrikeVect]]          (upward normal)

The cross product of the down-dip unit vector with the RHR strike unit vector
is the upward unit normal::

    n = (sin(dip) sin(dipaz), sin(dip) cos(dipaz), cos(dip))

That closed form is used here.  It is continuous for any real ``dip``: a "dip"
greater than 90 deg, or less than 0, describes the same plane as Cronin's
fold-over rule (dip' = 180 - dip, or |dip|, with the dip azimuth flipped by
180 deg), with the normal reversed.  Distances to the plane are taken in
absolute value, so the sign of the normal never matters.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def wrap360(a):
    return np.asarray(a, float) % 360.0


def make_vector(plunge, trend):
    """Unit vector from plunge/trend (Cronin ``makeVector``)."""
    p, t = np.radians(plunge), np.radians(trend)
    return np.array([np.cos(p) * np.sin(t), np.cos(p) * np.cos(t), -np.sin(p)])


def strike_from_dip_azimuth(dip_azimuth):
    """Right-hand-rule strike: dip direction minus 90 deg."""
    return wrap360(np.asarray(dip_azimuth, float) - 90.0)


def dip_azimuth_from_strike(strike):
    return wrap360(np.asarray(strike, float) + 90.0)


def plane_normal(dip, dip_azimuth):
    """Upward unit normal(s) of plane(s) with the given dip and dip azimuth.

    Broadcasts over array inputs; the vector components are on the last axis.
    """
    d, a = np.radians(dip), np.radians(dip_azimuth)
    d, a = np.broadcast_arrays(d, a)
    return np.stack([np.sin(d) * np.sin(a), np.sin(d) * np.cos(a), np.cos(d)], axis=-1)


def plane_normal_cronin(dip, dip_azimuth):
    """Normal built exactly as in Cronin's notebook (for cross-checking)."""
    dip_vec = make_vector(dip, dip_azimuth)
    strike = dip_azimuth + 270.0 if dip_azimuth < 90.0 else dip_azimuth - 90.0
    strike_vec = np.array([np.sin(np.radians(strike)), np.cos(np.radians(strike)), 0.0])
    n = np.cross(dip_vec, strike_vec)
    return n / np.linalg.norm(n)


def fold_plane(dip, dip_azimuth):
    """Normalise (dip, dip azimuth) to 0 <= dip <= 90 using Cronin's fold-over rule."""
    dip = float(dip)
    az = float(dip_azimuth) % 360.0
    if dip < 0.0:
        dip, az = -dip, (az + 180.0) % 360.0
    if dip > 90.0:
        dip, az = 180.0 - dip, (az + 180.0) % 360.0
    return dip, az


def slip_vector(strike, dip, rake):
    """Hanging-wall slip unit vector (x=E, y=N, z=up) from Aki & Richards angles."""
    s, d, r = np.radians([strike, dip, rake])
    # Aki & Richards (x=N, y=E, z=down) -> convert to (E, N, Up)
    north = np.cos(r) * np.cos(s) + np.sin(r) * np.cos(d) * np.sin(s)
    east = np.cos(r) * np.sin(s) - np.sin(r) * np.cos(d) * np.cos(s)
    down = -np.sin(r) * np.sin(d)
    return np.array([east, north, -down])


def auxiliary_plane(strike, dip, rake):
    """Second nodal plane (strike, dip, rake) of a double couple.

    Uses the normal/slip swap: the auxiliary plane's normal is the slip vector
    of the fault plane and vice versa.
    """
    n = plane_normal(dip, dip_azimuth_from_strike(strike))  # upward normal
    u = slip_vector(strike, dip, rake)
    # The auxiliary plane has normal u and slip n.  Orient its normal upward.
    n2, u2 = (u, n) if u[2] >= 0 else (-u, -n)
    dip2 = np.degrees(np.arccos(np.clip(n2[2], -1.0, 1.0)))
    dipaz2 = np.degrees(np.arctan2(n2[0], n2[1])) % 360.0
    strike2 = (dipaz2 - 90.0) % 360.0
    if dip2 < 1e-9:  # horizontal plane: strike undefined, choose slip direction
        strike2 = 0.0
    # rake: angle in the plane from strike to slip, anticlockwise looking at hanging wall
    s2 = np.radians(strike2)
    strike_vec = np.array([np.sin(s2), np.cos(s2), 0.0])
    dip_vec = make_vector(dip2, dipaz2)
    rake2 = np.degrees(np.arctan2(-np.dot(u2, dip_vec), np.dot(u2, strike_vec)))
    return float(strike2), float(dip2), float(rake2)


@dataclass
class NodalPlane:
    """One nodal plane, with 1-sigma-style uncertainties as published.

    ``dip_azimuth`` is measured from *true* north.  Uncertainties follow
    Cronin's input sheet: the dip-trend uncertainty equals the strike
    uncertainty (they are the same angle), and the dip uncertainty.
    """

    strike: float
    dip: float
    rake: float | None = None
    strike_uncertainty: float = 0.0
    dip_uncertainty: float = 0.0
    rake_uncertainty: float = 0.0
    label: str = ""

    @classmethod
    def from_dip_azimuth(cls, dip_azimuth, dip, rake=None, dip_azimuth_uncertainty=0.0,
                         dip_uncertainty=0.0, rake_uncertainty=0.0, label=""):
        return cls(float(strike_from_dip_azimuth(dip_azimuth)), dip, rake,
                   dip_azimuth_uncertainty, dip_uncertainty, rake_uncertainty, label)

    @property
    def dip_azimuth(self) -> float:
        return float(dip_azimuth_from_strike(self.strike))

    def auxiliary(self) -> "NodalPlane":
        if self.rake is None:
            raise ValueError("rake is required to compute the auxiliary plane")
        s, d, r = auxiliary_plane(self.strike, self.dip, self.rake)
        return NodalPlane(s, d, r, label=(self.label + " aux").strip())
