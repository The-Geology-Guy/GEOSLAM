"""Hypocentre uncertainty and the half-width of a SLAM swath.

Cronin (SLAMcode20210103.nb, ``findHalfWidth``) bounds each candidate fault
plane by two parallel planes tangent to the hypocentre uncertainty ellipsoid.
His module rotates the plane normal into the ellipsoid's principal frame
(axes along ``eh1Az``, perpendicular to it horizontally, and vertical), finds
the tangent point and returns its distance from the plane through the centre.
Algebraically that is the support function of the ellipsoid::

    h(n) = sqrt( (n . e1)^2 a1^2 + (n . e2)^2 a2^2 + (n . e3)^2 a3^2 ) = sqrt(n^T C n)

with semi-axes a_i along unit axes e_i and C = R diag(a^2) R^T.  This module
implements that general form, so a full 3-D QuakeML confidence ellipsoid
(azimuth, plunge, rotation) can be used, while Cronin's horizontal-ellipse
case (eh1, eh1Az, eh2, ez) is reproduced exactly.

Also provided, for history and comparison, the earlier box model of Bayliss
(2007, Eq. 7): ``w = ez cos(dip) + eh sin(dip)``.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .orientation import make_vector


@dataclass
class HypocenterUncertainty:
    """Confidence ellipsoid of a hypocentre, semi-axes in metres.

    ``axes`` is a 3x3 array whose rows are the unit principal axes (x=E, y=N,
    z=up) and ``semi_axes`` their lengths.  Use the constructors rather than
    building this directly.  ``confidence`` is informational (e.g. 0.90, 0.95).
    """

    semi_axes: np.ndarray
    axes: np.ndarray
    confidence: float | None = None
    source: str = ""
    covariance: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        self.semi_axes = np.asarray(self.semi_axes, float)
        self.axes = np.asarray(self.axes, float)
        if np.any(self.semi_axes < 0):
            raise ValueError("semi-axes must be non-negative")
        self.covariance = self.axes.T @ np.diag(self.semi_axes**2) @ self.axes

    # --- constructors -------------------------------------------------
    @classmethod
    def from_cronin(cls, eh1, eh1_azimuth, eh2, ez, confidence=None, source="cronin"):
        """Cronin's inputs: horizontal semi-axes eh1 (along ``eh1_azimuth``,
        grid azimuth in degrees) and eh2 (perpendicular), vertical ez; metres."""
        a = np.radians(eh1_azimuth)
        e1 = np.array([np.sin(a), np.cos(a), 0.0])
        e2 = np.array([-np.cos(a), np.sin(a), 0.0])
        e3 = np.array([0.0, 0.0, 1.0])
        return cls(np.array([eh1, eh2, ez]), np.vstack([e1, e2, e3]), confidence, source)

    @classmethod
    def isotropic_horizontal(cls, eh, ez, confidence=None, source="eh/ez"):
        return cls.from_cronin(eh, 0.0, eh, ez, confidence, source)

    @classmethod
    def from_quakeml_ellipsoid(cls, semi_major, semi_intermediate, semi_minor,
                               major_azimuth, major_plunge, major_rotation=0.0,
                               confidence=None, source="quakeml"):
        """QuakeML ``ConfidenceEllipsoid`` (metres, degrees).

        The major axis has the given azimuth/plunge.  The intermediate axis is
        taken perpendicular to it and horizontal-most before applying
        ``major_rotation`` about the major axis.  When rotation is 0 (the usual
        ComCat value) the convention ambiguity in the QuakeML spec disappears.
        Azimuths here must already be relative to grid north.
        """
        e1 = make_vector(major_plunge, major_azimuth)
        # horizontal vector perpendicular to e1
        h = np.cross(e1, [0.0, 0.0, 1.0])
        if np.linalg.norm(h) < 1e-12:  # vertical major axis
            h = np.array([1.0, 0.0, 0.0])
        h = h / np.linalg.norm(h)
        k = np.cross(e1, h)
        r = np.radians(major_rotation)
        e2 = np.cos(r) * h + np.sin(r) * k
        e3 = np.cross(e1, e2)
        return cls(np.array([semi_major, semi_intermediate, semi_minor]),
                   np.vstack([e1, e2, e3]), confidence, source)

    def rotated_azimuths(self, delta_deg: float) -> "HypocenterUncertainty":
        """Rotate the ellipsoid about the vertical by -delta (true->grid north)."""
        t = np.radians(-delta_deg)
        # azimuth increases clockwise; rotating azimuths by -delta is a
        # counter-clockwise rotation of vectors by delta in the (E, N) plane
        rz = np.array([[np.cos(t), np.sin(t), 0.0], [-np.sin(t), np.cos(t), 0.0], [0, 0, 1.0]])
        return HypocenterUncertainty(self.semi_axes.copy(), (rz @ self.axes.T).T,
                                     self.confidence, self.source)

    def scaled(self, factor: float) -> "HypocenterUncertainty":
        return HypocenterUncertainty(self.semi_axes * factor, self.axes.copy(),
                                     self.confidence, self.source)

    # --- evaluation ---------------------------------------------------
    def half_width(self, normal) -> np.ndarray:
        """Distance from the ellipsoid centre to the tangent planes with the
        given unit normal(s) (support function).  Vectorised on the last axis."""
        n = np.asarray(normal, float)
        return np.sqrt(np.einsum("...i,ij,...j->...", n, self.covariance, n))


def find_half_width_cronin(dip_az, dip_ang, h_er1, h_er1_az, h_er2, v_er):
    """Literal transcription of Cronin's ``findHalfWidth`` (for verification)."""
    var1 = h_er1_az - 90.0
    d, a = np.radians(dip_ang), np.radians(dip_az)
    var2 = np.array([np.cos(d) * np.sin(a), np.cos(d) * np.cos(a), -np.sin(d)])
    b = np.radians(dip_az - 90.0)
    var3 = np.array([np.sin(b), np.cos(b), 0.0])
    c = np.radians(var1)
    var4 = np.array([[np.cos(c), -np.sin(c), 0.0], [np.sin(c), np.cos(c), 0.0], [0.0, 0.0, 1.0]])
    var6 = var4 @ np.cross(var2, var3)
    var7 = (1 / h_er1**2) * (var6[0] * (h_er1**2 / 2)) ** 2
    var8 = (1 / h_er2**2) * (var6[1] * (h_er2**2 / 2)) ** 2
    var9 = (1 / v_er**2) * (var6[2] * (v_er**2 / 2)) ** 2
    var10 = np.sqrt(1 / (var7 + var8 + var9))
    var11 = np.array([var10 * var6[0] * h_er1**2 / 2, var10 * var6[1] * h_er2**2 / 2,
                      var10 * var6[2] * v_er**2 / 2])
    return float(var6 @ var11)


def box_half_width_bayliss(dip, eh, ez):
    """Bayliss (2007) Eq. 7 box-model half-width: ez cos(dip) + eh sin(dip)."""
    d = np.radians(dip)
    return ez * np.cos(d) + eh * np.sin(d)
