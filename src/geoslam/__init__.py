"""GEOSLAM - Python implementation of Vincent Cronin's Seismo-Lineament Analysis Method.

SLAM projects the nodal planes of a shallow earthquake's focal mechanism,
together with the uncertainty in hypocentre location and plane orientation, to
the ground surface (a DEM) to outline the swath - the *seismo-lineament* -
that should contain the surface trace of the seismogenic fault.

Method and original Mathematica code: Vince Cronin (Baylor University),
https://croninprojects.org/Vince/SLAM/ ; Cronin et al. (2008) Environmental &
Engineering Geoscience 14(3):199-219.  Python implementation: Luke Pajer.
"""
from .version import __version__
from .geodesy import UTMZone, grid_convergence, lonlat_to_utm, utm_to_lonlat, utm_zone
from .orientation import NodalPlane, auxiliary_plane, plane_normal
from .uncertainty import HypocenterUncertainty
from .dem import DEMGrid, GeographicDEM, read_geotiff
from .slam import Hypocenter, SLAMResult, seismo_lineament, suggest_half_extent

__all__ = [
    "__version__", "UTMZone", "grid_convergence", "lonlat_to_utm", "utm_to_lonlat", "utm_zone",
    "NodalPlane", "auxiliary_plane", "plane_normal", "HypocenterUncertainty", "DEMGrid",
    "GeographicDEM", "read_geotiff", "Hypocenter", "SLAMResult", "seismo_lineament",
    "suggest_half_extent",
]
