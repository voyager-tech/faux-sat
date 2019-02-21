# Utilized Modules
from matplotlib.path import Path
import sample.Requirements as Req
# http://econym.org.uk/gmap/states.xml


def satOverhead(sat_geo):
    """
    Determines if satellite is within a geometric boundary defined by
    an array of lattitude and longitude points.
        - Boundary Data from .npy file imported in Requirements.py

    Parameters
    ----------
    sat_geo : numpy matrix [3, 1]
        - Input data in Geodetic Coordinates ([[Lat], [Long], [Alt]])

    Returns
    -------
    in_out : int
        - 1 or 0 if sat is inside or outside the boundary, respectively

    See Also
    --------
    Sun_Contact_Times : Determine if a orbit vector is illuminated

    GS_Contact_Times : Determine if a orbit vector is in sight of a
    Ground Station
    """
    # Define path that encircles california through .npy coordinates
    CA_geo = Req.Target
    CA_bound = Path(CA_geo)

    # Check if inside or outside boundary
    in_out = CA_bound.contains_point((sat_geo[0, 0], sat_geo[1, 0]),
                                     transform=None, radius=0.0)
    return in_out
