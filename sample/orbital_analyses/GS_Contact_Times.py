# Utilized Modules
import numpy as np


def satContact(sat_R, gs_R):
    """
    Determines if satellite is within sight of a Ground Station

    Parameters
    ----------
    sat_R : numpy matrix [3, 1]
        - Input radius vector in Inertial System ([[X], [Y], [Y]])
    gs_R : numpy matrix [3, 1]
        - Input radius vector in Inertial System ([[X], [Y], [Y]])

    Returns
    -------
    inContact : int
        - 1 or 0 if sat is in sight our out of sight, respectively

    See Also
    --------
    Sun_Contact_Times : Determine if a orbit vector is illuminated

    Geo_Contact_Times : Determine if a orbit vector is in within a
    geometric boundary

    References
    ----------
    [1] D. Vallado, `Fundamentals of Astrodynamics and Applications`. 4th ed.,
    Microcosm Press, 2013.
        - Modified Alg. 35, pg. 308
    """
    # Simplifying equations
    mag_sat = np.linalg.norm(sat_R)
    mag_gs = np.linalg.norm(gs_R)
    dot_ss = np.dot(np.transpose(sat_R), gs_R)

    # Find minimum parametric value
    Tmin = (((mag_sat ** 2) - dot_ss) /
            ((mag_sat ** 2) + (mag_gs ** 2) - 2 * (dot_ss)))
    if Tmin < 0 or Tmin > 1:
        InContact = 1  # Satellite can see GS
    if Tmin > 0 and Tmin < 1:
        cTmin = (((1 - Tmin) * (mag_sat ** 2) +
                  (dot_ss * Tmin)) / (6378.137 ** 2))
        if cTmin > 1:
            InContact = 1  # Satellite can see GS
        if cTmin < 1:
            InContact = 0  # Satellite can't see GS
    return InContact
