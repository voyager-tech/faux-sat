# Utilized Modules
import numpy as np

from sample.orbital_analyses.Transform_State import TimeAdjust
from sample.orbital_analyses.Transform_Coordinate import FK5_Precession
# Sun Position Vector Determination


def sunPosition(GD_UTC):
    """
    Determines current lcation of the Sun
        - Used for determiing if satellite is illuminated

    Parameters
    ----------
    GD_UTC : numpy matrix [6, 1]
        - Input Gregorian Date

    Returns
    -------
    rad_km : numpy matrix [3, 1]
        - Return the inertial frame vector of the Sun

    See Also
    --------
    Sun_Contact_Times : Determine if a satellite is illuminated

    References
    ----------
    [1] D. Vallado, `Fundamentals of Astrodynamics and Applications`. 4th ed.,
    Microcosm Press, 2013.
        - Alg. 29, pg. 279-280
    """
    # Inaccurate to a small degree
    UTC, JD_UT1, TAI, TT = TimeAdjust(GD_UTC)  # UTC Input
    UT1 = ((JD_UT1 - 2451545.0) / 36525)
    mlon = (280.460 + (36000.771 * UT1))  # degrees
    # Reduce sp_mlon to within 360 degrees
    mlon_d = np.mod(mlon, 360, dtype=np.float64)  # Degrees

    TDB = UT1
    M = (357.5291092 + (35999.05034 * TDB))  # degrees

    # Reduce sp_M to within 360 degrees then convert to radians
    reduc2 = (M / 360)
    diff2 = reduc2 - np.floor(reduc2)
    M = (360 * diff2)  # degrees
    M_r = np.deg2rad(M)  # rad

    ecl = (mlon_d + 1.914666471 * np.sin(M_r) +
           0.019994643 * np.sin(2 * (M_r)))  # deg
    ecl_r = np.deg2rad(ecl)  # rad
    dist = (1.000140612 - 0.016708617 * np.cos(M_r) -
            0.000139589 * np.cos(2 * M_r))  # AU
    eps = (23.439291 - 0.0130042 * TDB)
    eps_r = (eps * (np.pi / 180))  # rad
    rad = np.zeros((3, 1))
    rad[0, 0] = (dist * np.cos(ecl_r))
    rad[1, 0] = (dist * np.cos(eps_r) * np.sin(ecl_r))
    rad[2, 0] = (dist * np.sin(eps_r) * np.sin(ecl_r))
    rad_kilo = (rad * 149597870.7)  # km

    # In Geocentric equitorial coordinates (MOD)
    # Conversion from MOD to GCRF to J2000
    # ASSUMES SUN IS IN MODEq, need additional transform for MODEc
    zero = np.zeros((3, 1), dtype=np.float64)
    rad_gcrf, vel_gcrf = FK5_Precession(rad_kilo, zero, GD_UTC)

    # Assuming GCRF and J2000 are closely alligned
    rad_km = rad_gcrf
    return rad_km


def satIlluminated(sat_R, gd_UTC):
    """
    Determines if satellite is illuminated by the Sun

    Parameters
    ----------
    sat_R : numpy matrix [3, 1]
        - Input radius vector in Inertial System ([[X], [Y], [Y]])
    gd_UTC : numpy matrix [6, 1]
        - Input Gregorian Date

    Returns
    -------
    inSun : int
        - 1 or 0 if sat is illuminated or in the dark, respectively

    See Also
    --------
    GS_Contact_Times : Determine if a orbit vector is in sight of a
    Ground Station

    Geo_Contact_Times : Determine if a orbit vector is in within a
    geometric boundary

    References
    ----------
    [1] D. Vallado, `Fundamentals of Astrodynamics and Applications`. 4th ed.,
    Microcosm Press, 2013.
        - Alg. 35 + Light Routine, pg. 309
    """
    # Determine Sight and Light Conditions
    # Determine sun position from sunPosition
    sun_R = sunPosition(gd_UTC)

    # Simplifying equations
    mag_sat = np.linalg.norm(sat_R)
    mag_sun = np.linalg.norm(sun_R)
    dot_ss = np.dot(np.transpose(sat_R), sun_R)

    # Find minimum parametric value
    Tmin = (((mag_sat ** 2) - dot_ss) /
            ((mag_sat ** 2) + (mag_sun ** 2) - 2 * (dot_ss)))
    if Tmin < 0 or Tmin > 1:
        InSun = 1  # Satellite is illuminated
    if Tmin > 0 and Tmin < 1:
        cTmin = (((1 - Tmin) * (mag_sat ** 2) +
                 (dot_ss * Tmin)) / (6378.137 ** 2))
        if cTmin > 1:
            InSun = 1  # Satellite is illuminated
        if cTmin < 1:
            InSun = 0  # Satellite is in the Dark
    return InSun
