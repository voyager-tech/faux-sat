# Utilized Modules
import numpy as np
from copy import deepcopy as dc

from orbital_analyses.Transform_State import TimeAdjust

# ********************** Coordinate Reference Sheet **********************
# *** Earth Based Systems ***
# Geocentric Equatorial Coordinate System (IJK) (ECI) (GEI) (GCI) - Inertial
# Geocentric Celestial Coordinate System (GCRF) - Inertial
# Body-Fixed Coordinate System (ITRF) (ECEF) - Body Fixed
# Topocentric Horizon Coordinate System (SEZ)
# Topocentric Equitorial Coordinate System (ItJtKt)

# *** Satellite Based Systems ***
# Perifocal Coordinate System (PQW)
# Satellite Radial System (RSW)
# Satellite Normal System (NTW)
# Equinoctial Coordinate System (EQW)

# *** Interplanetary Systems ***
# Heliocentric Coordinate System (XYZ)
# International Celestial Reference System (ICRF)
# ***********************************************************************


def ECFixed2Geodetic(rad):
    """
    Converts ECEF (Fixed, Earth Centered) coordinates to the Geodetic
    coordinate system (lat, long, alt).

    Parameters
    ----------
    rad : numpy matrix [3, 1] - [[X], [Y], [Z]]
        - Radius vector components defined in kilometers

    Returns
    -------
    lla : numpy matrix [3, 1] - [[Latitude], [Longitude], [Altitude]]
        - Latitude defined in decimal degrees
        - Longitude defined in decimal degrees
        - Altitude defined in kilometers from the GRS 80 reference ellipsoid

    See Also
    --------
    Geodetic2ECFixed : The reverse of this transform

    References
    ----------
    [1] EPSG,
    "Coordinate Conversions and Transformations including Formulas"
    Guidance Note 7-2. Section 2.2.1.
        - http://www.ihsenergy.com/epsg/guid7_2.pdf
    """
    # Using GRS 80 Geodetic Reference System for Earth's Ellipsoid Shape
    E_a = 6378137.0000  # m - semimajor axis
    E_b = 6356752.3141  # m - semiminor axis
    E_e2 = 0.00669438002290  # first eccentricity
    E_e2s = 0.00673949677548  # second eccentricity

    # Convert km to m
    rad_m = (rad * 1e3)  # m

    # Calculate intermediate values
    E_p = np.sqrt((rad_m[0, 0] ** 2) + (rad_m[1, 0] ** 2))
    E_q = np.arctan((rad_m[2, 0] * E_a) / (E_p * E_b))

    lat = np.arctan((rad_m[2, 0] + (E_e2s * E_b * (np.sin(E_q) ** 3))) /
                    (E_p - (E_e2 * E_a * (np.cos(E_q) ** 3))))  # rad
    E_v = (E_a / np.sqrt(1 - (E_e2 * (np.sin(lat) ** 2))))
    lon = np.arctan2(rad_m[1, 0], rad_m[0, 0])  # rad
    alt = ((E_p / np.cos(lat)) - E_v)  # m

    # Convert to Degrees
    lat_d = np.rad2deg(lat)  # Deg
    lon_d = np.rad2deg(lon)  # Deg
    if lon_d < 0:
        lon_d = (lon_d + 360)
    if lon_d >= 360:
        lon_d = (lon_d - 360)
    # Convert back to km
    alt_km = (alt * 1e-3)  # km

    # Set output values
    lla = np.asmatrix(np.zeros((3, 1), dtype=np.float64))
    lla[0, 0] = lat_d
    lla[1, 0] = lon_d
    lla[2, 0] = alt_km
    return lla


def Geodetic2ECFixed(lla):
    """
    Converts Geodetic coordinates to the ECEF (Fixed, Earth Centered)
    coordinate system.

    Parameters
    ----------
    lla : numpy matrix [3, 1] - [[Latitude], [Longitude], [Altitude]]
        - Latitude defined in decimal degrees
        - Longitude defined in decimal degrees
        - Altitude defined in kilometers from the GRS 80 reference ellipsoid

    Returns
    -------
    rad_ecef : numpy matrix [3, 1] - [[X], [Y], [Z]]
        - Radius vector components defined in kilometers

    See Also
    --------
    ECFixed2Geodetic : The reverse of this transform

    References
    ----------
    [1] EPSG,
    "Coordinate Conversions and Transformations including Formulas"
    Guidance Note 7-2. Section 2.2.1.
        - http://www.ihsenergy.com/epsg/guid7_2.pdf
    [2] H. Moritz, "Geodetic Reference System 1980"
        - ftp://athena.fsv.cvut.cz/ZFG/grs80-Moritz.pdf
    """
    # Using GRS 80 Geodetic Reference System for Earth's Ellipsoid Shape
    E_a = 6378.137  # km - semimajor axis
    E_e2 = 0.00669438002290  # first eccentricity

    # Ensure input variable is in numpy matrix dtype
    lla = np.asmatrix(lla)

    # Convert degrees to radians
    lat_r = np.deg2rad(lla[0, 0])
    lon_r = np.deg2rad(lla[1, 0])

    # Calculate intermediate values
    C = (E_a / (np.sqrt(1 - (E_e2 * (np.sin(lat_r) ** 2)))))
    S = ((E_a * (1 - E_e2)) / (np.sqrt(1 - E_e2 * (np.sin(lat_r) ** 2))))

    # Calculate radius vector components
    rad_arr = np.zeros((3, 1), dtype=np.float)
    rad_arr[0] = ((C + lla[2, 0]) * np.cos(lat_r) * np.cos(lon_r))  # km
    rad_arr[1] = ((C + lla[2, 0]) * np.cos(lat_r) * np.sin(lon_r))  # km
    rad_arr[2] = ((S + lla[2, 0]) * np.sin(lat_r))  # km
    rad_ecef = np.asmatrix(rad_arr)
    return rad_ecef


def Geodetic2Geocentric(lla_D):
    """
    Converts Geodetic coordinates to Geocentric coordinates.
    Only the lattitude value is changed in the transform.

    Parameters
    ----------
    lla_D : numpy matrix [3, 1] - [[Latitude], [Longitude], [Altitude]]
        - Latitude defined in decimal degrees normal to the ellipsoid
        - Longitude defined in decimal degrees
        - Altitude defined in kilometers from the GRS 80 reference ellipsoid

    Returns
    -------
    lla_C : numpy matrix [3, 1] - [[Latitude], [Longitude], [Altitude]]
        - Latitude defined in decimal degrees
        - Longitude defined in decimal degrees
        - Altitude defined in kilometers from the GRS 80 reference ellipsoid

    See Also
    --------
    Geodetic2ECFixed : Transform to the ECEF (Fixed, Earth Centered)
    coordinate system.

    References
    ----------
    [1] J. Clynch, "Geodetic Coordinate Conversions"
        - http://clynchg3c.com/Technote/geodesy/coordcvt.pdf
    [2] H. Moritz, "Geodetic Reference System 1980"
        - ftp://athena.fsv.cvut.cz/ZFG/grs80-Moritz.pdf
    """
    # Only works well close to earth (Ground -> LEO)

    # Using GRS 80 Geodetic Reference System for Earth's Ellipsoid Shape
    E_a = 6378.137  # km - semimajor axis
    E_e2 = 0.00669438002290  # first eccentricity

    # Convert to radians
    lat_r = np.deg2rad(lla_D[0, 0])

    # Calculate radius of curvature in the prime vertical
    Rn = (E_a / (np.sqrt(1 - (E_e2 * (np.sin(lat_r) ** 2)))))
    lat_C = np.arctan((1 - (E_e2 * (Rn / (Rn + lla_D[2, 0])))) * np.tan(lat_r))
    lla_C = np.matrix([[lat_C * (180 / np.pi)], [lla_D[1, 0]], [lla_D[2, 0]]])
    return lla_C

###############################################################################
###############################################################################


def IAU_PolarMotion(rad_itrf, vel_itrf, gd_UTC, Transpose):
    """
    Transforms vectors from ITRF to TIRS frame following IAU-2010 conventions

    Parameters
    ----------
    rad_itrf : numpy matrix [3, 1] - [[X], [Y], [Z]]
        - Radius vector components defined in kilometers in the ITRF frame
    vel_itrf : numpy matrix [3, 1] - [[VX], [VY], [VZ]]
        - Velocity vector components defined in kilometers in the ITRF frame
    gd_UTC : numpy matrix [6, 1] - [[Yr], [Mo], [Day], [Hr], [Min], [Sec]]
        - Gregorian Date
    Transpose : int of 0 or 1
        - Determines wether the transform is itrf->tirs (0) or tirs->itrf (1)

    Returns
    -------
    rad_tirs : numpy matrix [3, 1] - [[X], [Y], [Z]]
        - Radius vector components defined in kilometers in the TIRS frame
    vel_tirs : numpy matrix [3, 1] - [[VX], [VY], [VZ]]
        - Velocity vector components defined in kilometers in the TIRS frame

    See Also
    --------
    FK5_ECFixed2J2000 : Transform from an earth fixed frame (ECEF) to an
    earth based inertial frame (J2000) using the IAU-76/FK5 reduction

    FK5_J20002ECFixed : Transform from an earth based inertial frame (J2000)
    to an earth fixed frame (ECEF) using the IAU-76/FK5 reduction

    References
    ----------
    [1] D. Vallado, `Fundamentals of Astrodynamics and Applications`. 4th ed.,
    Microcosm Press, 2013.
        - Pg. 212
    """
    # Polar Motion (IAU-2006/2000, CIO Based) ITFR -> TIRS, Vallado Pg. 212

    # Raise errors before running script
    while Transpose not in [0, 1]:
        raise RuntimeError("Enter an int of 0; or 1 for the reverse transform")

    # Time Adjustments
    jd_UTC, jd_UT1, jd_TAI, jd_TT = TimeAdjust(gd_UTC)

    # Transform Constants
    # From Ast Almonac, 2006?:B76
    a_a = (0.12)  # Sexagesimal-Sec  # TODO: Replace with date specific data
    a_c = (0.26)  # Sexagesimal-Sec  # TODO: Replace with date specific data

    # These from IERS Earth Orientation Data - EOP 14 C04 (IAU2000A) - Latest
    # https://datacenter.iers.org/eop/-/somos/5Rgv/latest/223
    # find way to get site data and store it for entirety of the scipt run then delete afterwords
    Xp = np.deg2rad(-0.140682 * (1 / 3600))  # Sexagesimal-sec -> Radians # TODO: Replace with date specific data
    Yp = np.deg2rad(0.333309 * (1 / 3600))  # Sexagesimal-sec -> Radians # TODO: Replace with date specific data

    T_TT = np.linalg.norm((jd_TT - 2451545.0) / 36525)
    Sp = np.deg2rad((-0.0015 * (((a_c ** 2) / 1.2) + (a_a ** 2)) * T_TT) *
                    (1 / 3600))  # Sexag->Radian

    # Initialize rotation matrix [W]
    W = np.matrix(('0 0 0; 0 0 0; 0 0 0'), dtype=np.float64)
    W[0, 0] = (np.cos(Xp) * np.cos(Sp))
    W[0, 1] = ((-np.cos(Yp) * np.sin(Sp)) +
               (np.sin(Yp) * np.sin(Xp) * np.cos(Sp)))
    W[0, 2] = ((-np.sin(Yp) * np.sin(Sp)) -
               (np.cos(Yp) * np.sin(Xp) * np.cos(Sp)))
    W[1, 0] = (np.cos(Xp) * np.sin(Sp))
    W[1, 1] = ((np.cos(Yp) * np.cos(Sp)) +
               (np.sin(Yp) * np.sin(Xp) * np.sin(Sp)))
    W[1, 2] = ((np.sin(Yp) * np.cos(Sp)) -
               (np.cos(Yp) * np.sin(Xp) * np.sin(Sp)))
    W[2, 0] = (np.sin(Xp))
    W[2, 1] = (-np.sin(Yp) * np.cos(Xp))
    W[2, 2] = (np.cos(Yp) * np.cos(Xp))

    # Apply Transform based on Transpose value (Forward or Backward Transform)
    if Transpose == 0:
        rad_tirs = (W * rad_itrf)
        vel_tirs = ((W * vel_itrf))
    elif Transpose == 1:
        rad_tirs = (np.transpose(W) * rad_itrf)
        vel_tirs = (np.transpose(W) * vel_itrf)
    return rad_tirs, vel_tirs


def IAU_ERotationAngle(rad_tirs, vel_tirs, gd_UTC, Transpose):
    """
    Transforms vectors from TIRS to CIRS frame following IAU-2010 conventions

    Parameters
    ----------
    rad_itrf : numpy matrix [3, 1] - [[X], [Y], [Z]]
        - Radius vector components defined in kilometers in the TIRS frame
    vel_itrf : numpy matrix [3, 1] - [[VX], [VY], [VZ]]
        - Velocity vector components defined in kilometers in the TIRS frame
    gd_UTC : numpy matrix [6, 1] - [[Yr], [Mo], [Day], [Hr], [Min], [Sec]]
        - Gregorian Date
    Transpose : int of 0 or 1
        - Determines wether the transform is tirs->cirs (0) or cirs->tirs (1)

    Returns
    -------
    rad_tirs : numpy matrix [3, 1] - [[X], [Y], [Z]]
        - Radius vector components defined in kilometers in the CIRS frame
    vel_tirs : numpy matrix [3, 1] - [[VX], [VY], [VZ]]
        - Velocity vector components defined in kilometers in the CIRS frame

    See Also
    --------
    FK5_ECFixed2J2000 : Transform from an earth fixed frame (ECEF) to an
    earth based inertial frame (J2000) using the IAU-76/FK5 reduction

    FK5_J20002ECFixed : Transform from an earth based inertial frame (J2000)
    to an earth fixed frame (ECEF) using the IAU-76/FK5 reduction

    References
    ----------
    [1] D. Vallado, `Fundamentals of Astrodynamics and Applications`. 4th ed.,
    Microcosm Press, 2013.
        - Pg. 212-213
    """
    # Earth Rotation Angle (IAU-2006/2000, CIO Based) TIRS -> CIRS, Pg. 212

    # Time Adjustments
    jd_UTC, jd_UT1, jd_TAI, jd_TT = TimeAdjust(gd_UTC)

    # Transform Constants
    # These from IERS Earth Orientation Data - EOP 14 C04 (IAU2000A) - Latest
    # https://datacenter.iers.org/eop/-/somos/5Rgv/latest/223
    LOD = (0.0015308)  # sec # TODO: Replace with date specific data

    E_w = np.matrix([[0., 0., (7.292115146706979e-05 * (1 - (LOD / 86400)))]])
    T_era = np.linalg.norm((2 * np.pi) * (0.7790572732640 +
                           1.00273781191135448 * (jd_UT1 - 2451545.0)))  # rad
    # Reduce to within 2 pi radians
    T_ERA = np.mod(T_era, (2 * np.pi))
    # Convert from Seconds to Degrees

    # Initialize rotation matrix [ROT3]
    R = np.matrix(('0 0 0; 0 0 0; 0 0 0'), dtype=np.float64)
    R[0, 0] = np.cos(-T_ERA)
    R[0, 1] = np.sin(-T_ERA)
    R[0, 2] = 0.
    R[1, 0] = -np.sin(-T_ERA)
    R[1, 1] = np.cos(-T_ERA)
    R[1, 2] = 0.
    R[2, 0] = 0.
    R[2, 1] = 0.
    R[2, 2] = 1.

    # Apply Transform based on Transpose value (Forward or Backward Transform)
    if Transpose == 0:
        rad_cirs = (R * rad_tirs)
        vel_cirs = ((R * vel_tirs) +
                    np.transpose(np.cross(E_w, np.transpose(rad_tirs))))
    elif Transpose == 1:
        rad_cirs = (np.transpose(R) * rad_tirs)
        vel_cirs = (np.transpose(R) * vel_tirs -
                    np.transpose(np.cross(E_w, np.transpose(rad_tirs))))
    else:
        raise RuntimeError("Enter an int of 0, or 1 for the reverse transform")
    return rad_cirs, vel_cirs


###############################################################################
###############################################################################
###############################################################################
###############################################################################


def FK5_PolarMotion(rad, vel, gd_UTC):
    # Polar Motion (IAU-76/FK5) - ITRF -> PEF, Vallado pg. 223

    # Import radius and velocity
    rad_itrf = dc(rad)
    vel_itrf = dc(vel)

    # Time Adjustments
    jd_UTC, jd_UT1, jd_TAI, jd_TT = TimeAdjust(gd_UTC)
    # TODO: Need to use to determine constnts below from site at specific times
    # Call JD2Gregorian one you can read the below site

    # These from IERS Earth Orientation Data - EOP 14 C04 (IAU2000A) - Latest
    # https://datacenter.iers.org/eop/-/somos/5Rgv/latest/223

    # Xp_d = (0.100751 * (0.000028 / 0.1))  # Degrees
    # Xp = np.deg2rad(Xp_d)  # Radians
    # Yp_d = (0.447751 * (0.000028 / 0.1))  # Degrees
    # Yp = np.deg2rad(pm_Yp_d)  # Radians

    Xp = (-0.140682 * (1 / 3600) * (np.pi / 180))  # Radians
    Yp = (0.333309 * (1 / 3600) * (np.pi / 180))  # Radians

    ITRF2PEF = np.matrix(('0 0 0; 0 0 0; 0 0 0'), dtype=np.float)
    ITRF2PEF[0, 0] = (1)
    ITRF2PEF[0, 1] = (0)
    ITRF2PEF[0, 2] = (-Xp)
    ITRF2PEF[1, 0] = (0)
    ITRF2PEF[1, 1] = (1)
    ITRF2PEF[1, 2] = (Yp)
    ITRF2PEF[2, 0] = (Xp)
    ITRF2PEF[2, 1] = (-Yp)
    ITRF2PEF[2, 2] = (1)

    rad_pef = (ITRF2PEF * rad_itrf)
    vel_pef = (ITRF2PEF * vel_itrf)
    return rad_pef, vel_pef


def FK5_SiderealTime(rad, vel, gd_UTC):
    # Sidereal Time (IAU-76/FK5) - PEF -> TOD, Vallado pg. 224
    rad_pef = dc(rad)
    vel_pef = dc(vel)

    # Time Adjustments
    jd_UTC, jd_UT1, jd_TAI, jd_TT = TimeAdjust(gd_UTC)
    T_TT = ((jd_TT - 2451545.0) / 36525)
    T_UT1 = ((jd_UT1 - 2451545.0) / 36525)

    theta_s = (67310.54841 +
               (((876600 * 3600) + 8640184.812866) * T_UT1) +
               (0.093104 * (T_UT1 ** 2)) -
               (6.2e-06 * (T_UT1 ** 3)))  # Sec
    # Reduce to within 86400 seconds
    reduc = (theta_s / 86400)
    diff = reduc - np.floor(reduc)
    theta_sec = (86400 * diff)  # Seconds
    # Convert from Seconds to Degrees
    theta_d = (theta_sec * (1 / 240))  # Degrees
    if theta_d < 0:
        theta_gmst = (theta_d + 360)  # Degrees
    if theta_d > 0:
        theta_gmst = theta_d  # Degrees

    Nutation = np.load('orbital_analyses/Nutation.npy')

    # Need to calculate with sum on pg. 226

    M_moon_deg = (134.96298139 + (((1325 * 360) + 198.8673981) * T_TT) +
                  (0.0086972 * (T_TT ** 2)) + (1.78e-05 * (T_TT ** 3)))  # Deg
    # Reduce to within 360 degrees
    M_moon_d = np.mod(M_moon_deg, 360, dtype=np.float64)  # Deg
    M_moon = np.deg2rad(M_moon_d)  # Rad

    M_sun_deg = (357.52772333 + (((99 * 360) + 359.0503400) * T_TT) -
                 (0.0001603 * (T_TT ** 2)) - (3.3e-06 * (T_TT ** 3)))  # Deg
    # Reduce to within 360 degrees
    M_sun_d = np.mod(M_sun_deg, 360, dtype=np.float64)  # Deg
    M_sun = np.deg2rad(M_sun_d)  # Rad

    u_moon_deg = (93.27191028 + (((1342 * 360) + 82.0175381) * T_TT) -
                  (0.0036825 * (T_TT ** 2)) + (3.1e-06 * (T_TT ** 3)))  # Deg
    # Reduce to within 360 degrees
    u_moon_d = np.mod(u_moon_deg, 360, dtype=np.float64)  # Deg
    u_moon = np.deg2rad(u_moon_d)  # Rad

    D_sun_deg = (297.85036306 + (((1236 * 360) + 307.1114800) * T_TT) -
                 (0.0019142 * (T_TT ** 2)) + (5.3e-06 * (T_TT ** 3)))  # Deg
    # Reduce to within 360 degrees
    D_sun_d = np.mod(D_sun_deg, 360, dtype=np.float64)  # Deg
    D_sun = np.deg2rad(D_sun_d)  # Rad

    O_moon_deg = (125.04452222 - (((5 * 360) + 134.1362608) * T_TT) +
                  (0.0020708 * (T_TT ** 2)) + (2.2e-06 * (T_TT ** 3)))  # Deg
    # Reduce to within 360 degrees
    O_moon_d = np.mod(O_moon_deg, 360, dtype=np.float64)  # Deg
    O_moon = np.deg2rad(O_moon_d)  # Rad

    # Initialize summation variable for nutation in longitude
    psi_d = 0

    for i in range(106):
        api1 = ((Nutation[i, 0] * M_moon) + (Nutation[i, 1] * M_sun) +
                (Nutation[i, 2] * u_moon) + (Nutation[i, 3] * D_sun) +
                (Nutation[i, 4] * O_moon))
        psi_d = (psi_d + (((Nutation[i, 5] * (2.777778e-08 / 1)) +
                          ((Nutation[i, 6] * (2.777778e-08 / 1)) *
                           T_TT)) * np.sin(api1)))  # Deg

        # Initialize summation variable for nutation in the obliquity
    st_eps_d = 0

    for j in range(106):
        api2 = ((Nutation[j, 0] * M_moon) + (Nutation[j, 1] * M_sun) +
                (Nutation[j, 2] * u_moon) + (Nutation[j, 3] * D_sun) +
                (Nutation[j, 4] * O_moon))
        st_eps_d = (st_eps_d + (((Nutation[j, 7] * (2.777778e-08 / 1)) +
                                ((Nutation[j, 8] * (2.777778e-08 / 1)) *
                                 T_TT)) * np.cos(api2)))  # Deg

    psi_s = (psi_d * (3600 / 1))  # Sexagesimal

    eps_bar_d = (23.439291 - (0.0130042 * T_TT) - (1.64e-07 * (T_TT ** 2)) +
                 (5.04e-07 * (T_TT ** 3)))  # Deg
    eps_bar = (eps_bar_d * (np.pi / 180))  # Rad

    if jd_UTC < 2450506.500000:
        Eq_s = (psi_s * np.cos(eps_bar))  # Sexagesimal
        Eq_1982 = (Eq_s * (1 / 3600))  # Deg
    if jd_UTC > 2450506.500000:
        Eq_s = ((psi_s * np.cos(eps_bar)) + (0.00264 * np.sin(O_moon)) +
                (0.000063 * np.sin(2 * O_moon)))  # Sexagesimal
        Eq_1982 = (Eq_s * (1 / 3600))  # Deg

    theta_deg = (theta_gmst + Eq_1982)  # Deg
    theta_gast = (theta_deg * (np.pi / 180))  # Rad

    # Earth Angular Velocity
    ang = np.matrix(('0; 0; 0'), dtype=np.float)
    ang[0] = 0
    ang[1] = 0
    ang[2] = 0.0000729211585530

    PEF2TOD = np.matrix(('0 0 0; 0 0 0; 0 0 0'), dtype=np.float)
    PEF2TOD[0, 0] = (np.cos(-theta_gast))
    PEF2TOD[0, 1] = (np.sin(-theta_gast))
    PEF2TOD[0, 2] = (0)
    PEF2TOD[1, 0] = (-np.sin(-theta_gast))
    PEF2TOD[1, 1] = (np.cos(-theta_gast))
    PEF2TOD[1, 2] = (0)
    PEF2TOD[2, 0] = (0)
    PEF2TOD[2, 1] = (0)
    PEF2TOD[2, 2] = (1)

    rad_tod = (PEF2TOD * rad_pef)
    vel_tod = (PEF2TOD * (vel_pef + (np.cross(ang, rad_pef, axis=0))))
    return rad_tod, vel_tod


def FK5_Nutation(rad, vel, gd_UTC):
    # Nutation (IAU-76/FK5) - TOD -> MOD, Vallado pg. 224
    rad_tod = dc(rad)
    vel_tod = dc(vel)

    # Time Adjustments
    jd_UTC, jd_UT1, jd_TAI, jd_TT = TimeAdjust(gd_UTC)
    T_TT = ((jd_TT - 2451545.0) / 36525)

    Nutation = np.load('orbital_analyses/Nutation.npy')

    # These from IERS Earth Orientation Data - EOP 14 C04 (IAU2000A) - Latest
    # https://datacenter.iers.org/eop/-/somos/5Rgv/latest/223
#    nu_delt_psi = (-0.103204 * (0.000028 / 0.1) * (np.pi / 180))  # Radians
#    nu_delt_eps = (-0.012101 * (0.000028 / 0.1) * (np.pi / 180))  # Radians

    # Values are zero if converting to J2000 rather than GCRF
    delt_psi = 0  # Rad
    delt_eps = 0  # Rad

    eps_bar_d = (23.439291 - (0.0130042 * T_TT) - (1.64e-07 * (T_TT ** 2)) +
                 (5.04e-07 * (T_TT ** 3)))  # Deg
    eps_bar = np.deg2rad(eps_bar_d)  # Rad

    # Need to calculate with sum on pg. 226

    M_moon_deg = (134.96298139 + (((1325 * 360) + 198.8673981) * T_TT) +
                  (0.0086972 * (T_TT ** 2)) + (1.78e-05 * (T_TT ** 3)))  # Deg
    # Reduce to within 360 degrees
    M_moon_d = np.mod(M_moon_deg, 360, dtype=np.float64)  # Deg
    M_moon = np.deg2rad(M_moon_d)  # Rad

    M_sun_deg = (357.52772333 + (((99 * 360) + 359.0503400) * T_TT) -
                 (0.0001603 * (T_TT ** 2)) - (3.3e-06 * (T_TT ** 3)))  # Deg
    # Reduce to within 360 degrees
    M_sun_d = np.mod(M_sun_deg, 360, dtype=np.float64)  # Deg
    M_sun = np.deg2rad(M_sun_d)  # Rad

    u_moon_deg = (93.27191028 + (((1342 * 360) + 82.0175381) * T_TT) -
                  (0.0036825 * (T_TT ** 2)) + (3.1e-06 * (T_TT ** 3)))  # Deg
    # Reduce to within 360 degrees
    u_moon_d = np.mod(u_moon_deg, 360, dtype=np.float64)  # Deg
    u_moon = np.deg2rad(u_moon_d)  # Rad

    D_sun_deg = (297.85036306 + (((1236 * 360) + 307.1114800) * T_TT) -
                 (0.0019142 * (T_TT ** 2)) + (5.3e-06 * (T_TT ** 3)))  # Deg
    # Reduce to within 360 degrees
    D_sun_d = np.mod(D_sun_deg, 360, dtype=np.float64)  # Deg
    D_sun = np.deg2rad(D_sun_d)  # Rad

    O_moon_deg = (125.04452222 - (((5 * 360) + 134.1362608) * T_TT) +
                  (0.0020708 * (T_TT ** 2)) + (2.2e-06 * (T_TT ** 3)))  # Deg
    # Reduce to within 360 degrees
    O_moon_d = np.mod(O_moon_deg, 360, dtype=np.float64)  # Deg
    O_moon = np.deg2rad(O_moon_d)  # Rad

    # Initialize summation variable for nutation in longitude
    psi_d = 0

    for i in range(106):
        api1 = ((Nutation[i, 0] * M_moon) + (Nutation[i, 1] * M_sun) +
                (Nutation[i, 2] * u_moon) + (Nutation[i, 3] * D_sun) +
                (Nutation[i, 4] * O_moon))
        psi_d = (psi_d + (((Nutation[i, 5] * (2.777778e-08 / 1)) +
                          ((Nutation[i, 6] * (2.777778e-08 / 1)) *
                           T_TT)) * np.sin(api1)))  # Deg
    nu_psi_r = np.deg2rad(psi_d)  # Rad

    # Initialize summation variable for nutation in the obliquity
    eps_d = 0

    for j in range(106):
        api2 = ((Nutation[j, 0] * M_moon) + (Nutation[j, 1] * M_sun) +
                (Nutation[j, 2] * u_moon) + (Nutation[j, 3] * D_sun) +
                (Nutation[j, 4] * O_moon))
        eps_d = (eps_d + (((Nutation[j, 7] * (2.777778e-08 / 1)) +
                           ((Nutation[j, 8] * (2.777778e-08 / 1)) *
                            T_TT)) * np.cos(api2)))  # Deg

    # Determine values for matrix rotations
    eps_delt = ((eps_d * (np.pi / 180)) + delt_eps)  # Rad
    eps_1980 = (eps_bar + eps_delt)  # Rad
    psi_1980 = (nu_psi_r + delt_psi)  # Rad

    TOD2MOD_3 = np.matrix(('0 0 0; 0 0 0; 0 0 0'), dtype=np.float)
    TOD2MOD_3[0, 0] = (1)
    TOD2MOD_3[0, 1] = (0)
    TOD2MOD_3[0, 2] = (0)
    TOD2MOD_3[1, 0] = (0)
    TOD2MOD_3[1, 1] = (np.cos(eps_bar))
    TOD2MOD_3[1, 2] = (-np.sin(eps_bar))
    TOD2MOD_3[2, 0] = (0)
    TOD2MOD_3[2, 1] = (np.sin(eps_bar))
    TOD2MOD_3[2, 2] = (np.cos(eps_bar))

    TOD2MOD_2 = np.matrix(('0 0 0; 0 0 0; 0 0 0'), dtype=np.float)
    TOD2MOD_2[0, 0] = (np.cos(psi_1980))
    TOD2MOD_2[0, 1] = (np.sin(psi_1980))
    TOD2MOD_2[0, 2] = (0)
    TOD2MOD_2[1, 0] = (-np.sin(psi_1980))
    TOD2MOD_2[1, 1] = (np.cos(psi_1980))
    TOD2MOD_2[1, 2] = (0)
    TOD2MOD_2[2, 0] = (0)
    TOD2MOD_2[2, 1] = (0)
    TOD2MOD_2[2, 2] = (1)

    TOD2MOD_1 = np.matrix(('0 0 0; 0 0 0; 0 0 0'), dtype=np.float)
    TOD2MOD_1[0, 0] = (1)
    TOD2MOD_1[0, 1] = (0)
    TOD2MOD_1[0, 2] = (0)
    TOD2MOD_1[1, 0] = (0)
    TOD2MOD_1[1, 1] = (np.cos(eps_1980))
    TOD2MOD_1[1, 2] = (np.sin(eps_1980))
    TOD2MOD_1[2, 0] = (0)
    TOD2MOD_1[2, 1] = (-np.sin(eps_1980))
    TOD2MOD_1[2, 2] = (np.cos(eps_1980))

    rad_mod = (TOD2MOD_3 * TOD2MOD_2 * TOD2MOD_1 * rad_tod)
    vel_mod = (TOD2MOD_3 * TOD2MOD_2 * TOD2MOD_1 * vel_tod)
    return(rad_mod, vel_mod)


def FK5_Precession(rad, vel, gd_UTC):
    # Precession (IAU-76/FK5) - TOD -> MOD, Vallado pg. 226
    rad_mod = dc(rad)
    vel_mod = dc(vel)

    # Time Adjustments
    jd_UTC, jd_UT1, jd_TAI, jd_TT = TimeAdjust(gd_UTC)
    T_TT = ((jd_TT - 2451545.0) / 36525)

    # Determine values for matrix rotations
    sigma_s = ((2306.2181 * T_TT) + (0.30188 * (T_TT ** 2)) +
               (0.017998 * (T_TT ** 3)))  # Sexagesimal
    theta_s = ((2004.3109 * T_TT) - (0.42665 * (T_TT ** 2)) -
               (0.041833 * (T_TT ** 3)))  # Sexagesimal
    zeta_s = ((2306.2181 * T_TT) + (1.09468 * (T_TT ** 2)) +
              (0.018203 * (T_TT ** 3)))  # Sexagesimal

    sigma_r = (sigma_s * (1 / 3600) * (np.pi / 180))  # Rad
    theta_r = (theta_s * (1 / 3600) * (np.pi / 180))  # Rad
    zeta_r = (zeta_s * (1 / 3600) * (np.pi / 180))  # Rad

    # Set transformation matrix and fill with values
    MOD2GCRF = np.matrix(('0 0 0; 0 0 0; 0 0 0'), dtype=np.float)
    MOD2GCRF[0, 0] = ((np.cos(theta_r) * np.cos(zeta_r) *
                      np.cos(sigma_r)) - (np.sin(zeta_r) * np.sin(sigma_r)))
    MOD2GCRF[0, 1] = ((np.sin(zeta_r) * np.cos(theta_r) *
                      np.cos(sigma_r)) + (np.sin(sigma_r) * np.cos(zeta_r)))
    MOD2GCRF[0, 2] = (np.sin(theta_r) * np.cos(sigma_r))
    MOD2GCRF[1, 0] = ((-np.sin(sigma_r) * np.cos(theta_r) *
                      np.cos(zeta_r)) - (np.sin(zeta_r) * np.cos(sigma_r)))
    MOD2GCRF[1, 1] = ((-np.sin(zeta_r) * np.sin(sigma_r) *
                      np.cos(theta_r)) + (np.cos(zeta_r) * np.cos(sigma_r)))
    MOD2GCRF[1, 2] = (-np.sin(theta_r) * np.sin(sigma_r))
    MOD2GCRF[2, 0] = (-np.sin(theta_r) * np.cos(zeta_r))
    MOD2GCRF[2, 1] = (-np.sin(theta_r) * np.sin(zeta_r))
    MOD2GCRF[2, 2] = (np.cos(theta_r))

    rad_gcrf = (MOD2GCRF * rad_mod)
    vel_gcrf = (MOD2GCRF * vel_mod)

    return(rad_gcrf, vel_gcrf)


def FK5_ECFixed2J2000(rad, vel, gd_UTC):
    # Entire IAU-76/FK5 Reduction from ITFR -> GCRF/J2000
    fk_rad1, fk_vel1 = FK5_PolarMotion(rad, vel, gd_UTC)
    fk_rad2, fk_vel2 = FK5_SiderealTime(fk_rad1, fk_vel1, gd_UTC)
    fk_rad3, fk_vel3 = FK5_Nutation(fk_rad2, fk_vel2, gd_UTC)
    fk_rad4, fk_vel4 = FK5_Precession(fk_rad3, fk_vel3, gd_UTC)
    return fk_rad4, fk_vel4

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################


def FK5_Precession_T(rad, vel, gd_UTC):
    # Precession (IAU-76/FK5) - TOD -> MOD, Vallado pg. 226
    rad_mod = dc(rad)
    vel_mod = dc(vel)

    # Time Adjustments
    jd_UTC, jd_UT1, jd_TAI, jd_TT = TimeAdjust(gd_UTC)
    T_TT = ((jd_TT - 2451545.0) / 36525)

    # Determine values for matrix rotations
    sigma_s = ((2306.2181 * T_TT) + (0.30188 * (T_TT ** 2)) +
               (0.017998 * (T_TT ** 3)))  # Sexagesimal
    theta_s = ((2004.3109 * T_TT) - (0.42665 * (T_TT ** 2)) -
               (0.041833 * (T_TT ** 3)))  # Sexagesimal
    zeta_s = ((2306.2181 * T_TT) + (1.09468 * (T_TT ** 2)) +
              (0.018203 * (T_TT ** 3)))  # Sexagesimal

    sigma_r = (sigma_s * (1 / 3600) * (np.pi / 180))  # Radians
    theta_r = (theta_s * (1 / 3600) * (np.pi / 180))  # Radians
    zeta_r = (zeta_s * (1 / 3600) * (np.pi / 180))  # Radians

    # Set transformation matrix and fill with values
    PREC = np.matrix(('0 0 0; 0 0 0; 0 0 0'), dtype=np.float)
    PREC[0, 0] = ((np.cos(theta_r) * np.cos(zeta_r) *
                   np.cos(sigma_r)) - (np.sin(zeta_r) * np.sin(sigma_r)))
    PREC[0, 1] = ((np.sin(zeta_r) * np.cos(theta_r) *
                   np.cos(sigma_r)) + (np.sin(sigma_r) * np.cos(zeta_r)))
    PREC[0, 2] = (np.sin(theta_r) * np.cos(sigma_r))
    PREC[1, 0] = ((-np.sin(sigma_r) * np.cos(theta_r) *
                   np.cos(zeta_r)) - (np.sin(zeta_r) * np.cos(sigma_r)))
    PREC[1, 1] = ((-np.sin(zeta_r) * np.sin(sigma_r) *
                   np.cos(theta_r)) + (np.cos(zeta_r) * np.cos(sigma_r)))
    PREC[1, 2] = (-np.sin(theta_r) * np.sin(sigma_r))
    PREC[2, 0] = (-np.sin(theta_r) * np.cos(zeta_r))
    PREC[2, 1] = (-np.sin(theta_r) * np.sin(zeta_r))
    PREC[2, 2] = (np.cos(theta_r))

    PREC_T = np.transpose(PREC)
    rad_gcrf = (PREC_T * rad_mod)
    vel_gcrf = (PREC_T * vel_mod)
    return(rad_gcrf, vel_gcrf)


def FK5_Nutation_T(rad, vel, gd_UTC):
    # Nutation (IAU-76/FK5) - TOD -> MOD, Vallado pg. 224
    rad_tod = dc(rad)
    vel_tod = dc(vel)

    # Time Adjustments
    jd_UTC, jd_UT1, jd_TAI, jd_TT = TimeAdjust(gd_UTC)
    T_TT = ((jd_TT - 2451545.0) / 36525)

    Nutation = np.load('orbital_analyses/Nutation.npy')

    # These from IERS Earth Orientation Data - EOP 14 C04 (IAU2000A) - Latest
    # https://datacenter.iers.org/eop/-/somos/5Rgv/latest/223
#    nu_delt_psi = (-0.052242 * (0.000028 / 0.1) * (np.pi / 180))  # Rad
#    nu_delt_eps = (-0.003951 * (0.000028 / 0.1) * (np.pi / 180))  # Rad

    # Values are zero if converting to J2000 rather than GCRF
    delt_psi = 0  # Rad
    delt_eps = 0  # Rad

    eps_bar_d = (23.439291 - (0.0130042 * T_TT) - (1.64e-07 * (T_TT ** 2)) +
                 (5.04e-07 * (T_TT ** 3)))  # Deg
    eps_bar = np.deg2rad(eps_bar_d)  # Rad

    # Need to calculate with sum on pg. 226

    M_moon_deg = (134.96298139 + (((1325 * 360) + 198.8673981) * T_TT) +
                  (0.0086972 * (T_TT ** 2)) + (1.78e-05 * (T_TT ** 3)))  # Deg
    # Reduce to within 360 degrees
    M_moon_d = np.mod(M_moon_deg, 360, dtype=np.float64)  # Deg
    M_moon = np.deg2rad(M_moon_d)  # Rad

    M_sun_deg = (357.52772333 + (((99 * 360) + 359.0503400) * T_TT) -
                 (0.0001603 * (T_TT ** 2)) - (3.3e-06 * (T_TT ** 3)))  # Deg
    # Reduce to within 360 degrees
    M_sun_d = np.mod(M_sun_deg, 360, dtype=np.float64)  # Deg
    M_sun = np.deg2rad(M_sun_d)  # Rad

    u_moon_deg = (93.27191028 + (((1342 * 360) + 82.0175381) * T_TT) -
                  (0.0036825 * (T_TT ** 2)) + (3.1e-06 * (T_TT ** 3)))  # Deg
    # Reduce to within 360 degrees
    u_moon_d = np.mod(u_moon_deg, 360, dtype=np.float64)  # Deg
    u_moon = np.deg2rad(u_moon_d)  # Rad

    D_sun_deg = (297.85036306 + (((1236 * 360) + 307.1114800) * T_TT) -
                 (0.0019142 * (T_TT ** 2)) + (5.3e-06 * (T_TT ** 3)))  # Deg
    # Reduce to within 360 degrees
    D_sun_d = np.mod(D_sun_deg, 360, dtype=np.float64)  # Deg
    D_sun = np.deg2rad(D_sun_d)  # Rad

    O_moon_deg = (125.04452222 - (((5 * 360) + 134.1362608) * T_TT) +
                  (0.0020708 * (T_TT ** 2)) + (2.2e-06 * (T_TT ** 3)))  # Deg
    # Reduce to within 360 degrees
    O_moon_d = np.mod(O_moon_deg, 360, dtype=np.float64)  # Deg
    O_moon = np.deg2rad(O_moon_d)  # Rad

    # Initialize summation variable for nutation in longitude
    psi_d = 0

    for i in range(106):
        api1 = ((Nutation[i, 0] * M_moon) + (Nutation[i, 1] * M_sun) +
                (Nutation[i, 2] * u_moon) + (Nutation[i, 3] * D_sun) +
                (Nutation[i, 4] * O_moon))
        psi_d = (psi_d + (((Nutation[i, 5] * (2.777778e-08 / 1)) +
                           ((Nutation[i, 6] * (2.777778e-08 / 1)) *
                            T_TT)) * np.sin(api1)))  # Deg
    psi_r = np.deg2rad(psi_d)  # Rad

    # Initialize summation variable for nutation in the obliquity
    eps_d = 0

    for j in range(106):
        api2 = ((Nutation[j, 0] * M_moon) + (Nutation[j, 1] * M_sun) +
                (Nutation[j, 2] * u_moon) + (Nutation[j, 3] * D_sun) +
                (Nutation[j, 4] * O_moon))
        eps_d = (eps_d + (((Nutation[j, 7] * (2.777778e-08 / 1)) +
                           ((Nutation[j, 8] * (2.777778e-08 / 1)) *
                            T_TT)) * np.cos(api2)))  # Deg

    # Determine values for matrix rotations
    eps_delt = ((eps_d * (np.pi / 180)) + delt_eps)  # Rad
    eps_1980 = (eps_bar + eps_delt)  # Rad
    psi_1980 = (psi_r + delt_psi)  # Rad

    NUT_3 = np.matrix(('0 0 0; 0 0 0; 0 0 0'), dtype=np.float)
    NUT_3[0, 0] = (1)
    NUT_3[0, 1] = (0)
    NUT_3[0, 2] = (0)
    NUT_3[1, 0] = (0)
    NUT_3[1, 1] = (np.cos(eps_bar))
    NUT_3[1, 2] = (-np.sin(eps_bar))
    NUT_3[2, 0] = (0)
    NUT_3[2, 1] = (np.sin(eps_bar))
    NUT_3[2, 2] = (np.cos(eps_bar))

    NUT_2 = np.matrix(('0 0 0; 0 0 0; 0 0 0'), dtype=np.float)
    NUT_2[0, 0] = (np.cos(psi_1980))
    NUT_2[0, 1] = (np.sin(psi_1980))
    NUT_2[0, 2] = (0)
    NUT_2[1, 0] = (-np.sin(psi_1980))
    NUT_2[1, 1] = (np.cos(psi_1980))
    NUT_2[1, 2] = (0)
    NUT_2[2, 0] = (0)
    NUT_2[2, 1] = (0)
    NUT_2[2, 2] = (1)

    NUT_1 = np.matrix(('0 0 0; 0 0 0; 0 0 0'), dtype=np.float)
    NUT_1[0, 0] = (1)
    NUT_1[0, 1] = (0)
    NUT_1[0, 2] = (0)
    NUT_1[1, 0] = (0)
    NUT_1[1, 1] = (np.cos(eps_1980))
    NUT_1[1, 2] = (np.sin(eps_1980))
    NUT_1[2, 0] = (0)
    NUT_1[2, 1] = (-np.sin(eps_1980))
    NUT_1[2, 2] = (np.cos(eps_1980))

    NUT_T = np.transpose(NUT_1 * NUT_2 * NUT_3)
    rad_mod = (NUT_T * rad_tod)
    vel_mod = (NUT_T * vel_tod)
    return(rad_mod, vel_mod)


def FK5_SiderealTime_T(rad, vel, gd_UTC):
    # Sidereal Time (IAU-76/FK5) - PEF -> TOD, Vallado pg. 224
    rad_pef = dc(rad)
    vel_pef = dc(vel)

    # Time Adjustments
    jd_UTC, jd_UT1, jd_TAI, jd_TT = TimeAdjust(gd_UTC)
    T_TT = ((jd_TT - 2451545.0) / 36525)
    T_UT1 = ((jd_UT1 - 2451545.0) / 36525)

    theta_s = (67310.54841 + (((876600 * 3600) + 8640184.812866) * T_UT1) +
               (0.093104 * (T_UT1 ** 2)) - (6.2e-06 * (T_UT1 ** 3)))  # Sec
    # Reduce to within 86400 seconds
    reduc = (theta_s / 86400)
    diff = reduc - np.floor(reduc)
    theta_sec = (86400 * diff)  # Sec
    # Convert from Seconds to Degrees
    theta_d = (theta_sec * (1 / 240))  # Deg
    if theta_d < 0:
        theta_gmst = (theta_d + 360)  # Deg
    if theta_d > 0:
        theta_gmst = theta_d  # Deg

    Nutation = np.load('orbital_analyses/Nutation.npy')

    # Need to calculate with sum on pg. 226

    M_moon_deg = (134.96298139 + (((1325 * 360) + 198.8673981) * T_TT) +
                  (0.0086972 * (T_TT ** 2)) + (1.78e-05 * (T_TT ** 3)))  # Deg
    # Reduce to within 360 degrees
    M_moon_d = np.mod(M_moon_deg, 360, dtype=np.float64)  # Deg
    M_moon = np.deg2rad(M_moon_d)  # Rad

    M_sun_deg = (357.52772333 + (((99 * 360) + 359.0503400) * T_TT) -
                 (0.0001603 * (T_TT ** 2)) - (3.3e-06 * (T_TT ** 3)))  # Deg
    # Reduce to within 360 degrees
    M_sun_d = np.mod(M_sun_deg, 360, dtype=np.float64)  # Deg
    M_sun = np.deg2rad(M_sun_d)  # Rad

    u_moon_deg = (93.27191028 + (((1342 * 360) + 82.0175381) * T_TT) -
                  (0.0036825 * (T_TT ** 2)) + (3.1e-06 * (T_TT ** 3)))  # Deg
    # Reduce to within 360 degrees
    u_moon_d = np.mod(u_moon_deg, 360, dtype=np.float64)  # Deg
    u_moon = np.deg2rad(u_moon_d)  # Rad

    D_sun_deg = (297.85036306 + (((1236 * 360) + 307.1114800) * T_TT) -
                 (0.0019142 * (T_TT ** 2)) + (5.3e-06 * (T_TT ** 3)))  # Deg
    # Reduce to within 360 degrees
    D_sun_d = np.mod(D_sun_deg, 360, dtype=np.float64)  # Deg
    D_sun = np.deg2rad(D_sun_d)  # Rad

    O_moon_deg = (125.04452222 - (((5 * 360) + 134.1362608) * T_TT) +
                  (0.0020708 * (T_TT ** 2)) + (2.2e-06 * (T_TT ** 3)))  # Deg
    # Reduce to within 360 degrees
    O_moon_d = np.mod(O_moon_deg, 360, dtype=np.float64)  # Deg
    O_moon = np.deg2rad(O_moon_d)  # Rad

    # Initialize summation variable for nutation in longitude
    psi_d = 0

    for i in range(106):
        api1 = ((Nutation[i, 0] * M_moon) + (Nutation[i, 1] * M_sun) +
                (Nutation[i, 2] * u_moon) + (Nutation[i, 3] * D_sun) +
                (Nutation[i, 4] * O_moon))
        psi_d = (psi_d + (((Nutation[i, 5] * (2.777778e-08 / 1)) +
                          ((Nutation[i, 6] * (2.777778e-08 / 1)) *
                           T_TT)) * np.sin(api1)))  # Deg

        # Initialize summation variable for nutation in the obliquity
    eps_d = 0

    for j in range(106):
        api2 = ((Nutation[j, 0] * M_moon) + (Nutation[j, 1] * M_sun) +
                (Nutation[j, 2] * u_moon) + (Nutation[j, 3] * D_sun) +
                (Nutation[j, 4] * O_moon))
        eps_d = (eps_d + (((Nutation[j, 7] * (2.777778e-08 / 1)) +
                          ((Nutation[j, 8] * (2.777778e-08 / 1)) *
                           T_TT)) * np.cos(api2)))  # Deg

    psi_s = (psi_d * (3600 / 1))  # Sexagesimal

    eps_bar_d = (23.439291 - (0.0130042 * T_TT) - (1.64e-07 * (T_TT ** 2)) +
                 (5.04e-07 * (T_TT ** 3)))  # Deg
    eps_bar = np.deg2rad(eps_bar_d)  # Rad

    if jd_UTC < 2450506.500000:
        Eq_s = (psi_s * np.cos(eps_bar))  # Sexagesimal
        Eq_1982 = (Eq_s * (1 / 3600))  # Deg
    if jd_UTC > 2450506.500000:
        Eq_s = ((psi_s * np.cos(eps_bar)) + (0.00264 * np.sin(O_moon)) +
                (0.000063 * np.sin(2 * O_moon)))  # Sexagesimal
        Eq_1982 = (Eq_s * (1 / 3600))  # Deg

    theta_deg = (theta_gmst + Eq_1982)  # Deg
    theta_gast = np.deg2rad(theta_deg)  # Rad

    # Earth Angular Velocity
    ang = np.matrix(([0.], [0.], [0.0000729211585530]), dtype=np.float)

    ST = np.matrix(('0 0 0; 0 0 0; 0 0 0'), dtype=np.float)
    ST[0, 0] = (np.cos(-theta_gast))
    ST[0, 1] = (np.sin(-theta_gast))
    ST[0, 2] = (0)
    ST[1, 0] = (-np.sin(-theta_gast))
    ST[1, 1] = (np.cos(-theta_gast))
    ST[1, 2] = (0)
    ST[2, 0] = (0)
    ST[2, 1] = (0)
    ST[2, 2] = (1)

    ST_T = np.transpose(ST)
    rad_tod = (ST_T * rad_pef)
    vel_tod = (ST_T * (vel_pef + (np.cross(ang, rad_pef, axis=0))))
    return (rad_tod, vel_tod)


def FK5_PolarMotion_T(rad, vel, gd_UTC):
    # Polar Motion (IAU-76/FK5) - ITRF -> PEF, Vallado pg. 223
    rad_itrf = dc(rad)
    vel_itrf = dc(vel)

    # Time Adjustments
    jd_UTC, jd_UT1, jd_TAI, jd_TT = TimeAdjust(gd_UTC)
    # Need to use to determine constnts below from site at specific times
    # Call JD2Gregorian one you can read the below site

    # These from IERS Earth Orientation Data - EOP 14 C04 (IAU2000A) - Latest
    # https://datacenter.iers.org/eop/-/somos/5Rgv/latest/223

    # Xp_d = (0.100751 * (0.000028 / 0.1))  # Deg
    # Xp = np.deg2rad(pm_Xp_d)  # Rad
    # Yp_d = (0.447751 * (0.000028 / 0.1))  # Deg
    # Yp = np.deg2rad(pm_Yp_d)  # Rad

    Xp = (-0.140865 * (1 / 3600) * (np.pi / 180))  # Rad
    Yp = (0.333414 * (1 / 3600) * (np.pi / 180))  # Rad

    PM = np.matrix(('0 0 0; 0 0 0; 0 0 0'), dtype=np.float)
    PM[0, 0] = (1)
    PM[0, 1] = (0)
    PM[0, 2] = (-Xp)
    PM[1, 0] = (0)
    PM[1, 1] = (1)
    PM[1, 2] = (Yp)
    PM[2, 0] = (Xp)
    PM[2, 1] = (-Yp)
    PM[2, 2] = (1)

    PM_T = np.transpose(PM)
    rad_pef = (PM_T * rad_itrf)
    vel_pef = (PM_T * vel_itrf)

    return (rad_pef, vel_pef)


# Radius Transform is not good
# Velocity Transform is not good
def FK5_J20002ECFixed(rad, vel, gd_UTC):
    # Entire IAU-76/FK5 Reduction from GCRF/J2000 -> ITFR (ECEF)
    fk_rad1, fk_vel1 = FK5_Precession_T(rad, vel, gd_UTC)  # GOOD
    fk_rad2, fk_vel2 = FK5_Nutation_T(fk_rad1, fk_vel1, gd_UTC)
    fk_rad3, fk_vel3 = FK5_SiderealTime_T(fk_rad2, fk_vel2, gd_UTC)
    fk_rad4, fk_vel4 = FK5_PolarMotion_T(fk_rad3, fk_vel3, gd_UTC)
    return fk_rad4, fk_vel4  # Results in about 200 m error across X, Y, Z


# =============================================================================
# import numpy as np
# from Transform_Coordinate import FK5_Precession_T
# from Transform_Coordinate import FK5_Nutation_T
# from Transform_Coordinate import FK5_SiderealTime_T
# from Transform_Coordinate import FK5_PolarMotion_T
#
# rad = np.matrix([[-3112.262153644471],
#                  [7512.28484046077],
#                  [3745.374697833162]])
# vel = np.matrix([[-5.468896292796726],
#                  [-2.664564053439451],
#                  [3.67784024548162]])
#
# GD_UTC = np.zeros((1, 6))
# GD_UTC[0, 0] = (2004)
# GD_UTC[0, 1] = (4)
# GD_UTC[0, 2] = (6)
# GD_UTC[0, 3] = (0)
# GD_UTC[0, 4] = (0)
# GD_UTC[0, 5] = (0)
# gd_UTC = np.matrix([[GD_UTC[0, 0]], [GD_UTC[0, 1]], [GD_UTC[0, 2]],
#                     [GD_UTC[0, 3]], [GD_UTC[0, 4]], [GD_UTC[0, 5]]])
# =============================================================================
