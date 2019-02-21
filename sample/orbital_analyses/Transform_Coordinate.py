# Utilized Modules
import numpy as np
from copy import deepcopy as dc
from datetime import datetime, timedelta
import os
from urllib.request import urlopen
from bs4 import BeautifulSoup
import re

from sample.orbital_analyses.Transform_State import convert_time
from sample.orbital_analyses.Transform_State import JD2Gregorian
from sample.orbital_analyses.Transform_State import gregorian2julian_date

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


# TODO: add handing for acceleration transforms
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


# TODO: add handing for acceleration transforms
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


# TODO: add handing for acceleration transforms
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


# TODO: add handing for acceleration transforms
def IAU_PolarMotion(rad_itrf, vel_itrf, gd_UTC, Transpose):
    """
    Transform vectors between ITRF to TIRS frame following IAU-2010 conventions

    Parameters
    ----------
    rad_itrf : numpy matrix [3, 1] - [[X], [Y], [Z]]
        - Radius vector components defined in kilometers in the ITRF frame
    vel_itrf : numpy matrix [3, 1] - [[VX], [VY], [VZ]]
        - Velocity vector components defined in kilometers in the ITRF frame
    gd_UTC : numpy matrix [6, 1] - [[Yr], [Mo], [Day], [Hr], [Min], [Sec]]
        - Gregorian Date
    Transpose : int of 0 or 1
        - Determines wether the transform is ITRF->TIRS (0) or TIRS->ITRF (1)

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
    jd_UTC, jd_UT1, jd_TAI, jd_TT = convert_time(gd_UTC)

    # Get Modified Julian Date from input GD
    MJD_UTC = (jd_UTC - 2400000.5)
    # Scrape IERS Data to get dUT1
    if os.path.exists(r'orbital_analyses\EOPCO4.npy'):
        EOPCO4 = np.load(r'orbital_analyses\EOPCO4.npy')
    elif os.path.exists(r'EOPCO4.npy'):
        EOPCO4 = np.load(r'EOPCO4.npy')
    else:
        EOP_scrape = "https://datacenter.iers.org/eop/-/somos/5Rgv/latest/224"  # TODO: Check to see if correct final data table
        EOP_page = urlopen(EOP_scrape)
        EOP_soup = BeautifulSoup(EOP_page, "lxml")
        EOP_string = str(EOP_soup.body.p.string)
        EOP_list = re.split(r'\n+', EOP_string.rstrip('\n'))
        del EOP_list[:10]  # Delete initial description lines
        EOPCO4 = np.matrix(np.zeros((np.size(EOP_list), 16), dtype=float))
        for n in range(np.size(EOP_list)):  # convert to numpy matrix
            EOPCO4[n, :] = np.fromstring(EOP_list[n], dtype=float, sep=" ")
        np.save("EOPCO4.npy", EOPCO4)
    EOP_index = np.searchsorted(np.ravel(EOPCO4[:, 3]), MJD_UTC)

    # Check if date is exactly or greater than index value and edit if needed
    if EOP_index == np.size(EOPCO4[:, 3]):
        EOP_index = (EOP_index - 1)
    elif MJD_UTC != EOPCO4[EOP_index, 3]:
        EOP_index = (EOP_index - 1)

    # Transform Constants
    # From Ast Almonac, 2006?:B76
    a_a = (0.12)  # Arcseconds  # TODO: Replace with date specific data
    a_c = (0.26)  # Arcseconds  # TODO: Replace with date specific data

    # These from IERS Earth Orientation Data - EOP 14 C04 (IAU2000A) - Latest
    Xp = np.deg2rad((np.asscalar(EOPCO4[EOP_index, 4])) * (1 / 3600))
    Yp = np.deg2rad((np.asscalar(EOPCO4[EOP_index, 5])) * (1 / 3600))

    T_TT = np.linalg.norm((jd_TT - 2451545.0) / 36525)
    Sp = np.deg2rad((-0.0015 * (((a_c ** 2) / 1.2) + (a_a ** 2)) * T_TT) *
                    (1 / 3600))  # Arcseconds->Radian

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


# TODO: add handing for acceleration transforms
def IAU_ERotationAngle(rad_tirs, vel_tirs, gd_UTC, Transpose):
    """
    Transform vectors between TIRS & CIRS frame following IAU-2010 conventions

    Parameters
    ----------
    rad_tirs : numpy matrix [3, 1] - [[X], [Y], [Z]]
        - Radius vector components defined in kilometers in the TIRS frame
    vel_tirs : numpy matrix [3, 1] - [[VX], [VY], [VZ]]
        - Velocity vector components defined in kilometers in the TIRS frame
    gd_UTC : numpy matrix [6, 1] - [[Yr], [Mo], [Day], [Hr], [Min], [Sec]]
        - Gregorian Date
    Transpose : int of 0 or 1
        - Determines wether the transform is TIRS->CIRS (0) or CIRS->TIRS (1)

    Returns
    -------
    rad_cirs : numpy matrix [3, 1] - [[X], [Y], [Z]]
        - Radius vector components defined in kilometers in the CIRS frame
    vel_cirs : numpy matrix [3, 1] - [[VX], [VY], [VZ]]
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

    # Raise errors before running script
    while Transpose not in [0, 1]:
        raise RuntimeError("Enter an int of 0; or 1 for the reverse transform")

    # Time Adjustments
    jd_UTC, jd_UT1, jd_TAI, jd_TT = convert_time(gd_UTC)

    # Get Modified Julian Date from input GD
    MJD_UTC = (jd_UTC - 2400000.5)
    # Scrape IERS Data to get dUT1
    if os.path.exists(r'orbital_analyses\EOPCO4.npy'):
        EOPCO4 = np.load(r'orbital_analyses\EOPCO4.npy')
    elif os.path.exists(r'EOPCO4.npy'):
        EOPCO4 = np.load(r'EOPCO4.npy')
    else:
        EOP_scrape = "https://datacenter.iers.org/eop/-/somos/5Rgv/latest/224"  # TODO: Check to see if correct final data table
        EOP_page = urlopen(EOP_scrape)
        EOP_soup = BeautifulSoup(EOP_page, "lxml")
        EOP_string = str(EOP_soup.body.p.string)
        EOP_list = re.split(r'\n+', EOP_string.rstrip('\n'))
        del EOP_list[:10]  # Delete initial description lines
        EOPCO4 = np.matrix(np.zeros((np.size(EOP_list), 16), dtype=float))
        for n in range(np.size(EOP_list)):  # convert to numpy matrix
            EOPCO4[n, :] = np.fromstring(EOP_list[n], dtype=float, sep=" ")
        np.save("EOPCO4.npy", EOPCO4)
    EOP_index = np.searchsorted(np.ravel(EOPCO4[:, 3]), MJD_UTC)

    # Check if date is exactly or greater than index value and edit if needed
    if EOP_index == np.size(EOPCO4[:, 3]):
        EOP_index = (EOP_index - 1)
    elif MJD_UTC != EOPCO4[EOP_index, 3]:
        EOP_index = (EOP_index - 1)

    # This from IERS Earth Orientation Data - EOP 14 C04 (IAU2000A) - Latest
    LOD = np.asscalar(EOPCO4[EOP_index, 7])  # sec

    E_w = np.matrix([[0., 0., (7.292115146706979e-05 * (1 - (LOD / 86400)))]])
    T_era = np.linalg.norm((2 * np.pi) * (0.7790572732640 +
                           1.00273781191135448 * (jd_UT1 - 2451545.0)))  # rad
    # Reduce to within 2 pi radians
    T_ERA = np.mod(T_era, (2 * np.pi))

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
        vel_cirs = (R * (vel_tirs +
                    np.transpose(np.cross(E_w, np.transpose(rad_tirs)))))
    elif Transpose == 1:
        rad_cirs = (np.transpose(R) * rad_tirs)
        vel_cirs = (np.transpose(R) * (vel_tirs -
                    np.transpose(np.cross(E_w, np.transpose(rad_tirs)))))
    return rad_cirs, vel_cirs


# TODO: add handing for acceleration transforms
def IAU_PrecessionNutation(rad_cirs, vel_cirs, gd_UTC, Transpose):
    """
    Transform vectors between CIRS & GCRF frame following IAU-2010 conventions

    Parameters
    ----------
    rad_cirs : numpy matrix [3, 1] - [[X], [Y], [Z]]
        - Radius vector components defined in kilometers in the CIRS frame
    vel_cirs : numpy matrix [3, 1] - [[VX], [VY], [VZ]]
        - Velocity vector components defined in kilometers in the CIRS frame
    gd_UTC : numpy matrix [6, 1] - [[Yr], [Mo], [Day], [Hr], [Min], [Sec]]
        - Gregorian Date
    Transpose : int of 0 or 1
        - Determines wether the transform is CIRS->GCRF (0) or GCRF->CIRS (1)

    Returns
    -------
    rad_gcrf : numpy matrix [3, 1] - [[X], [Y], [Z]]
        - Radius vector components defined in kilometers in the GCRF frame
    vel_gcrf : numpy matrix [3, 1] - [[VX], [VY], [VZ]]
        - Velocity vector components defined in kilometers in the GCRF frame

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
        - Pg. 213-215, 1045
    [2] IAU 2006/2000A expression of the X coordinate of the CIP in the GCRS
        - http://maia.usno.navy.mil/conventions/2010/2010_official/chapter5/additional_info/tab5.2a.txt
    [3] IAU 2006/2000A expression of the Y coordinate of the CIP in the GCRS
        - http://maia.usno.navy.mil/conventions/2010/2010_official/chapter5/additional_info/tab5.2b.txt
    [4] IAU 2006/2000A expression of the quantity s(t) + XY/2
        - http://maia.usno.navy.mil/conventions/2010/2010_official/chapter5/additional_info/tab5.2d.txt
    """
    # Precession-Nutation (IAU-2006/2000, CIO Based) CIRS -> GCRF, Pg. 213

    # Raise errors before running script
    while Transpose not in [0, 1]:
        raise RuntimeError("Enter an int of 0; or 1 for the reverse transform")

    # Time Adjustments
    jd_UTC, jd_UT1, jd_TAI, jd_TT = convert_time(gd_UTC)
    T_TT = np.linalg.norm((jd_TT - 2451545.0) / 36525)
    gd_TT = JD2Gregorian(jd_TT)
    # Initial time for TDB determination
    t_ms = np.asmatrix(np.zeros((7, 1), dtype=np.float64))
    t_ms[0:5, 0] = gd_TT[0:5, 0]
    t_ms[5, 0] = np.floor(gd_TT[5, 0])
    t_ms[6, 0] = (np.mod(gd_TT[5, 0], 1) * 1e6)
    t_ms = t_ms.astype(int)
    TDB_i = datetime(t_ms[0, 0], t_ms[1, 0], t_ms[2, 0],
                     t_ms[3, 0], t_ms[4, 0], t_ms[5, 0], t_ms[6, 0])
    # Calculate time delta to adjust time from TT to TDB
    # TODO: Double check this time conversion
    add_sec = (0.001657 * np.sin(628.3076 * T_TT + 6.2401) +
               (0.000022 * np.sin(575.3385 * T_TT + 4.2970)) +
               (0.000014 * np.sin(1256.6152 * T_TT + 6.1969)) +
               (0.000005 * np.sin(606.9777 * T_TT + 4.0212)) +
               (0.000005 * np.sin(52.9691 * T_TT + 0.4444)) +
               (0.000002 * np.sin(21.3299 * T_TT + 5.5431)) +
               (0.000010 * T_TT * np.sin(628.3076 * T_TT + 4.2490)))  # sec
    TDB_add = timedelta(seconds=add_sec)
    TDB_f = TDB_i + TDB_add
    # Convert back to numpy array
    TDB_np = np.datetime64(TDB_f)
    gd_TDB = np.asmatrix(np.zeros((6, 1), dtype=np.float64))
    gd_TDB[0, 0] = TDB_np.astype(object).year
    gd_TDB[1, 0] = TDB_np.astype(object).month
    gd_TDB[2, 0] = TDB_np.astype(object).day
    gd_TDB[3, 0] = TDB_np.astype(object).hour
    gd_TDB[4, 0] = TDB_np.astype(object).minute
    gd_TDB[5, 0] = (TDB_np.astype(object).second +
                    (TDB_np.astype(object).microsecond * 1e-6))
    # Convert to JD
    jd_TDB = gregorian2julian_date(gd_TDB)
    # Determine julian centuries for TDB
    T_TDB = np.linalg.norm((jd_TDB - 2451545.0) / 36525)
    # Additional data from online
    # These from IERS Earth Orientation Data - EOP 14 C04 (IAU2000A) - Latest
    # https://datacenter.iers.org/eop/-/somos/5Rgv/latest/223
    dX = -0.000205  # sexagesimal
    dY = -0.000136  # sexagesimal

    # Transform Constants
    # Determine angles for Earth's Nutation (r = 360 deg)
    Nu_term = np.zeros((14, 1), dtype=np.float)
    Nu_term[0, 0] = ((485868.249036 + (1717915923.2178 * T_TT) +
                      (31.8792 * (T_TT ** 2)) + (0.051635 * (T_TT ** 3)) -
                      (0.00024470 * (T_TT ** 4))) * (1 / 3600))  # deg
    Nu_term[1, 0] = ((1287104.79305 + (129596581.0481 * T_TT) -
                      (0.5532 * (T_TT ** 2)) + (0.000136 * (T_TT ** 3)) -
                      (0.00001149 * (T_TT ** 4))) * (1 / 3600))  # deg
    Nu_term[2, 0] = ((335779.526232 + (1739527262.8478 * T_TT) -
                      (12.7512 * (T_TT ** 2)) - (0.001037 * (T_TT ** 3)) +
                      (0.00000417 * (T_TT ** 4))) * (1 / 3600))  # deg
    Nu_term[3, 0] = ((1072260.70369 + (1602961601.2090 * T_TT) -
                      (6.3706 * (T_TT ** 2)) + (0.006593 * (T_TT ** 3)) -
                      (0.00003169 * (T_TT ** 4))) * (1 / 3600))  # deg
    Nu_term[4, 0] = ((450160.398036 - (6962890.5431 * T_TT) +
                      (7.4722 * (T_TT ** 2)) + (0.007702 * (T_TT ** 3)) -
                      (0.00005939 * (T_TT ** 4))) * (1 / 3600))  # deg
    # Determine planetary nutation values
    Nu_term[5, 0] = (252.250905494 + (149472.6746358 * T_TDB))  # deg
    Nu_term[6, 0] = (181.979800853 + (58517.8156748 * T_TDB))  # deg
    Nu_term[7, 0] = (100.466448494 + (35999.3728521 * T_TDB))  # deg
    Nu_term[8, 0] = (355.433274605 + (19140.299314 * T_TDB))  # deg
    Nu_term[9, 0] = (34.351483900 + (3034.90567464 * T_TDB))  # deg
    Nu_term[10, 0] = (50.0774713998 + (1222.11379404 * T_TDB))  # deg
    Nu_term[11, 0] = (314.055005137 + (428.466998313 * T_TDB))  # deg
    Nu_term[12, 0] = (304.348665499 + (218.486200208 * T_TDB))  # deg
    Nu_term[13, 0] = ((1.39697137214 * T_TDB) + (0.0003086 * (T_TDB ** 2)))
    # Reduce all values to withing 360 deg
    for j in range(np.size(Nu_term)):
        Nu_term[j, 0] = np.mod(Nu_term[j, 0], (360))  # deg

    # Import constant matrices
    Xcon = np.load('orbital_analyses/maia_tab5.2a.npy')
    Ycon = np.load('orbital_analyses/maia_tab5.2b.npy')
    Scon = np.load('orbital_analyses/maia_tab5.2d.npy')

    # Calculate all necessary a_p coefficients
    a_pX = np.asmatrix(np.zeros((np.size(Xcon, axis=0), 1), dtype=np.float))
    a_pY = np.asmatrix(np.zeros((np.size(Ycon, axis=0), 1), dtype=np.float))
    a_pS = np.asmatrix(np.zeros((np.size(Scon, axis=0), 1), dtype=np.float))

    # Calculate the X, Y, & s arguments for sin & cos
    for i in range(np.size(Xcon, axis=0)):  # X paramater
        a_pX[i, 0] = sum(x * y for x, y in zip(Nu_term[:, 0], Xcon[i, 3:17]))
        # Reduce withing 360 deg and convert to radians
        a_pX[i, 0] = np.deg2rad(np.mod(a_pX[i, 0], (360)))  # rad
    for i in range(np.size(Ycon, axis=0)):  # Y paramater
        a_pY[i, 0] = sum(x * y for x, y in zip(Nu_term[:, 0], Ycon[i, 3:17]))
        # Reduce withing 360 deg and convert to radians
        a_pY[i, 0] = np.deg2rad(np.mod(a_pY[i, 0], (360)))  # rad
    for i in range(np.size(Scon, axis=0)):  # s paramater
        a_pS[i, 0] = sum(x * y for x, y in zip(Nu_term[:, 0], Scon[i, 3:17]))
        # Reduce withing 360 deg and convert to radians
        a_pS[i, 0] = np.deg2rad(np.mod(a_pS[i, 0], (360)))  # rad

    # X paramater
    Xsum = np.asmatrix(np.zeros((np.size(Xcon, axis=0), 1), dtype=np.float))
    x = 0
    y = 0
    for j in [1306, 253, 36, 4, 1]:
        Xint = 0
        for i in range(0, j):
            Xint = Xint + (((Xcon[i + x, 1] * 1e-06 * np.sin(a_pX[i + x, 0])) +
                           (Xcon[i + x, 2] * 1e-06 * np.cos(a_pX[i + x, 0]))) *
                           (T_TT ** y))
        Xsum[i + x, 0] = Xint
        x += j
        y += 1
    # Y paramater
    Ysum = np.asmatrix(np.zeros((np.size(Ycon, axis=0), 1), dtype=np.float))
    x = 0
    y = 0
    for j in [962, 277, 30, 5, 1]:
        Yint = 0
        for i in range(0, j):
            Yint = Yint + (((Ycon[i + x, 1] * 1e-06 * np.sin(a_pY[i + x, 0])) +
                           (Ycon[i + x, 2] * 1e-06 * np.cos(a_pY[i + x, 0]))) *
                           (T_TT ** y))
        Ysum[i + x, 0] = Yint
        x += j
        y += 1
    # s paramater
    Ssum = np.asmatrix(np.zeros((np.size(Scon, axis=0), 1), dtype=np.float))
    x = 0
    y = 0
    for j in [33, 3, 25, 4, 1]:
        Sint = 0
        for i in range(0, j):
            Sint = Sint + (((Scon[i + x, 1] * 1e-06 * np.sin(a_pS[i + x, 0])) +
                           (Scon[i + x, 2] * 1e-06 * np.cos(a_pS[i + x, 0]))) *
                           (T_TT ** y))
        Ssum[i + x, 0] = Sint
        x += j
        y += 1

    # X, Y, & s final calculations
    X = (-0.016617 + (2004.191898 * T_TT) - (0.4297829 * (T_TT ** 2)) -
         (0.19861834 * (T_TT ** 3)) + (0.000007578 * (T_TT ** 4)) +
         (0.0000059285 * (T_TT ** 5)) + np.sum(Xsum[0:1306, 0]) +
         np.sum(Xsum[1306:1559, 0]) + np.sum(Xsum[1559:1595, 0]) +
         np.sum(Xsum[1595:1599, 0]) + np.sum(Xsum[1599:1600, 0]) +
         dX)  # arcseconds
    Y = (-0.006951 - (0.025896 * T_TT) - (22.4072747 * (T_TT ** 2)) +
         (0.00190059 * (T_TT ** 3)) + (0.001112526 * (T_TT ** 4)) +
         (0.0000001358 * (T_TT ** 5)) + np.sum(Ysum[0:962, 0]) +
         np.sum(Ysum[962:1239, 0]) + np.sum(Ysum[1239:1269, 0]) +
         np.sum(Ysum[1269:1274, 0]) + np.sum(Ysum[1274:1275, 0]) +
         dY)  # arcseconds
    s = (-((X * Y * 1e-06) / 2) + 0.000094 + (0.00380865 * T_TT) -
         (0.00012268 * (T_TT ** 2)) - (0.07257411 * (T_TT ** 3)) +
         (0.00002798 * (T_TT ** 4)) + (0.00001562 * (T_TT ** 5)) +
         np.sum(Ssum[0:33, 0]) + np.sum(Ssum[33:36, 0]) +
         np.sum(Ssum[36:61, 0]) + np.sum(Ssum[61:65, 0]) +
         np.sum(Ssum[65:66, 0]))  # arcseconds

    # Determine a from X and Y
    Xr = np.deg2rad(X * (1 / 3600))  # Rad
    Yr = np.deg2rad(Y * (1 / 3600))  # Rad
    d = np.arctan(np.sqrt(((Xr ** 2) + (Yr ** 2)) /
                          (1 - (Xr ** 2) - (Yr ** 2))))  # deg
    a = np.rad2deg(1 / (1 + np.cos(np.deg2rad(d))))  # deg
    # Convert s and a to radians before plugging into matrices
    ar = np.deg2rad(a)  # rad
    sr = np.deg2rad(s * (1 / 3600))  # rad

    # Initialize rotation matrix [P]
    P = np.matrix(('0 0 0; 0 0 0; 0 0 0'), dtype=np.float64)
    P[0, 0] = (1 - ar * (Xr ** 2))
    P[0, 1] = (-ar * Xr * Yr)
    P[0, 2] = (Xr)
    P[1, 0] = (-ar * Xr * Yr)
    P[1, 1] = (1 - ar * (Yr ** 2))
    P[1, 2] = (Yr)
    P[2, 0] = (-Xr)
    P[2, 1] = (-Yr)
    P[2, 2] = (1 - ar * ((Yr ** 2) + (Xr ** 2)))
    # Initialize rotation matrix [N (ROT3)]
    N = np.matrix(('0 0 0; 0 0 0; 0 0 0'), dtype=np.float64)
    N[0, 0] = np.cos(sr)
    N[0, 1] = np.sin(sr)
    N[0, 2] = 0.
    N[1, 0] = -np.sin(sr)
    N[1, 1] = np.cos(sr)
    N[1, 2] = 0.
    N[2, 0] = 0.
    N[2, 1] = 0.
    N[2, 2] = 1.

    # Apply Transform based on Transpose value (Forward or Backward Transform)
    if Transpose == 0:
        rad_gcrf = (P * N * rad_cirs)
        vel_gcrf = (P * N * vel_cirs)
    elif Transpose == 1:
        # TODO: Check this transform
        rad_gcrf = (np.transpose(P) * np.transpose(N) * rad_cirs)
        vel_gcrf = (np.transpose(P) * np.transpose(N) * vel_cirs)
    return rad_gcrf, vel_gcrf


# TODO: add handing for acceleration transforms
def IAU_2000Reduction(rad, vel, gd_UTC, Transpose):
    """
    Transform vectors between ITRF & GCRF frame following IAU-2010 conventions

    Parameters
    ----------
    rad : numpy matrix [3, 1] - [[X], [Y], [Z]]
        - Radius vector defined in kilometers in the ITRF or GCRF frame
    vel : numpy matrix [3, 1] - [[VX], [VY], [VZ]]
        - Velocity vector defined in kilometers in the ITRF or GCRF frame
    gd_UTC : numpy matrix [6, 1] - [[Yr], [Mo], [Day], [Hr], [Min], [Sec]]
        - Gregorian Date
    Transpose : int of 0 or 1
        - Determines wether the transform is ITRF->GCRF (0) or GCRF->ITRF (1)

    Returns
    -------
    rad : numpy matrix [3, 1] - [[X], [Y], [Z]]
        - Radius vector defined in kilometers in the ITRF or GCRF frame
        depending upon value of "Transpose"
    vel : numpy matrix [3, 1] - [[VX], [VY], [VZ]]
        - Velocity vector defined in kilometers in the ITRF or GCRF frame
        depending upon value of "Transpose"

    See Also
    --------
    FK5_ECFixed2J2000 : Transform from an earth fixed frame (ECEF) to an
    earth based inertial frame (J2000) using the IAU-76/FK5 reduction

    FK5_J20002ECFixed : Transform from an earth based inertial frame (J2000)
    to an earth fixed frame (ECEF) using the IAU-76/FK5 reduction
    """
    # Entire IAU-2000 Reduction from ITRF -> GCRF/J2000
    if Transpose == 0:
        rad1, vel1 = IAU_PolarMotion(rad, vel, gd_UTC, Transpose)
        rad2, vel2 = IAU_ERotationAngle(rad1, vel1, gd_UTC, Transpose)
        rad3, vel3 = IAU_PrecessionNutation(rad2, vel2, gd_UTC, Transpose)
    elif Transpose == 1:
        rad1, vel1 = IAU_PrecessionNutation(rad, vel, gd_UTC, Transpose)
        rad2, vel2 = IAU_ERotationAngle(rad1, vel1, gd_UTC, Transpose)
        rad3, vel3 = IAU_PolarMotion(rad2, vel2, gd_UTC, Transpose)
    return rad3, vel3

###############################################################################
###############################################################################


# TODO: add handing for acceleration transforms
def FK5_PolarMotion(rad, vel, gd_UTC):
    # Polar Motion (IAU-76/FK5) - ITRF -> PEF, Vallado pg. 223

    # Import radius and velocity
    rad_itrf = dc(rad)
    vel_itrf = dc(vel)

    # Time Adjustments
    jd_UTC, jd_UT1, jd_TAI, jd_TT = convert_time(gd_UTC)

    # Get Modified Julian Date from input GD
    MJD_UTC = (jd_UTC - 2400000.5)
    # Scrape IERS Data to get dUT1
    if os.path.exists(r'orbital_analyses\EOPCO4.npy'):
        EOPCO4 = np.load(r'orbital_analyses\EOPCO4.npy')
    elif os.path.exists(r'EOPCO4.npy'):
        EOPCO4 = np.load(r'EOPCO4.npy')
    else:
        EOP_scrape = "https://datacenter.iers.org/eop/-/somos/5Rgv/latest/224"  # TODO: Check to see if correct final data table
        EOP_page = urlopen(EOP_scrape)
        EOP_soup = BeautifulSoup(EOP_page, "lxml")
        EOP_string = str(EOP_soup.body.p.string)
        EOP_list = re.split(r'\n+', EOP_string.rstrip('\n'))
        del EOP_list[:10]  # Delete initial description lines
        EOPCO4 = np.matrix(np.zeros((np.size(EOP_list), 16), dtype=float))
        for n in range(np.size(EOP_list)):  # convert to numpy matrix
            EOPCO4[n, :] = np.fromstring(EOP_list[n], dtype=float, sep=" ")
        np.save("EOPCO4.npy", EOPCO4)
    EOP_index = np.searchsorted(np.ravel(EOPCO4[:, 3]), MJD_UTC)

    # Check if date is exactly or greater than index value and edit if needed
    if EOP_index == np.size(EOPCO4[:, 3]):
        EOP_index = (EOP_index - 1)
    elif MJD_UTC != EOPCO4[EOP_index, 3]:
        EOP_index = (EOP_index - 1)

    # These from IERS Earth Orientation Data - EOP 14 C04 (IAU2000A) - Latest
    Xp = np.deg2rad((np.asscalar(EOPCO4[EOP_index, 4])) * (1 / 3600))
    Yp = np.deg2rad((np.asscalar(EOPCO4[EOP_index, 5])) * (1 / 3600))

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


# TODO: add handing for acceleration transforms
def FK5_SiderealTime(rad, vel, gd_UTC):
    # Sidereal Time (IAU-76/FK5) - PEF -> TOD, Vallado pg. 224
    rad_pef = dc(rad)
    vel_pef = dc(vel)

    # Time Adjustments
    jd_UTC, jd_UT1, jd_TAI, jd_TT = convert_time(gd_UTC)
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


# TODO: add handing for acceleration transforms
def FK5_Nutation(rad, vel, gd_UTC):
    # Nutation (IAU-76/FK5) - TOD -> MOD, Vallado pg. 224
    rad_tod = dc(rad)
    vel_tod = dc(vel)

    # Time Adjustments
    jd_UTC, jd_UT1, jd_TAI, jd_TT = convert_time(gd_UTC)
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


# TODO: add handing for acceleration transforms
def FK5_Precession(rad, vel, gd_UTC):
    # Precession (IAU-76/FK5) - TOD -> MOD, Vallado pg. 226
    rad_mod = dc(rad)
    vel_mod = dc(vel)

    # Time Adjustments
    jd_UTC, jd_UT1, jd_TAI, jd_TT = convert_time(gd_UTC)
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


# TODO: add handing for acceleration transforms
    # TODO: make like other transform series and remove all duplicate code below
def FK5_ECFixed2J2000(rad, vel, gd_UTC):
    # Entire IAU-76/FK5 Reduction from ITRF -> GCRF/J2000
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
    jd_UTC, jd_UT1, jd_TAI, jd_TT = convert_time(gd_UTC)
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
    jd_UTC, jd_UT1, jd_TAI, jd_TT = convert_time(gd_UTC)
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
    jd_UTC, jd_UT1, jd_TAI, jd_TT = convert_time(gd_UTC)
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
    jd_UTC, jd_UT1, jd_TAI, jd_TT = convert_time(gd_UTC)
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
