# Utilized Modules
import numpy as np
from copy import deepcopy as dc
import os
from urllib.request import urlopen
from bs4 import BeautifulSoup
import re
from datetime import datetime, timedelta

# ALL ASSUMED TO BE IN Stated Coordinate Systems
# Must convert to appropriate coordinates before using these functions
# Use functions in Transform_Coordinate.py for that


def JD2Gregorian(JD):
    """
    Converts Julian Date to Gregorian Date Format.

    Parameters
    ----------
    JD : array_like [1, n] or float / int input
        - Can handle array_like of n Julian Dates as input
        - Does not handle ModJulian (Modified Julian Date)

    Returns
    -------
    GD : numpy matrix [6, n] - [[Yr], [Mo], [Day], [Hr], [Min], [Sec]]
        - Returns values as floats due to Sec having decimal values

    See Also
    --------
    Gregorian2JD : The reverse of this conversion

    References
    ----------
    [1] D. Vallado, `Fundamentals of Astrodynamics and Applications`. 4th ed.,
    Microcosm Press, 2013.
        - Alg. 22, pg. 202
    """
    # Initialize Vectors
    JD = np.asmatrix(JD)
    length = np.size(JD)
    GD = np.zeros((6, length), dtype=float)

    for i in range(length):
        T_1900 = ((JD[0, i] - 2415019.5) / 365.25)
        # Determine year and deal with any leap year issues
        GD[0, i] = (1900 + np.trunc(T_1900))
        leap = (np.trunc((GD[0, i] - 1900 - 1) * 0.25))
        day = ((JD[0, i] - 2415019.5) - (((GD[0, i] - 1900) * 365) + leap))
        if day < 1.0:
            GD[0, i] = (GD[0, i] - 1)
            leap = (np.trunc((GD[0, i] - 1900 - 1) * 0.25))
            day = ((JD[0, i] - 2415019.5) - (((GD[0, i] - 1900) * 365) + leap))
        # Matrix of month lengths
        ML = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
        # Account for extra day in Feb in leap years
        if np.mod(GD[0, i], 4) == 0:
            ML[1] = 29
        dayC = np.trunc(day)
        # Determine month
        mon_add = 0
        for k in range(12):
            if mon_add < dayC:
                mon_add = (mon_add + ML[k])
            if mon_add >= dayC:
                GD[1, i] = (k + 1)
                mon_add = (mon_add - ML[k])
                break
        # Determine Day
        GD[2, i] = (dayC - mon_add)
        # Determine Hour
        tau = ((day - dayC) * 24)
        GD[3, i] = np.trunc(tau)
        # Determine Minute
        GD[4, i] = np.trunc((tau - GD[3, i]) * 60)
        # Determine Second
        GD[5, i] = ((tau - GD[3, i] - (GD[4, i] / 60)) * 3600)
    return np.asmatrix(GD)


def Gregorian2JD(GD):
    """
    Converts Julian Date to Gregorian Date Format.
        - Accurate for dates between year 1900 - 2100

    Parameters
    ----------
    GD : array_like [6, n]
        - Can handle a numpy matrix of n Gregorian Dates as an input

    Returns
    -------
    JD : numpy matrix [1, n]
        - Returns values as floats due to Sec having decimal values

    See Also
    --------
    JD2Gregorian : The reverse of this conversion

    References
    ----------
    [1] D. Vallado, `Fundamentals of Astrodynamics and Applications`. 4th ed.,
    Microcosm Press, 2013.
        - Alg. 14, pg. 183
    """
    # Initialize Vectors
    GD = np.matrix(GD)
    length = np.size(GD, axis=1)
    JD = np.zeros((length, 1))

    for j in range(length):
        gj_1 = (367 * GD[0, j])
        gj_2 = np.int((7 * (GD[0, j] + np.int((GD[1, j] + 9) / 12))) / 4)
        gj_3 = np.int((275 * GD[1, j]) / 9)
        gj_4 = (GD[2, j])
        gj_5 = (((((GD[5, j] / 60) + GD[4, j]) / 60) + GD[3, j]) / 24)

        JD[j] = np.asmatrix(np.float64((gj_1) - (gj_2) + (gj_3) +
                            (gj_4) + (1721013.5) + (gj_5)))
    return JD


def TimeAdjust(GD):
    """
    Converts UTC Gregorian Date to UT1, TAI, and TT timeframes.
    Outputs are all Julian Dates
        - Input timeframe must be UTC to maintain accuracy

    Parameters
    ----------
    GD : array_like [6, 1]
        - Single Gregorian Date input. Timeframe in UTC

    Returns
    -------
    JD_UTC : float - Julian Date

    JD_UT1 : float - Julian Date

    JD_TAI : float - Julian Date

    JD_TT : float - Julian Date

    References
    ----------
    [1] D. Vallado, `Fundamentals of Astrodynamics and Applications`. 4th ed.,
    Microcosm Press, 2013.
    Alg. 16, pg. 195
    """
    # Get Modified Julian Date from input GD
    JD_UTC = np.asscalar(Gregorian2JD(np.asmatrix(GD)))
    MJD_UTC = (JD_UTC - 2400000.5)
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

    # Scrape IERS Data to get dAT
    if os.path.exists(r'orbital_analyses\deltaTA.npy'):
        deltaTA = np.load(r'orbital_analyses\deltaTA.npy')
    elif os.path.exists(r'deltaTA.npy'):
        deltaTA = np.load(r'deltaTA.npy')
    dTA_index = np.searchsorted(np.ravel(deltaTA[:, 3]), JD_UTC)

    # Check if date is exactly or greater than index value and edit if needed
    if EOP_index == np.size(EOPCO4[:, 3]):
        EOP_index = (EOP_index - 1)
    elif MJD_UTC != EOPCO4[EOP_index, 3]:
        EOP_index = (EOP_index - 1)

    if dTA_index == np.size(deltaTA[:, 3]):
        dTA_index = (dTA_index - 1)
    elif JD_UTC != deltaTA[dTA_index, 3]:
        dTA_index = (dTA_index - 1)

    # TODO: Currently date outside data range takes last value in data, change?
    dUT1 = np.asscalar(EOPCO4[EOP_index, 6])  # Seconds
    dAT = np.asscalar(deltaTA[dTA_index, 4])  # Seconds
    dTT = (32.184)  # Seconds  # TODO: Does this change?
    deltas = np.array([dUT1, dAT, (dAT + dTT)])

    # Add time to GD_UTC
    t_ms = np.asmatrix(np.zeros((7, 1), dtype=np.float64))
    t_ms[0:5, 0] = GD[0:5, 0]
    t_ms[5, 0] = np.floor(GD[5, 0])  # Seconds
    t_ms[6, 0] = (np.mod(GD[5, 0], 1) * 1e6)  # Milliseconds
    t_ms = t_ms.astype(int)
    GD_utc = datetime(t_ms[0, 0], t_ms[1, 0], t_ms[2, 0], t_ms[3, 0],
                      t_ms[4, 0], t_ms[5, 0], t_ms[6, 0])  # To datetime
    # Calculate time delta to adjust time from TT to TDB
    JD_deltas = np.zeros((3), dtype=float)
    for i in range(np.size(deltas)):
        add_sec = deltas[i]
        GD_add = timedelta(seconds=add_sec)
        GD_f = GD_utc + GD_add
        # Convert back to numpy array
        GD_np = np.datetime64(GD_f)
        GD_new = np.asmatrix(np.zeros((6, 1), dtype=np.float64))
        GD_new[0, 0] = GD_np.astype(object).year
        GD_new[1, 0] = GD_np.astype(object).month
        GD_new[2, 0] = GD_np.astype(object).day
        GD_new[3, 0] = GD_np.astype(object).hour
        GD_new[4, 0] = GD_np.astype(object).minute
        GD_new[5, 0] = (GD_np.astype(object).second +
                        (GD_np.astype(object).microsecond * 1e-6))
        JD_deltas[i] = np.asscalar(Gregorian2JD(GD_new))

    return JD_UTC, JD_deltas[0], JD_deltas[1], JD_deltas[2]


def Keplerian2Perifocal(kep):
    """
    Converts Keplerian orbital elements to Rad / Vel in the Pericocal System

    Parameters
    ----------
    kep : array_like [6, n]
        - Can handle a numpy matrix of n Keplerian Elements as an input

    Returns
    -------
    rad : numpy matrix [3, n]

    vel : numpy matrix [3, n]

    See Also
    --------
    Perifocal2Keplerian : The reverse of this conversion
        - In development

    References
    ----------
    [1] D. Vallado, `Fundamentals of Astrodynamics and Applications`. 4th ed.,
    Microcosm Press, 2013.
        - Part of Alg. 10, pg. 118
    """
    # Initialize Vectors
    kep = np.asmatrix(kep)
    if np.size(kep, axis=0) != 6:
        kep = np.transpose(kep)
    length = np.size(kep, axis=1)

    rad = np.asmatrix(np.zeros((3, length), dtype=np.float64))
    vel = np.asmatrix(np.zeros((3, length), dtype=np.float64))
    E_Mu = 3.986004418e5

    for k in range(length):
        # Determine vriables for special case orbits
        lop = kep[3, k] + kep[4, k]             # omega true
        aol = kep[4, k] + kep[5, k]             # u
        tl = kep[3, k] + kep[4, k] + kep[5, k]  # lambda true
        if (kep[1, k] < 1e-6 and kep[2, k] == 0):
            kep[3, k], kep[4, k], kep[5, k] = 0, 0, tl
        if (kep[1, k] < 1e-6 and kep[2, k] != 0):
            kep[4, k], kep[5, k] = 0, aol
        if (kep[1, k] > 1e-6 and kep[2, k] == 0):
            kep[3, k], kep[4, k] = 0, lop

        # Convert from degrees to radians
        kep[2:, k] = np.deg2rad(kep[2:, k])

        # Solve for Perifocal System Coordinates
        rad[0, k] = ((kep[0, k] * np.cos(kep[5, k])) /
                     (1 + kep[1, k] * np.cos(kep[5, k])))
        rad[1, k] = ((kep[0, k] * np.sin(kep[5, k])) /
                     (1 + kep[1, k] * np.cos(kep[5, k])))
        rad[2, k] = 0.
        vel[0, k] = (-(np.sqrt(E_Mu / kep[0, k])) * np.sin(kep[5, k]))
        vel[1, k] = ((np.sqrt(E_Mu / kep[0, k])) *
                     (kep[1, k] + np.cos(kep[5, k])))
        vel[2, k] = 0.
    return rad, vel


def Perifocal2Keplerian(rad, vel):
    return


def Keplerian2Cartesian(kep):
    """
    Converts Keplerian orbital elements to Rad / Vel in the Cartesian System

    Parameters
    ----------
    kep : array_like [6, n]
        - Can handle a numpy matrix of n Keplerian Elements as an input

    Returns
    -------
    rad : numpy matrix [3, n]

    vel : numpy matrix [3, n]

    See Also
    --------
    Cartesian2Keplerian : The reverse of this conversion

    References
    ----------
    [1] D. Vallado, `Fundamentals of Astrodynamics and Applications`. 4th ed.,
    Microcosm Press, 2013.
        - Alg. 10, pg. 118-119
    """
    # Initialize Vectors
    if np.size(kep, axis=0) != 6:
        kep = np.asmatrix(np.transpose(dc(kep)))
    else:
        kep = np.asmatrix(dc(kep))

    length = np.size(kep, axis=1)
    E_Mu = 3.986004418e5

    r = np.asmatrix(np.zeros((3, length), dtype=np.float64))
    v = np.asmatrix(np.zeros((3, length), dtype=np.float64))
    rad = np.asmatrix(np.zeros((3, length), dtype=np.float64))
    vel = np.asmatrix(np.zeros((3, length), dtype=np.float64))

    for l in range(length):
        # Convert sma to paramater
        kep[0, l] = (kep[0, l] * (1 - (kep[1, l] ** 2)))

        # Determine vriables for special case orbits
        lop = (dc(kep[3, l]) + dc(kep[4, l]))                   # omega true
        aol = (dc(kep[4, l]) + dc(kep[5, l]))                   # u
        tl = (dc(kep[3, l]) + dc(kep[4, l]) + dc(kep[5, l]))    # lambda true
        if (kep[1, l] < 1e-10 and kep[2, l] == 0):
            kep[3, l], kep[4, l], kep[5, l] = 0, 0, dc(tl)
        if (kep[1, l] < 1e-10 and kep[2, l] != 0):
            kep[4, l], kep[5, l] = 0, dc(aol)
        if (kep[1, l] > 1e-10 and kep[2, l] == 0):
            kep[3, l], kep[4, l] = 0, dc(lop)

        # Convert from degrees to radians
        kep[2:, l] = np.deg2rad(kep[2:, l])

        # Solve for Perifocal System Coordinates
        r[0] = ((kep[0, l] * np.cos(kep[5, l])) /
                (1 + kep[1, l] * np.cos(kep[5, l])))
        r[1] = ((kep[0, l] * np.sin(kep[5, l])) /
                (1 + kep[1, l] * np.cos(kep[5, l])))
        r[2] = 0.
        v[0] = (-(np.sqrt(E_Mu / kep[0, l])) * np.sin(kep[5, l]))
        v[1] = ((np.sqrt(E_Mu / kep[0, l])) * (kep[1, l] + np.cos(kep[5, l])))
        v[2] = 0.

        # Rotate Perifocal System Coordinates to Geocentric Equitorial System
        IJK2PQW = np.matrix(('0 0 0; 0 0 0; 0 0 0'), dtype=np.float)
        IJK2PQW[0, 0] = ((np.cos(kep[3, l]) * np.cos(kep[4, l])) -
                         (np.sin(kep[3, l]) * np.sin(kep[4, l]) *
                          np.cos(kep[2, l])))
        IJK2PQW[1, 0] = ((np.sin(kep[3, l]) * np.cos(kep[4, l])) +
                         (np.cos(kep[3, l]) * np.sin(kep[4, l]) *
                          np.cos(kep[2, l])))
        IJK2PQW[2, 0] = (np.sin(kep[4, l]) * np.sin(kep[2, l]))
        IJK2PQW[0, 1] = (-(np.cos(kep[3, l]) * np.sin(kep[4, l])) -
                         (np.sin(kep[3, l]) * np.cos(kep[4, l]) *
                          np.cos(kep[2, l])))
        IJK2PQW[1, 1] = (-(np.sin(kep[3, l]) * np.sin(kep[4, l])) +
                         (np.cos(kep[3, l]) * np.cos(kep[4, l]) *
                          np.cos(kep[2, l])))
        IJK2PQW[2, 1] = (np.cos(kep[4, l]) * np.sin(kep[2, l]))
        IJK2PQW[0, 2] = (np.sin(kep[3, l]) * np.sin(kep[2, l]))
        IJK2PQW[1, 2] = (-(np.cos(kep[3, l]) * np.sin(kep[2, l])))
        IJK2PQW[2, 2] = (np.cos(kep[2, l]))

        # Transform coordinates
        rad[:, l] = (IJK2PQW * r)
        vel[:, l] = (IJK2PQW * v)
    return rad, vel


def Cartesian2Keplerian(rad, vel):
    """
    Converts Rad / Vel in the Cartesian System to Keplerian orbital elements

    Parameters
    ----------
    rad : array_like [3, n]
        - Can handle a numpy matrix of n radius vectors as an input
    vel : array_like [3, n]
        - Can handle a numpy matrix of n velocity vectors as an input

    Returns
    -------
    kep : numpy matrix [6, n]

    See Also
    --------
    Keplerian2Cartesian : The reverse of this conversion

    References
    ----------
    [1] D. Vallado, `Fundamentals of Astrodynamics and Applications`. 4th ed.,
    Microcosm Press, 2013.
        - Alg. 9, pg. 113-114
    """
    # Initialize Vectors
    rad = np.asmatrix(rad)
    if np.size(rad, axis=0) != 3:
        rad = np.transpose(rad)
    vel = np.asmatrix(vel)
    if np.size(vel, axis=0) != 3:
        vel = np.transpose(vel)

    length = np.size(rad, axis=1)
    # From Fundamentals of Astrydynamics & Applications: Vallado - pg.113-114
    E_Mu = 3.986004418e5
    kep = np.asmatrix(np.zeros((6, length), dtype=np.float64))

    for m in range(length):
        h_v = np.cross(np.transpose(rad[:, m]), np.transpose(vel[:, m]))
        h = np.linalg.norm(h_v)
        n_v = np.cross([0, 0, 1], h_v)
        n = np.linalg.norm(n_v)
        r = np.linalg.norm(rad[:, m])
        v = np.linalg.norm(vel[:, m])
        ecc_v = (((((v ** 2) - (E_Mu / r)) * rad[:, m]) -
                  (vel[:, m] * np.dot((np.transpose(rad[:, m])),
                                      (vel[:, m])))) / E_Mu)
        kep[1, m] = np.linalg.norm(ecc_v)  # Unitless
        me = (((v ** 2) / 2) - (E_Mu / r))
        kep[0, m] = (-(E_Mu / (2 * me)))  # km
        kep[2, m] = np.arccos(h_v[0, 2] / h) * (180 / np.pi)  # deg

        # Check quadrant resolution and undefined values for angled values
        kep[3, m] = np.arccos(n_v[0, 0] / n) * (180 / np.pi)  # deg
        if np.isnan(kep[3, m]):
            kep[3, m] = 0
        if n_v[0, 1] < 0:
            kep[3, m] = ((360) - kep[3, m])

        kep[4, m] = np.arccos(np.dot(n_v, ecc_v) /
                              (n * kep[1, m])) * (180 / np.pi)  # deg
        if np.isnan(kep[4, m]):
            kep[4, m] = 0
        if ecc_v[2, 0] < 0:
            kep[4, m] = ((360) - kep[4, m])

        kep[5, m] = np.arccos(np.dot(np.transpose(ecc_v), rad[:, m]) /
                              (kep[1, m] * r)) * (180 / np.pi)  # deg
        if np.dot(np.transpose(rad[:, m]), vel[:, m]) < 0:
            kep[5, m] = ((360) - kep[5, m])
    return kep


def Cartesian2Perifocal(rad, vel):
    return


def Perifocal2Cartesian(rad, vel):
    return
