# Utilized Modules
import os
import re
import numpy as np
from bs4 import BeautifulSoup
from copy import deepcopy as dc
from urllib.request import urlopen, Request
from datetime import datetime, timedelta
from sample.decorators import gregorian_date_validation

# ALL ASSUMED TO BE IN Stated Coordinate Systems
# Must convert to appropriate coordinates before using these functions
# Use functions in Transform_Coordinate.py for that


# %% GOOD on arrays
@gregorian_date_validation
def increment_gregorian_date(gregorian_date, delta):
    """
    Increments an array of gregorian date(s) by the delta(s) provided
        - Multiple deltas can be inputted

    Parameters
    ----------
    gregorian_date : array_like (n, 6)
        - Imput an numpy ndarray of n Gregorian Dates
        - Arranged in the format np.array([year, month, day, hour, min, sec])
    delta : array_like (n,)
        - Input an numpy ndarray of n time deltas
        - Time deltas must be given in seconds, either positive or negative

    Returns
    -------
    incremented_gregorian_date : tuple (n)
        - Returns a tuple of size n, n defined by number of gregorian dates
        - For each entry, an (m, 6) array, where m is defined by deltas given
    """
    # Determine number of date entries in gregorian_date
    gd = np.array(gregorian_date)
    try:
        gd.shape[1]
    except IndexError:
        entries = 1
    else:
        entries = gd.shape[0]
    # Determine number of delta values to evaluate for each gregorian date
    # TODO: add size validation here?
    deltas = np.array(delta).shape[0]

    # Initialize output tuple
    incremented_gregorian_date = []
    # Loop to iterate through gregorian dates
    for i in range(entries):
        # Handling if only one gregorian date supplied
        if entries == 1:
            # Add time to GD_UTC
            time_ms = np.zeros((7, 1), dtype=float)
            time_ms[0:5, 0] = gd[0:5]
            time_ms[5, 0] = np.floor(gd[5])  # Seconds
            time_ms[6, 0] = (np.mod(gd[5], 1) * 1e6)  # Milliseconds
            time_ms = time_ms.astype(int)
            gd_UTC = datetime(time_ms[0, 0], time_ms[1, 0], time_ms[2, 0],
                              time_ms[3, 0], time_ms[4, 0], time_ms[5, 0],
                              time_ms[6, 0])  # Convert to datetime
            # Using time delta, initilize array to put final gd values into
            gd_deltas = np.zeros((deltas, 6), dtype=float)
            # Loop to iterate through deltas
        else:
            # Add time to GD_UTC
            time_ms = np.zeros((7, 1), dtype=float)
            time_ms[0:5, 0] = gd[i, 0:5]
            time_ms[5, 0] = np.floor(gd[i, 5])  # Seconds
            time_ms[6, 0] = (np.mod(gd[i, 5], 1) * 1e6)  # Milliseconds
            time_ms = time_ms.astype(int)
            gd_UTC = datetime(time_ms[0, 0], time_ms[1, 0], time_ms[2, 0],
                              time_ms[3, 0], time_ms[4, 0], time_ms[5, 0],
                              time_ms[6, 0])  # Convert to datetime
            # Using time delta, initilize array to put final gd values into
            gd_deltas = np.zeros((deltas, 6), dtype=float)
            # Loop to iterate through deltas
        for j in range(deltas):
            # Add seconds from each user defined delta and increment the date
            add_sec = delta[j]
            gd_add = timedelta(seconds=add_sec)
            gd_final = gd_UTC + gd_add
            # Convert back to numpy array
            gd_np = np.datetime64(gd_final)
            gd_deltas[j, 0] = gd_np.astype(object).year
            gd_deltas[j, 1] = gd_np.astype(object).month
            gd_deltas[j, 2] = gd_np.astype(object).day
            gd_deltas[j, 3] = gd_np.astype(object).hour
            gd_deltas[j, 4] = gd_np.astype(object).minute
            gd_deltas[j, 5] = (gd_np.astype(object).second +
                               (gd_np.astype(object).microsecond * 1e-6))
        incremented_gregorian_date.append(gd_deltas)

    return tuple(incremented_gregorian_date)


# %% # TODO: Needs Transform to arrays from matrices
def JD2Gregorian(JD):  # TODO: julian2gregorian_date()
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
    gregorian2julian_date : The reverse of this conversion

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
        # TODO: Check JD = 2453140.1965277777053 does not get 60 sec
        # Determine Second
        GD[5, i] = ((tau - GD[3, i] - (GD[4, i] / 60)) * 3600)
    return np.asmatrix(GD)


# %% GOOD on arrays
@gregorian_date_validation
def gregorian2julian_date(gregorian_date):
    """
    Converts Julian Date to Gregorian Date Format.
        - Accurate for dates between year 1900 - 2100

    Parameters
    ----------
    gregorian_date : array_like (n, 6)
        - Imput an numpy ndarray of n Gregorian Dates
        - Arranged in the format np.array([year, month, day, hour, min, sec])

    Returns
    -------
    julian_date : numpy matrix (n,)
        - Returns numpy ndarray of n Julian Dates as floats

    See Also
    --------
    JD2Gregorian : The reverse of this conversion

    References
    ----------
    [1] D. Vallado, `Fundamentals of Astrodynamics and Applications`. 4th ed.,
    Microcosm Press, 2013.
        - Alg. 14, pg. 183
    """
    # Determine number of date entries in gregorian_date
    gd = np.array(gregorian_date)
    try:
        gd.shape[1]
    except IndexError:
        entries = 1
    else:
        entries = gd.shape[0]
    # Convert single date or array of dates to julian dates
    if entries == 1:
        gj_1 = (367 * gd[0])
        gj_2 = np.int((7 * (gd[0] + np.int((gd[1] + 9) / 12))) / 4)
        gj_3 = np.int((275 * gd[1]) / 9)
        gj_4 = (gd[2])
        gj_5 = (((((gd[5] / 60) + gd[4]) / 60) + gd[3]) / 24)
        julian_date = np.asscalar(np.float64((gj_1) - (gj_2) + (gj_3) +
                                             (gj_4) + (1721013.5) + (gj_5)))
    else:
        julian_date = np.zeros(entries, dtype=float)
        for j in range(entries):
            gj_1 = (367 * gd[j, 0])
            gj_2 = np.int((7 * (gd[j, 0] + np.int((gd[j, 1] + 9) / 12))) / 4)
            gj_3 = np.int((275 * gd[j, 1]) / 9)
            gj_4 = (gd[j, 2])
            gj_5 = (((((gd[j, 5] / 60) + gd[j, 4]) / 60) + gd[j, 3]) / 24)
            julian_date[j] = np.asscalar(np.float64((gj_1) - (gj_2) +
                                                    (gj_3) + (gj_4) +
                                                    (1721013.5) + (gj_5)))
    return julian_date


# %% GOOD on arrays
@gregorian_date_validation
def convert_time(gregorian_date):
    """
    Converts UTC Gregorian Date to UT1, TAI, and TT timeframes.
    Outputs are all Julian Dates
        - Input timeframe must be UTC +0:00

    Parameters
    ----------
    gregorian_date : array_like (n, 6)
        - Imput an numpy ndarray of n Gregorian Dates
        - Arranged in the format np.array([year, month, day, hour, min, sec])

    Returns
    -------
    julian_date_conv : tuple (n)(4,) each entry is an ndarray of 4 julian dates
        - julian_date_conv[n][0] = julian_date_UTC
            - Julian Date in Universal Coordinated Time (UTC)
        - julian_date_conv[n][1] = julian_date_UT1
            - Julian Date in Universal Time (UT1)
        - julian_date_conv[n][2] = julian_date_TAI
            - Julian Date in International Atomic Time (TAI)
        - julian_date_conv[n][3] = julian_date_TT
            - Julian Date in Terrestrial Time (TT)

    References
    ----------
    [1] D. Vallado, `Fundamentals of Astrodynamics and Applications`. 4th ed.,
    Microcosm Press, 2013.
        - Alg. 16, pg. 195
    [2] https://datacenter.iers.org
    """
    # Initialize Vectors
    # Determine how many date entries are imputted
    try:
        gregorian_date.shape[1]
    except IndexError:
        entries = 1
    else:
        entries = gregorian_date.shape[0]

    # Get IERS data from online if local file doesn't exist
    if os.path.exists(r'sample\orbital_analyses\EOPCO4.npy'):
        EOPCO4 = np.load(r'sample\orbital_analyses\EOPCO4.npy')
    elif os.path.exists(r'sample\EOPCO4.npy'):
        EOPCO4 = np.load(r'sample\EOPCO4.npy')
    else:
        # EOP 14 C04 (IAU2000A) [May need occasional updating to propper link]
        EOP_scrape = "https://datacenter.iers.org/data/latestVersion/224_EOP_C04_14.62-NOW.IAU2000A224.txt"  # TODO: Check if correct data table
        hdr = {'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8'}
        EOP_request = Request(url=EOP_scrape, headers=hdr)
        EOP_page = urlopen(EOP_request)
        EOP_soup = BeautifulSoup(EOP_page, "lxml")
        EOP_string = str(EOP_soup.body.p.string)
        EOP_list = re.split(r'\n+', EOP_string.rstrip('\n'))
        del EOP_list[:10]  # Delete initial description lines
        EOPCO4 = np.matrix(np.zeros((np.size(EOP_list), 16), dtype=float))
        for n in range(np.size(EOP_list)):  # convert to numpy matrix
            EOPCO4[n, :] = np.fromstring(EOP_list[n], dtype=float, sep=" ")
        np.save("EOPCO4.npy", EOPCO4)
    # Get IERS data from local file
    if os.path.exists(r'sample\orbital_analyses\deltaTA.npy'):
        deltaTA = np.load(r'sample\orbital_analyses\deltaTA.npy')
    elif os.path.exists(r'deltaTA.npy'):
        deltaTA = np.load(r'deltaTA.npy')

    # Initalize array for all ndarrays to be placed into for output
    julian_date_conv = []
    for i in range(entries):
        # Handling if only one gregorian date is supplied
        if entries == 1:
            julian_date_UTC = (gregorian2julian_date(gregorian_date))
            modjulian_date_UTC = (julian_date_UTC - 2400000.5)
        else:
            julian_date_UTC = (gregorian2julian_date(gregorian_date[i, :]))
            modjulian_date_UTC = (julian_date_UTC - 2400000.5)
        # Scrape IERS Data files to get dUT1 & dAT
        EOP_index = np.searchsorted(np.ravel(EOPCO4[:, 3]), modjulian_date_UTC)
        dTA_index = np.searchsorted(np.ravel(deltaTA[:, 3]), julian_date_UTC)
        # Check if date is >= the index value and edit if needed
        if EOP_index == np.size(EOPCO4[:, 3]):
            EOP_index = (EOP_index - 1)
        elif modjulian_date_UTC != EOPCO4[EOP_index, 3]:
            EOP_index = (EOP_index - 1)
        if dTA_index == np.size(deltaTA[:, 3]):
            dTA_index = (dTA_index - 1)
        elif julian_date_UTC != deltaTA[dTA_index, 3]:
            dTA_index = (dTA_index - 1)
        # TODO: Date outside data range takes last value in data, change?
        dUT1 = np.asscalar(EOPCO4[EOP_index, 6])  # Seconds
        dAT = np.asscalar(deltaTA[dTA_index, 4])  # Seconds
        dTT = (32.184)  # Seconds
        deltas = np.array([dUT1, dAT, (dAT + dTT)])
        print(deltas)

        # Handling outputs if only one gregorian date is supplied
        if entries == 1:
            deltad = increment_gregorian_date(gregorian_date, deltas)
            julian_d_conv = np.array([julian_date_UTC, 0, 0, 0])
            for j in range(deltas.size):
                julian_d_conv[j + 1] = (gregorian2julian_date(deltad[0][j]))
        else:
            deltad = increment_gregorian_date(gregorian_date[i, :], deltas)
            julian_d_conv = np.array([julian_date_UTC, 0, 0, 0])
            for j in range(deltas.size):
                julian_d_conv[j + 1] = (gregorian2julian_date(deltad[0][j]))
        julian_date_conv.append(julian_d_conv)
    return tuple(julian_date_conv)



# %% # TODO: Needs Transform to arrays from matrices
# TODO: add handing for acceleration transforms
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

# %% # TODO: Needs Transform to arrays from matrices
# TODO: add handing for acceleration transforms
def Perifocal2Keplerian(rad, vel):
    return

# %% # TODO: Needs Work
# TODO: add handing for acceleration transforms
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

# %% # TODO: Needs Transform to arrays from matrices
# TODO: add handing for acceleration transforms
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

# %% # TODO: Needs Transform to arrays from matrices
# TODO: add handing for acceleration transforms
# TODO: Make this exist please
def Cartesian2Perifocal(rad, vel):
    return

# %% # TODO: Needs Transform to arrays from matrices
# TODO: add handing for acceleration transforms
def Perifocal2Cartesian(rad, vel):
    return
