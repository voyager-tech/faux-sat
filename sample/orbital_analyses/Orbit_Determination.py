# Utilized Modules
import numpy as np
import struct
from copy import deepcopy as dc
from orbital_analyses import u
from skopt import dummy_minimize
from joblib import Parallel
import time

import Requirements as Req
from parts_list import Parts_List
from orbital_analyses.Transform_State import TimeAdjust
from orbital_analyses.Propagation import Prop_Cowell
from orbital_analyses.Transform_State import Keplerian2Cartesian
from orbital_analyses.Transform_Coordinate import Geodetic2ECFixed
from orbital_analyses.Transform_Coordinate import FK5_ECFixed2J2000

getBin = lambda x: x > 0 and str(bin(x))[2:] or str(bin(x))[3:]


def floatToBinary(value):
    val = struct.unpack('Q', struct.pack('d', value))[0]
    return getBin(val)


def binaryToFloat(value):
    hx = hex(int(value, 2))
    return struct.unpack("d", struct.pack("q", int(hx, 16)))[0]


def Prop_xyz(x):
    rad, vel, jd, gd = Prop_Cowell(Sat_Rad_inrt, Sat_Vel_inrt,
                                   Req.GD_UTC, 2, (x[0] * u.sec))
    Sat_Ri = rad[:, 1]
    Sat_gd = gd[:, 1]
    # Determine Rad/Vel of the GS in the Inertial Frame
    GS_Rad_inrt = np.asmatrix(np.zeros((3, 1), dtype=np.float64))
    GS_Vel_inrt = np.asmatrix(np.zeros((3, 1), dtype=np.float64))
    GS_fix = Geodetic2ECFixed(GS_lla[:, s])
    GS_zero = np.matrix([[0.0], [0.0], [0.0]])
    GS_Rad_inrt[:, 0], GS_Vel_inrt[:, 0] = FK5_ECFixed2J2000(GS_fix, GS_zero,
                                                             Sat_gd)
    # Determine difference in Rad - variable to minimize
    Diff_lon = np.sqrt(((Sat_Ri[0, 0] - GS_Rad_inrt[0, 0]) ** 2) +
                       ((Sat_Ri[1, 0] - GS_Rad_inrt[1, 0]) ** 2) +
                       ((Sat_Ri[2, 0] - GS_Rad_inrt[2, 0]) ** 2))
    return (Diff_lon)


# Pg. 891 Vallado
# =============================================================================
# Limitations of Method:
# Can only be based on 1 ground location
# Only passes the ground location at the periapsis
# Limitations based on view range of sensor
#
# Guarentees of Method:
# Orbit will initially pass over the ground loaction
#
#
# What can we Vary:
# Length of the repeat period (days)
# Inclination (2 choices) (ex: 43 deg or -43 deg)
#
# What Variations cause:
# Change in Eccentricity
# Change in SMA (Apoapsis height)
# Changes in RAAN / AoP (if INC is changed)
#
# =============================================================================

# %%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Determine Initial Version of Orbit
GS_lat = Req.ex_gs_ll[0, 0]
GS_lon = Req.ex_gs_ll[1, 0]
GS_alt = Req.ex_gs_a
cam_alt = Parts_List['Payload']['operational_alt']
cam_fov = (20 * u.deg)  # .plus_minus(15)
cam_fov_rang = 10 * u.deg

# Initialize Constants
R_Earth = 6378.1363 * u.km
Mu_Earth = 398600.4418 * ((u.km ** 3) / (u.sec ** 2))
e2 = 0.006694385
w0_Earth = 0.0000729211585530 * (u.rad / u.sec)
J2_Earth = 0.0010826267
f_Earth = (1 / 298.257)

# Determine ellipsoidal Earth's coordinates
C_Earth = (R_Earth / (np.sqrt(1 - (e2 * (np.sin(GS_lat) ** 2)))))
S_Earth = C_Earth * (1 - e2)
R_delt = ((C_Earth + GS_alt) * np.cos(GS_lat))
R_k = ((S_Earth + GS_alt) * np.sin(GS_lat))
R_mag = np.sqrt((R_delt ** 2) + (R_k ** 2))

# Determine radius of the Perigee
R_per = (R_mag + cam_alt)

# Max amount of drift before location is out of sight
gamma = np.arcsin(((cam_alt + R_Earth) * np.sin(cam_fov)) / R_Earth).to(u.deg)
# TODO: Confirm this is the only check needed
if gamma.magnitude <= 90:
    gamma = (180 * u.deg) - gamma

rho = (R_Earth * np.cos(gamma) + (cam_alt + R_Earth) * np.cos(cam_fov))

GR_max = np.arcsin((rho * np.sin(cam_fov)) / R_Earth)
GR_max_km = (R_Earth * GR_max.magnitude)

# Determine equator crossing value
# TODO: JANK - Determine best day repeat period
GR3 = np.ceil(GR_max_km * 3.25)
K_rev2rep = np.floor((2 * np.pi * R_Earth) / GR3) * u.turn
K_day2rep = Req.GS_repeat

# First Iteration of Orbit Variables
n0 = (K_rev2rep / K_day2rep).to(u.rad / u.sec)


# Determine SMA
a0 = ((Mu_Earth * ((1 / n0.to(u.rad / u.sec)) ** 2)) ** (1 / 3)).to(u.km)
# Determine ECC to keep same value of R_per determined earlier
e0 = ((a0 - R_per) / a0)
# Set INC value as initial guess
i0 = GS_lat - (15 * u.deg)

RAAN_i = (-((3 * n0 * (R_Earth ** 2) * J2_Earth) / (2 * (a0 ** 2))) *
          np.cos(i0))
AoP_i = (((3 * n0 * (R_Earth ** 2) * J2_Earth) / (4 * (a0 ** 2))) *
         (4 - (5 * (np.cos(i0) ** 2))))

# Initial Guess
init = ((n0 / w0_Earth).magnitude)
Iterate = 3
frac = init
PQ_loop = np.asmatrix(np.zeros((Iterate, 2), dtype=np.int))

for n in range(Iterate):
    alphan = np.int(frac)
    bn = (1 / np.mod(frac, 1))
    if n == 0:
        PQ_loop[n, 0] = alphan
        PQ_loop[n, 1] = 1
        frac = dc(bn)
    if n == 1:
        PQ_loop[n, 0] = (1 + (PQ_loop[0, 0] * alphan))
        PQ_loop[n, 1] = alphan
        frac = dc(bn)
    if n > 1:
        PQ_loop[n, 0] = ((alphan * PQ_loop[n - 1, 0]) + PQ_loop[n - 2, 0])
        PQ_loop[n, 1] = ((alphan * PQ_loop[n - 1, 1]) + PQ_loop[n - 2, 1])
        frac = dc(bn)

Orbit_val = np.asmatrix(np.zeros((7, Iterate + 1), dtype=np.float64))
Orbit_val[0, 0] = n0.magnitude
Orbit_val[1, 0] = a0.magnitude
Orbit_val[2, 0] = e0.magnitude
Orbit_val[3, 0] = i0.magnitude

for j in range(Iterate):
    # This loop for iterating through PQ pairs
    P_i = PQ_loop[j, 0]
    Q_i = PQ_loop[j, 1]
    # Ititalize loop variables
    n_i = n0
    sma_i = a0
    ecc_i = e0
    inc_i = i0
    delta_sma = 1.0  # Initial value
    k = 0
    while k <= 5:
        # This loop for iterating through orbits of a PQ pair
        # Calculate intermediate values (pg.649-652) Vallado
        RAAN_i = (-((3 * n_i * (R_Earth ** 2) * J2_Earth) /
                    (2 * (sma_i ** 2))) * np.cos(inc_i))
        AoP_i = (((3 * n_i * (R_Earth ** 2) * J2_Earth) /
                  (4 * (sma_i ** 2))) * (4 - (5 * (np.cos(inc_i) ** 2))))
        # Calculate updated Orbital Variables
        n_j = ((P_i / Q_i) * (w0_Earth - RAAN_i) *
               (1 - ((3 * J2_Earth) / 2) *
                ((R_Earth / a0) ** 2) * (3 - (4 * (np.sin(i0) ** 2)))))
        sma_j = ((Mu_Earth * ((1 / n_j) ** 2)) ** (1 / 3)).to(u.km)
        ecc_j = ((sma_j - R_per) / sma_j)
        # Update while loop condition
        delta_sma = (np.abs(n_j - n_i)).magnitude
        # Update Loop Variables
        n_i = n_j
        sma_i = sma_j
        ecc_i = ecc_j
        k = k + 1
        # Check if SMA is smaller than R_per (Incompatibe Orbit)
        if sma_i.magnitude < R_per.magnitude:
            delta_sma = 0
            n_i = 0
            sma_i = 0 * u.dimensionless
            ecc_i = 0 * u.dimensionless
            inc_i = 0 * u.dimensionless
            break
    Orbit_val[0, j + 1] = n_i.magnitude
    Orbit_val[1, j + 1] = sma_i.magnitude
    Orbit_val[2, j + 1] = ecc_i.magnitude
    Orbit_val[3, j + 1] = inc_i.magnitude

for l in range(Iterate):
    if Orbit_val[1, l + 1] == 0:
        Orbit_val[4, l + 1] = 0.0
        Orbit_val[5, l + 1] = 0.0
        Orbit_val[6, l + 1] = 0.0
    else:
        # Intermediate values for OMEGA_0
        beta = np.arcsin((np.cos(Orbit_val[3, l + 1] * u.deg)) /
                         (np.cos(GS_lat))).to(u.deg)
        Lam_u = np.arccos((np.cos(beta)) /
                          (np.sin(Orbit_val[3, l + 1] * u.deg))).to(u.deg)
        # Determine Greenwich Sidereal Time - ALG 15 - Vallado
        jd_UTC, jd_UT1, jd_TAI, jd_TT = TimeAdjust(Req.GD_UTC)
        T_UT1 = ((jd_UT1 - 2451545.0) / 36525)
        Theta_s = (67310.54841 +
                   (((876600 * 3600) + 8640184.812866) * T_UT1) +
                   (0.093104 * (T_UT1 ** 2)) - (6.2e-06 * (T_UT1 ** 3)))  # Sec
        # Reduce to within 86400 seconds
        reduc = (Theta_s / 86400)
        diff = reduc - np.floor(reduc)
        Theta_sec = (86400 * diff)  # Seconds
        # Convert from Seconds to Degrees
        Theta_d = (Theta_sec * (1 / 240))  # Degrees
        if Theta_d < 0:
            GMST_d = (Theta_d + 360) * u.deg
        if Theta_d > 0:
            GMST_d = Theta_d * u.deg
            # Determine Remaining Orbit Values
            RAAN_0 = (GMST_d - Lam_u + GS_lon)
            AoP_0 = np.arcsin((np.sin(GS_lat)) /
                              (np.sin(Orbit_val[3, l + 1] *
                                      u.deg))).to(u.deg)
            M_0 = (0.0 * u.deg)
            # Input into final matrix as l+1
            Orbit_val[4, l + 1] = (RAAN_0.to(u.deg)).magnitude
            if Orbit_val[4, l + 1] < 0:
                Orbit_val[4, l + 1] = Orbit_val[4, l + 1] + 360.
            Orbit_val[5, l + 1] = (AoP_0.to(u.deg)).magnitude
            if Orbit_val[5, l + 1] < 0:
                Orbit_val[5, l + 1] = Orbit_val[5, l + 1] + 360.
            # TODO: swap out the added value
            Orbit_val[6, l + 1] = (M_0.magnitude + 180)
            if Orbit_val[6, l + 1] < 0:
                Orbit_val[6, l + 1] = Orbit_val[6, l + 1] + 360.

# %%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Genetic Algorithm for determining best orbit
# https://www.researchgate.net/profile/D_Mortari/publication/254199757_Reconnaissance_Problem_Using_Genetic_Algorithms/links/5538e83d0cf247b8587d7edf/Reconnaissance-Problem-Using-Genetic-Algorithms.pdf

# TODO: Remove and fix above here
# Initial Orbit Guess
Orbit_val = np.asmatrix(np.zeros((7, 3), dtype=np.float64))
Orbit_val[0, 2] = 0.1
Orbit_val[1, 2] = 10000.
Orbit_val[2, 2] = 0.00001
Orbit_val[3, 2] = 0.
Orbit_val[4, 2] = 0.
Orbit_val[5, 2] = 0.
Orbit_val[6, 2] = 0.

# Constants
E_Rad = 6378.137 * u.km                              # Earth Radius
E_w0 = 0.0000729211585530 * (u.rad / u.sec)          # Earth angular velocity
E_Mu = 398600.4418 * ((u.km ** 3) / (u.sec ** 2))    # Earth Mu
steps = Req.Steps
step_size = Req.Ssize

# Number of Ground Points to consider
Stations = 4
GS_lla = np.asmatrix(np.zeros((3, Stations), dtype=np.float64))
GS_lla[0, 0] = (Req.gs_ll[0, 0].to(u.deg)).magnitude                # Deg
GS_lla[1, 0] = (Req.gs_ll[1, 0].to(u.deg)).magnitude                # Deg
GS_lla[2, 0] = (Req.gs_a.to(u.km)).magnitude                        # km
GS_lla[0, 1] = 40
GS_lla[1, 1] = 238
GS_lla[2, 1] = 0.600
GS_lla[0, 2] = 34
GS_lla[1, 2] = 244
GS_lla[2, 2] = 1.5
GS_lla[0, 3] = -5
GS_lla[1, 3] = 300
GS_lla[2, 3] = .4

# Determine initial Rad value for Satellite orbit
Kep = dc(Orbit_val[1:, 2])
Sat_Rad_inrt, Sat_Vel_inrt = Keplerian2Cartesian(Kep)

# %%###########################################################################
Population = 7
Member = 8

# List all 8 population members (1 is initial orbit from above, 7 random)
# TODO: Make all of these limits more dynamic / set via requirement limitations
Pop_sma = np.asmatrix(np.zeros((Population + 1, Member), dtype=np.float64))
Pop_sma[0, 0] = Kep[0, 0]
Pop_sma[0, 1:] = np.random.uniform(E_Rad.magnitude + 300.,
                                   E_Rad.magnitude + 2000., 7)  # km

Pop_ecc = np.asmatrix(np.zeros((Population + 1, Member), dtype=np.float64))
Pop_ecc[0, 0] = Kep[1, 0]
Pop_ecc[0, 1:] = np.random.uniform(0., 0.2, 7)

Pop_inc = np.asmatrix(np.zeros((Population + 1, Member), dtype=np.float64))
Pop_inc[0, 0] = Kep[2, 0]
Pop_inc[0, 1:] = np.random.uniform(0., 180., 7)  # deg

Pop_raan = np.asmatrix(np.zeros((Population + 1, Member), dtype=np.float64))
Pop_raan[0, 0] = Kep[3, 0]
Pop_raan[0, 1:] = np.random.uniform(0., 360., 7)  # deg

Pop_aop = np.asmatrix(np.zeros((Population + 1, Member), dtype=np.float64))
Pop_aop[0, 0] = Kep[4, 0]
Pop_aop[0, 1:] = np.random.uniform(0., 360., 7)  # deg

Pop_ta = np.asmatrix(np.zeros((Population + 1, Member), dtype=np.float64))
Pop_ta[0, 0] = Kep[5, 0]
Pop_ta[0, 1:] = np.random.uniform(0., 360., 7)  # deg

# Determine initial Rad value for Satellite orbit
Fit_Terr = np.asmatrix(np.zeros((Population, Member), dtype=np.float64))
for p in range(Population):
    print("p=%s" % p)  # TODO: Get Rid of Eventually
    # Evaluate all Random orbits to see if they are viable
    # TODO: Make less random fixes here
    Flag = 0
    while Flag == 0:
        Sflag = 0
        for l in range(Member):
            # Check if orbit is viable
            Peri = (Pop_sma[p, l] * (1 - Pop_ecc[p, l]))  # km
            if Peri <= (E_Rad.magnitude + 160):  # If orbit below decay range
                Pop_sma[p, l] = np.random.uniform(E_Rad.magnitude + 300.,
                                                  E_Rad.magnitude + 2000.)
                Pop_ecc[p, l] = np.random.uniform(0., 0.2)
                Sflag = Sflag + 1
        if Sflag == 0:
            Flag = 1
#    start = time.time()  # TODO: Get Rid of Eventually
    # Run Fitness simulation
    for m in range(Member):
        # Determine initial Rad value for Satellite orbit
        Keplerian = np.matrix([[Pop_sma[p, m]], [Pop_ecc[p, m]],
                               [Pop_inc[p, m]], [Pop_raan[p, m]],
                               [Pop_aop[p, m]], [Pop_ta[p, m]]])
        Sat_Rad_inrt, Sat_Vel_inrt = Keplerian2Cartesian(Keplerian)
        # Evaluate Orbit with Fit Function
        Fit_err = np.asmatrix(np.zeros((1, Stations), dtype=np.float64))
        Fit_Step = np.asmatrix(np.zeros((1, Stations), dtype=np.float64))
        # Determine Rad of satellite to compare location with ground stations
        for s in range(Stations):
            # TODO: Figure out how to do parallel processing
            res = dummy_minimize(Prop_xyz, [(1.0, 86400.0)], n_calls=50)
            Fit_err[:, s] = res['fun']
#            Fit_Step[:, s] = float(res['x'][0])  # TODO: Get Rid of Eventually
        Fit_Terr[p, m] = np.sum(Fit_err[0, :])
#    end = time.time()  # TODO: Get Rid of Eventually
#    print(end - start)  # TODO: Get Rid of Eventually
    # Using fitness values, select and create next generation
    h = 4  # Half of Population
    Pop_fit = np.array([Fit_Terr[p, 0], Fit_Terr[p, 1], Fit_Terr[p, 2],
                        Fit_Terr[p, 3], Fit_Terr[p, 4], Fit_Terr[p, 5],
                        Fit_Terr[p, 6], Fit_Terr[p, 7]])
    idx = np.argpartition(Pop_fit, h)  # 0:3 hold indices of 4 lowest values
    Pop_new = np.matrix([Fit_Terr[p, idx[0]], Fit_Terr[p, idx[1]],
                         Fit_Terr[p, idx[2]], Fit_Terr[p, idx[3]], 0, 0, 0, 0])
    # Convert orbital Paramaters to Binary Strings for the Genetic Algorithm
    # TODO: assumed that all converted keplerian values are always positive
    B_sma = []
    for k in range(h):
        B_sma.append(floatToBinary(Pop_sma[p, idx[k]]))

    B_ecc = []
    for k in range(h):
        B_ecc.append(floatToBinary(Pop_ecc[p, idx[k]]))

    B_inc = []
    for k in range(h):
        B_inc.append(floatToBinary(Pop_inc[p, idx[k]]))

    B_raan = []
    for k in range(h):
        B_raan.append(floatToBinary(Pop_raan[p, idx[k]]))

    B_aop = []
    for k in range(h):
        B_aop.append(floatToBinary(Pop_aop[p, idx[k]]))

    B_ta = []
    for k in range(h):
        B_ta.append(floatToBinary(Pop_ta[p, idx[k]]))
    # Create 2nd half of population through Recombination (IEEE 754 - 64 Bit)
    # TODO: Could use some actual reasoning here in the randomization....
    # TODO: Consider adjusting alteration ranges for later populations
    # SMA
    pair = np.random.choice(range(4), 4, replace=False)
    ran = np.random.choice(range(7), 4)
    B_sma.append(B_sma[pair[0]][:(14 + ran[0])] +
                 B_sma[pair[1]][(14 + ran[0]):])
    B_sma.append(B_sma[pair[1]][:(14 + ran[1])] +
                 B_sma[pair[0]][(14 + ran[1]):])
    B_sma.append(B_sma[pair[2]][:(14 + ran[2])] +
                 B_sma[pair[3]][(14 + ran[2]):])
    B_sma.append(B_sma[pair[3]][:(14 + ran[3])] +
                 B_sma[pair[2]][(14 + ran[3]):])
    # ECC
    pair = np.random.choice(range(4), 4, replace=False)
    ran = np.random.choice(range(7), 4)
    B_ecc.append(B_ecc[pair[0]][:(14 + ran[0])] +
                 B_ecc[pair[1]][(14 + ran[0]):])
    B_ecc.append(B_ecc[pair[1]][:(14 + ran[1])] +
                 B_ecc[pair[0]][(14 + ran[1]):])
    B_ecc.append(B_ecc[pair[2]][:(14 + ran[2])] +
                 B_ecc[pair[3]][(14 + ran[2]):])
    B_ecc.append(B_ecc[pair[3]][:(14 + ran[3])] +
                 B_ecc[pair[2]][(14 + ran[3]):])
    # INC
    pair = np.random.choice(range(4), 4, replace=False)
    ran = np.random.choice(range(7), 4)
    B_inc.append(B_inc[pair[0]][:(14 + ran[0])] +
                 B_inc[pair[1]][(14 + ran[0]):])
    B_inc.append(B_inc[pair[1]][:(14 + ran[1])] +
                 B_inc[pair[0]][(14 + ran[1]):])
    B_inc.append(B_inc[pair[2]][:(14 + ran[2])] +
                 B_inc[pair[3]][(14 + ran[2]):])
    B_inc.append(B_inc[pair[3]][:(14 + ran[3])] +
                 B_inc[pair[2]][(14 + ran[3]):])
    # RAAN
    pair = np.random.choice(range(4), 4, replace=False)
    ran = np.random.choice(range(7), 4)
    B_raan.append(B_raan[pair[0]][:(14 + ran[0])] +
                  B_raan[pair[1]][(14 + ran[0]):])
    B_raan.append(B_raan[pair[1]][:(14 + ran[1])] +
                  B_raan[pair[0]][(14 + ran[1]):])
    B_raan.append(B_raan[pair[2]][:(14 + ran[2])] +
                  B_raan[pair[3]][(14 + ran[2]):])
    B_raan.append(B_raan[pair[3]][:(14 + ran[3])] +
                  B_raan[pair[2]][(14 + ran[3]):])
    # AOP
    pair = np.random.choice(range(4), 4, replace=False)
    ran = np.random.choice(range(7), 4)
    B_aop.append(B_aop[pair[0]][:(14 + ran[0])] +
                 B_aop[pair[1]][(14 + ran[0]):])
    B_aop.append(B_aop[pair[1]][:(14 + ran[1])] +
                 B_aop[pair[0]][(14 + ran[1]):])
    B_aop.append(B_aop[pair[2]][:(14 + ran[2])] +
                 B_aop[pair[3]][(14 + ran[2]):])
    B_aop.append(B_aop[pair[3]][:(14 + ran[3])] +
                 B_aop[pair[2]][(14 + ran[3]):])
    # TA
    pair = np.random.choice(range(4), 4, replace=False)
    ran = np.random.choice(range(7), 4)
    B_ta.append(B_ta[pair[0]][:(14 + ran[0])] + B_ta[pair[1]][(14 + ran[0]):])
    B_ta.append(B_ta[pair[1]][:(14 + ran[1])] + B_ta[pair[0]][(14 + ran[1]):])
    B_ta.append(B_ta[pair[2]][:(14 + ran[2])] + B_ta[pair[3]][(14 + ran[2]):])
    B_ta.append(B_ta[pair[3]][:(14 + ran[3])] + B_ta[pair[2]][(14 + ran[3]):])
    # Add in 1 random mutation per element population
    mut = np.random.choice(range(50), 6)
    string = np.random.choice(range(8), 6)
    # SMA
    mid = B_sma[string[0]]
    if mid[(12 + mut[0])] == '1':
        new = '0'
    elif mid[(12 + mut[0])] == '0':
        new = '1'
    mutated = (mid[:(12 + mut[0])] + new + mid[((12 + mut[0]) + 1):])
    B_sma[string[0]] = mutated
    if len(B_sma[string[0]]) > 63:
        B_sma[string[0]] = B_sma[string[0]][:-1]
    # ECC
    mid = B_ecc[string[1]]
    if mid[(12 + mut[1])] == '1':
        new = '0'
    elif mid[(12 + mut[1])] == '0':
        new = '1'
    mutated = (mid[:(12 + mut[1])] + new + mid[((12 + mut[1]) + 1):])
    B_ecc[string[1]] = mutated
    if len(B_ecc[string[1]]) > 63:
        B_ecc[string[1]] = B_ecc[string[1]][:-1]
    # INC
    mid = B_inc[string[2]]
    if mid[(12 + mut[2])] == '1':
        new = '0'
    elif mid[(12 + mut[2])] == '0':
        new = '1'
    mutated = (mid[:(12 + mut[2])] + new + mid[((12 + mut[2]) + 1):])
    B_inc[string[2]] = mutated
    if len(B_inc[string[2]]) > 63:
        B_inc[string[2]] = B_inc[string[2]][:-1]
    # RAAN
    mid = B_raan[string[3]]
    if mid[(12 + mut[3])] == '1':
        new = '0'
    elif mid[(12 + mut[3])] == '0':
        new = '1'
    mutated = (mid[:(12 + mut[3])] + new + mid[((12 + mut[3]) + 1):])
    B_raan[string[3]] = mutated
    if len(B_raan[string[3]]) > 63:
        B_raan[string[3]] = B_raan[string[3]][:-1]
    # AOP
    mid = B_aop[string[4]]
    if mid[(12 + mut[4])] == '1':
        new = '0'
    elif mid[(12 + mut[4])] == '0':
        new = '1'
    mutated = (mid[:(12 + mut[4])] + new + mid[((12 + mut[4]) + 1):])
    B_aop[string[4]] = mutated
    if len(B_aop[string[4]]) > 63:
        B_aop[string[4]] = B_aop[string[4]][:-1]
    # TA
    mid = B_ta[string[5]]
    if mid[(12 + mut[5])] == '1':
        new = '0'
    elif mid[(12 + mut[5])] == '0':
        new = '1'
    mutated = (mid[:(12 + mut[5])] + new + mid[((12 + mut[5]) + 1):])
    B_ta[string[5]] = mutated
    if len(B_ta[string[5]]) > 63:
        B_ta[string[5]] = B_ta[string[5]][:-1]
    # Convert back to Intgers and place into population matrix for next iter.
    # TODO: Add limit chaecks when converting back and add 180 or 360 as needed
    for c in range(8):
        # SMA
        Pop_sma[p + 1, c] = binaryToFloat(B_sma[c])
        # ECC
        Pop_ecc[p + 1, c] = binaryToFloat(B_ecc[c])
        while Pop_ecc[p + 1, c] > 1:
            Pop_ecc[p + 1, c] = Pop_ecc[p + 1, c] * 0.1
        # INC
        Pop_inc[p + 1, c] = binaryToFloat(B_inc[c])
        while Pop_inc[p + 1, c] > 180:
            Pop_inc[p + 1, c] = Pop_inc[p + 1, c] - 180
        # RAAN
        Pop_raan[p + 1, c] = binaryToFloat(B_raan[c])
        while Pop_raan[p + 1, c] > 360:
            Pop_raan[p + 1, c] = Pop_raan[p + 1, c] - 360
        # AOP
        Pop_aop[p + 1, c] = binaryToFloat(B_aop[c])
        while Pop_aop[p + 1, c] > 360:
            Pop_aop[p + 1, c] = Pop_aop[p + 1, c] - 360
        # TA
        Pop_ta[p + 1, c] = binaryToFloat(B_ta[c])
        while Pop_ta[p + 1, c] > 360:
            Pop_ta[p + 1, c] = Pop_ta[p + 1, c] - 360
# Fnd lowest fit value of final population and return orbital elements
Pop_best = Fit_Terr[Population - 1, :]
Best = np.array([Pop_best[0, 0], Pop_best[0, 1], Pop_best[0, 2],
                 Pop_best[0, 3], Pop_best[0, 4], Pop_best[0, 5],
                 Pop_best[0, 6], Pop_best[0, 7]])
index = np.argpartition(Best, 0)  # 0:3 hold indices of 4 lowest values
In = 1
for n in range(8):
    if index[n] == 0:
        In = n
        break
    else:
        In = 1
Kep_Final = np.matrix([[Pop_sma[Population, In]], [Pop_ecc[Population, In]],
                       [Pop_inc[Population, In]], [Pop_raan[Population, In]],
                       [Pop_aop[Population, In]], [Pop_ta[Population, In]]])

# finsh with gradient optimization using final guess from GA




















# %%###########################################################################
###############################################################################

# Thoughts
# Don't mess with the exponent values
# Recombination will occur only within the fraction portion

# from 3 past decimal to 10

B_sma[0][:11]  # 10000001011
B_sma[0][11:]  # 1110010001001100110110000000100101110000000101111110
          ||     ||
100000010111110010001001100110110000000100101110000000101111110
100000010111101111010000010101000100111001011100111010010010000

binaryToFloat('100000010111110010010000010101000100111001011100111010010010000')

binaryToFloat('100000010111110010001001100110110000000100101110000000101111110')

# floats are represented by IEEE 754 floating-point format which are
# 64 bits long (not 32 bits)

0b0100000000101000010011101111011111101111111001100010110001101111

12.15423536
sign = 0  # + or -
Exponent = 10000000010  # 2^3 - 11 points long
Mantissa = 1000010011101111011111101111111001100010110001101111  # Decimal - 52 points long, ()

# finsh with gradient optimization using final guess from GA
