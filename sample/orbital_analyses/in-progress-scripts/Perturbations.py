# Modules Utilized
import numpy as np
from sympy import Symbol, S
import pickle

from poliastro.twobody.perturbations import atmospheric_drag
from poliastro.twobody.perturbations import J2_perturbation

from orbital_analyses.Transform_Coordinate import ECFixed2Geodetic
from orbital_analyses.Transform_Coordinate import J20002FK5_ECFixed
from orbital_analyses.Transform_Coordinate import Geodetic2Geocentric


# DONT KNOW WHAT REFERENCE FRAME THIS IS IN
def J2_Perturbation(t0, state, k, J2, R):
    # Perturbations from the oblateness of Earth in km/sec^2
    # Howard Curtis, (12.30)
    j2_rad_v = state[:3]
    j2_rad = np.linalg.norm(j2_rad_v)

    j2_fac = (((3.0 / 2.0) * k * J2 * (R ** 2)) / (j2_rad ** 5))
    j2_a_x = (((5.0 * (j2_rad_v[2] ** 2)) / (j2_rad ** 2)) - 1)
    j2_a_y = (((5.0 * (j2_rad_v[2] ** 2)) / (j2_rad ** 2)) - 1)
    j2_a_z = (((5.0 * (j2_rad_v[2] ** 2)) / (j2_rad ** 2)) - 3)
    j2_accel = np.concatenate((j2_a_x, j2_a_y, j2_a_z), axis=0)

    # Output in array not matrix
    j2_pert = ((np.array([(j2_accel[0, 0] * j2_rad_v[0, 0]),
                          (j2_accel[1, 0] * j2_rad_v[1, 0]),
                          (j2_accel[2, 0] * j2_rad_v[2, 0])]))) * (j2_fac)

    return(j2_pert)


# def Atmosphere_Drag(t0, state, k, R, C_D, A, m, H0, rho0):
#    # Perturbations from earth's atmosphere in km/sec^2
#
#    H = norm(state[:3])
#
#    v_vec = state[3:]
#    v = norm(v_vec)
#    B = C_D * A / m
#    rho = rho0 * np.exp(-(H - R) / H0)
#
#    return -(1.0 / 2.0) * rho * B * v * v_vec


def Earth_Gravity(eg_rad, eg_vel, eg_GD):
    # Projections of acceleration from earth gravitation force in ECEF frame
    # Degree - n, Order - m

    # Need inputs of radius, latuitude, and longitude | Store in eg_LLR
    # Other constants for this sim below:
    eg_GM = 398600.4418  # km^3/sec^2
    eg_R = 6378.1370  # km

    # Convert Inertial to ECEF
    eg_ECEF, eg_vel = J20002FK5_ECFixed(eg_rad, eg_vel, eg_GD)

    # Convert ECEF to Rad Lat Lon
    eg_llr = ECFixed2Geodetic(eg_ECEF)

    # Normalized Associated Legendre Polynomial
    # n = var('n')  # Degree
    # m = var('m')  # Order
    # Legrende Function

###############################################################################
    # EX
#    eg_LLR = np.matrix([[-81.5158], [100], [27.6648]], dtype=np.float64)

    # Inititlize variables
    eg_n = np.linspace(2.0, 200, num=199, dtype=np.int)
    u = Symbol('u')

    # Initialize output matrices
    eg_fun = np.asmatrix(np.zeros((201, 1), dtype=np.float))
    eg_dif = np.asmatrix(np.zeros((201, 1), dtype=np.float))

    # Calculate eg_u for each step and substitute into general formulas
    with open('orbital_analyses/General_FNALFs/DO_200', 'rb') as fp:
        DO_200 = pickle.load(fp)
    with open('orbital_analyses/General_FNALFs/DO_200_dif', 'rb') as fp:
        DO_200_dif = pickle.load(fp)
    FNALF_200 = DO_200[-201:]
    FNALF_200_dif = DO_200_dif[-201:]

    EGM2008_200 = np.load('orbital_analyses/EGM2008_200.npy')
    C_nm = np.transpose(np.asmatrix(EGM2008_200[:, 2]))
    S_nm = np.transpose(np.asmatrix(EGM2008_200[:, 3]))
    C_nm = C_nm[-201:, 0].copy()
    S_nm = S_nm[-201:, 0].copy()

    for i in range(201):
        # Calculate cos value for step's geocentric lat / convert to geocentric
        eg_LLR = Geodetic2Geocentric(eg_llr)
        eg_u = np.cos(eg_LLR[0, 0] * (np.pi / 180))

        # Replace u with calculated value above after Sympyfing
        eg_fun[i, 0] = (S(FNALF_200[i])).subs(u, eg_u)
        eg_dif[i, 0] = (S(FNALF_200_dif[i])).subs(u, eg_u)

###############################################################################
    # https://www.hindawi.com/journals/ijap/2014/903026/
    # Legrende Polynomial (Defined by n)

    eg_n = np.linspace(2.0, 200, num=199, dtype=np.int)
    eg_m = np.linspace(0.0, 200, num=201, dtype=np.int)

    # Calculate trig values for step's geocentric latitude
    eg_c_lat = np.cos(eg_LLR[0, 0] * (np.pi / 180))
    eg_s_lat = np.sin(eg_LLR[0, 0] * (np.pi / 180))

    # Initialize vectors
    drad_msum = np.asmatrix(np.zeros((201, 1), dtype=np.float64))
    dlat_msum = np.asmatrix(np.zeros((201, 1), dtype=np.float64))
    dlon_msum = np.asmatrix(np.zeros((201, 1), dtype=np.float64))
    drad_nsum = np.asmatrix(np.zeros((199, 1), dtype=np.float64))
    dlat_lon_nsum = np.asmatrix(np.zeros((199, 1), dtype=np.float64))

    for n in eg_n:
        drad_nsum[n - 2, 0] = ((n + 1) *
                               ((eg_R / (eg_R + eg_LLR[2, 0])) ** (n + 2)))
        dlat_lon_nsum[n - 2, 0] = ((eg_R / (eg_R + eg_LLR[2, 0])) ** n)

    # Unitless
    drad_n = np.sum(drad_nsum)
    dlat_lon_n = np.sum(dlat_lon_nsum)

    for m in eg_m:
        eg_cm = np.cos(eg_LLR[1, 0] * (np.pi / 180) * m)
        eg_sm = np.sin(eg_LLR[1, 0] * (np.pi / 180) * m)
        # satrts 20097 at n=200, m=0 -> 200

        drad_msum[m, 0] = (eg_fun[m, 0] * eg_s_lat *
                           ((C_nm[m, 0] * eg_cm) + (S_nm[m, 0] * eg_sm)))
        dlat_msum[m, 0] = (eg_dif[m, 0] * eg_c_lat *
                           ((C_nm[m, 0] * eg_cm) + (S_nm[m, 0] * eg_sm)))
        dlon_msum[m, 0] = (m * eg_fun[m, 0] * eg_s_lat *
                           ((-C_nm[m, 0] * eg_cm) + (S_nm[m, 0] * eg_sm)))

    drad_m = np.sum(drad_msum)
    dlat_m = np.sum(dlat_msum)
    dlon_m = np.sum(dlon_msum)

    # Combine all sums into final differential equations
    eg_dU_drad = ((-eg_GM / eg_R) * (1 + (drad_n * drad_m)))
    eg_dU_dlat = ((eg_GM / eg_R) * (dlat_lon_n * dlat_m))
    eg_dU_dlon = ((eg_GM / eg_R) * (dlat_lon_n * dlon_m))

    # Calculate remaining components
    eg_r = np.sqrt((eg_ECEF[0, 0] ** 2) + (eg_ECEF[1, 0] ** 2) +
                   (eg_ECEF[2, 0] ** 2))
    eg_rho = np.sqrt((eg_ECEF[0, 0] ** 2) + (eg_ECEF[1, 0] ** 2))
    eg_x = eg_ECEF[0, 0]
    eg_y = eg_ECEF[1, 0]
    eg_z = eg_ECEF[2, 0]

    # Accelerations in ECEF
    eg_Ax = ((eg_dU_drad * (eg_x / eg_r)) -
             (eg_dU_dlat * ((eg_x * eg_z) / (eg_rho * (eg_r ** 2)))) -
             (eg_dU_dlon * (eg_y / (eg_rho ** 2))))

    eg_Ay = ((eg_dU_drad * (eg_y / eg_r)) -
             (eg_dU_dlat * ((eg_y * eg_z) / (eg_rho * (eg_r ** 2)))) +
             (eg_dU_dlon * (eg_x / (eg_rho ** 2))))

    eg_Az = ((eg_dU_drad * (eg_z / eg_r)) +
             (eg_dU_dlon * (eg_rho / (eg_r ** 2))))

    return eg_Ax, eg_Ay, eg_Az

###############################################################################
# Algorithm to access location in degree / order list

c = 0
for a in eg_n:
    for b in range(a + 1):
        c = c + 1
print(c)
