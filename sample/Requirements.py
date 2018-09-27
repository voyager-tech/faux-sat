# Utilized Modules
import numpy as np
from orbital_analyses import u
# from parts_list import Parts_List
from orbital_analyses.Transform_State import Gregorian2JD

# TODO: Add uncertainties in with pint
# import uncertainties

# Scripts contains all variables that correspond to various requirements

# TODO: Iterations
# Initial iteration should be completly with user defined variables aside from any user added components.
# Following iterations should aim to use more part data as the iterations continue.
# Continue iterating until a solution is found where all requirements are set by data from real parts. 

# %% Initial Iteration

###############################################################################
# Testing values (ITRF)
Rad = np.matrix([[-1033.4793830],
                 [7901.2952754],
                 [6380.3565958]])
Vel = np.matrix([[-3.225636520],
                 [-2.872451450],
                 [5.531924446]])
# Initial Gregorian Date ([year; month; day; hour; minute; second])
GD_UTC = np.zeros((6, 1), dtype=np.float64)
GD_UTC[0, 0] = (2004)
GD_UTC[1, 0] = (4)
GD_UTC[2, 0] = (6)
GD_UTC[3, 0] = (7)
GD_UTC[4, 0] = (51)
GD_UTC[5, 0] = (28.386009)
GD_UTC = np.asmatrix(GD_UTC)


# Mission Requirements
###############################################################################

## Initial Radius and Velocity Vectors (INERTIAL) (J2000)
#Rad = np.matrix([[-2374.333753542197],
#                 [5963.296226761164],
#                 [3105.173918117055]]) * u.km
#Vel = np.matrix([[-5.658163685962852],
#                 [-3.49897296995646],
#                 [4.064517100619342]]) * (u.km * u.s)
#
## Initial Gregorian Date ([year; month; day; hour; minute; second])
#GD_UTC = np.zeros((6, 1), dtype=np.float64)
#GD_UTC[0, 0] = (2004)
#GD_UTC[1, 0] = (4)
#GD_UTC[2, 0] = (6)
#GD_UTC[3, 0] = (0)
#GD_UTC[4, 0] = (0)
#GD_UTC[5, 0] = (0)
#GD_UTC = np.asmatrix(GD_UTC)
#
#JD_UTC = Gregorian2JD(GD_UTC)

# Propagation info
Ssize = 1.0 * u.min
Steps = 600
Stime = (Steps * Ssize)
# Initial step size
# Max and min step size, max step attempts, accuracy

# Unchanging values

# Ground Station Requirements
# deg/deg/m - Madrid, Spain - For DSN
gs_ll = np.matrix([[40.416775], [356.29621]]) * u.deg
gs_a = np.matrix([[648.0]]) * u.m
GS = np.matrix([[gs_ll], [gs_a]])
GS_repeat = 4 * u.day

# EXAMPLE Ground Station Requirements
# deg/deg/m - Sacramento, CA
ex_gs_ll = np.matrix([[38.575764], [-121.478851]]) * u.deg
ex_gs_a = np.matrix([[7]]) * u.m
ex_GS = np.matrix([[ex_gs_ll], [ex_gs_a]])

# Target Requirements
# Define path that encircles california through lat/lon coordinates
Target = np.load('orbital_analyses/CA_Coords.npy')

# Subsystem Requirements
###############################################################################

# General
Mass = 1.33 * u.kg
Cost = (500000)
Vol = 1000000 * (u.mm ** 3)

# Power System
Batt_min_percent = 0.2
Batt_max_capacity = 38.0 * u.Wh

# C&DH System
Data_max_capacity = 2.0 * u.Gbyte
Data_min_store = 20.0 * u.Mbyte  # minimum data threshold to downlink
Data_max_pass = 0.0 * u.byte  # We don't know yet

# %% 2nd Iteration and on...

# Mission Requirements
###############################################################################
# TODO: Machine Learning :o


# Subsystem Requirements
###############################################################################
# Replaced by part data
# General
#IVol = Parts_List['Chassis']['volume']


#Data_max_capacity = Parts_List['Payload']['storage'] + Parts_List['Processor']['storage']














