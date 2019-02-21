# Utilized Modules
import numpy as np
import time
import os

from sample.orbital_analyses.Propagation import Prop_Cowell

from sample.orbital_analyses.Transform_Coordinate import Geodetic2ECFixed
from sample.orbital_analyses.Transform_Coordinate import FK5_J20002ECFixed
from sample.orbital_analyses.Transform_Coordinate import FK5_ECFixed2J2000
from sample.orbital_analyses.Transform_Coordinate import ECFixed2Geodetic

from sample.orbital_analyses.GS_Contact_Times import satContact
from sample.orbital_analyses.Sun_Contact_Times import satIlluminated
from sample.orbital_analyses.Geo_Contact_Times import satOverhead

from sample.con_ops import component_operations

import sample.Requirements as Req
# from orbital_analyses import u

t0 = time.time()
# %%  Orbit Propagation & Analysis

# Initialize matrices to hold calculated rad, vel, JD values
orbit_r = np.asmatrix(np.zeros((3, Req.Steps), dtype=float))
orbit_v = np.asmatrix(np.zeros((3, Req.Steps), dtype=float))
orbit_jd = np.asmatrix(np.zeros((1, Req.Steps), dtype=float))
orbit_gd = np.asmatrix(np.zeros((6, Req.Steps), dtype=float))

# Propagate Orbit by defined number of steps using Cowell's Method
[orbit_r, orbit_v, orbit_jd, orbit_gd] = Prop_Cowell(Req.Rad, Req.Vel,
                                                     Req.GD_UTC, Req.Steps,
                                                     Req.Ssize)

# Illumination Simulation
# INPUT: Inertial
inSun = np.asmatrix(np.zeros((1, Req.Steps), dtype=np.float64))
for i in range(Req.Steps):
    inSun[:, i] = satIlluminated(orbit_r[:, i], orbit_gd[:, i])

# Convert Ground Station to fixed
gs_zero = np.asmatrix(np.zeros((3, 1), dtype=np.float64))
gs_geo = np.matrix([[Req.gs_ll[0, 0].magnitude],
                    [Req.gs_ll[1, 0].magnitude],
                    [Req.gs_a.magnitude]])
gs_fix = Geodetic2ECFixed(gs_geo)

# Ground Station Contact Simulation
# INPUT: Inertial
inContact = np.asmatrix(np.zeros((1, Req.Steps), dtype=int))
for j in range(Req.Steps):
    # Convert ground station to inertial
    gs_inert, gs_vel = FK5_ECFixed2J2000(gs_fix, gs_zero, orbit_gd[:, j])
    inContact[:, j] = satContact(orbit_r[:, j], gs_inert)

# Transform from J2000 to Fixed Frame (EFEC)
# INPUT: Inertial, OUTPUT: Fixed
fixed_r = np.asmatrix(np.zeros((3, Req.Steps), dtype=np.float64))
fixed_v = np.asmatrix(np.zeros((3, Req.Steps), dtype=np.float64))
for i in range(Req.Steps):
    fixed_r[:, i], fixed_v[:, i] = FK5_J20002ECFixed(orbit_r[:, i],
                                                     orbit_v[:, i],
                                                     orbit_gd[:, i])

# Calculate Geodetic Coordinates of the satellite (Lat / Lon / Alt)
# INPUT: Fixed, OUTPUT: Geodetic
sat_geo = np.asmatrix(np.zeros((3, Req.Steps), dtype=int))
for j in range(Req.Steps):
    sat_geo[:, j] = ECFixed2Geodetic(fixed_r[:, j])

# Calculate if Sat is over CA
# INPUT: Geodetic
inCA = np.asmatrix(np.zeros((1, Req.Steps), dtype=int))
for k in range(Req.Steps):
    inCA[:, k] = satOverhead(sat_geo[:, k])

# Run Con-Ops
# Convert Ground Station to Inertial coordinates across steps
gs_zero = np.asmatrix(np.zeros((3, Req.Steps), dtype=np.float64))
sc_zero = np.asmatrix(np.zeros((3, 1), dtype=np.float64))
gs_inert = np.asmatrix(np.zeros((3, Req.Steps), dtype=np.float64))
for l in range(Req.Steps):
    gs_inert[:, l], gs_zero[:, l] = FK5_ECFixed2J2000(gs_fix, sc_zero,
                                                      orbit_gd[:, l])

sc_state, sc_batt_percent, sc_data_percent, boo, sc_vol, sc_mass = component_operations(fixed_r, fixed_v, orbit_jd, inSun, inContact, inCA, gs_inert)

if os.path.exists(r'orbital_analyses\EOPCO4.npy'):
    os.remove(r'orbital_analyses\EOPCO4.npy')
elif os.path.exists(r'EOPCO4.npy'):
    os.remove(r'EOPCO4.npy')

t1 = time.time()
print(t1 - t0)
