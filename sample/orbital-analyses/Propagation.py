# Utilized Modules
import numpy as np
from poliastro.twobody.propagation import cowell
from poliastro.twobody import Orbit
from poliastro.bodies import Earth
from astropy import units as uu
from datetime import datetime, timedelta

# Utilized Scripts
from Orbital_Analyses.Transform_State import Gregorian2JD

from Orbital_Analyses import u


def Prop_Cowell(rad, vel, time, step, Ssize):
    # Input Radius and Velocity vectors in a fixed reference frame (ECEF)
    # If you have Keplerian elements, must convert before using this script
    # Use Keplerian2Cartesian for that

    # Determine if necessary
    # TODO: raised exception for non 3x1 rad and vel entries
    # TODO: raised exception for non 6x1 GD entries or 1x1 JD entries

    # Time can be in either Julian or Gregorian Date (will check and convert)
    # Step is defined in minutes (for now)

    # Make sure vectors are in matrix format not arrays
    pc_rad = np.matrix.tolist(np.asmatrix(rad)) * uu.km
    pc_vel = np.matrix.tolist(np.asmatrix(vel)) * (uu.km / uu.s)

    # Convert time to matrix format
    pc_time = np.asmatrix(time)
    # Check what type of starting time and convert if necessary to Julian Date
    # Time defined in UTC (for now)
    if np.size(pc_time, axis=0) == 1:
        if np.size(pc_time, axis=1) == 1:  # time = Julian
            pc_JD = pc_time
        if np.size(pc_time, axis=1) != 1:  # time = Gregorian (need transpose)
            pc_GD = np.transpose(pc_time)
            pc_JD = Gregorian2JD(pc_GD)
    if np.size(pc_time, axis=0) == 6:  # time = Gregorian (no transpose)
        pc_GD = pc_time
        pc_JD = Gregorian2JD(pc_GD)
    pc_JD = np.asmatrix(pc_JD)
    pc_GD = np.asmatrix(pc_GD)

    # Initialize Arrays to put calculated values into
    if step == 0:
        orbit_r = np.asmatrix(np.zeros((3, 1), dtype=np.float64))
        orbit_v = np.asmatrix(np.zeros((3, 1), dtype=np.float64))
        orbit_jd = np.asmatrix(np.zeros((1, 1), dtype=np.float64))
        orbit_gd = np.asmatrix(np.zeros((6, 1), dtype=np.float64))
    else:
        orbit_r = np.asmatrix(np.zeros((3, step), dtype=np.float64))
        orbit_v = np.asmatrix(np.zeros((3, step), dtype=np.float64))
        orbit_jd = np.asmatrix(np.zeros((1, step), dtype=np.float64))
        orbit_gd = np.asmatrix(np.zeros((6, step), dtype=np.float64))

    # Define initial orbit state for propagator
    pc_initial = Orbit.from_vectors(Earth, pc_rad, pc_vel, pc_JD)

    # Convert step size to seconds
    pc_secs = (Ssize.to(u.sec)).magnitude

    # Define initial date with datetime
    time_ms = np.asmatrix(np.zeros((7, 1), dtype=np.float64))
    time_ms[0:5, 0] = time[0:5, 0]
    time_ms[5, 0] = np.floor(time[5, 0])
    time_ms[6, 0] = (np.mod(time[5, 0], 1) * 1e6)
    time_ms = time_ms.astype(int)
    Itime = datetime(time_ms[0, 0], time_ms[1, 0],
                     time_ms[2, 0], time_ms[3, 0],
                     time_ms[4, 0], time_ms[5, 0], time_ms[6, 0])

    # TODO: Prop from step 3-4 then 4-5 rather than 1-4, 1-5
    # Second method takes very long time but doesn't propagate errors
    # first method has error propagation, but not sure to what degree
    # TODO: ALter nsteps dependent upon tof value

    # Propagate for defined step length
    if step == 0:
        orbit_r[:, step] = (pc_rad.value) * u.km
        orbit_v[:, step] = (pc_vel.value) * (u.km / u.sec)
        orbit_gd[:, step] = pc_GD
        orbit_jd[:, step] = pc_JD
    else:
        for t in range(step):
            # Propagate using Cowell's Method
            if t == 0:
                orbit_r[:, t] = pc_rad.value
                orbit_v[:, t] = pc_vel.value
                orbit_gd[:, t] = pc_GD
                orbit_jd[:, t] = pc_JD
            else:
                orbit_val = cowell(pc_initial, (t * pc_secs))
                # Set rad, vel values into pre-initialized vectors
                # orbit_r / orbit_v are defined in ECEF
                orbit_r[:, t] = orbit_val[0]
                orbit_v[:, t] = orbit_val[1]
                # Elapsed time of current step past the starting time
                elap = (t * pc_secs)  # time of step in elapsed seconds
                # Add elapsed time to initial for current step time (Stime)
                add_sec = timedelta(seconds=elap)
                Stime = Itime + add_sec
                np_Stime = np.datetime64(Stime)
                GD_Stime = np.asmatrix(np.zeros((6, 1), dtype=np.float64))
                GD_Stime[0, 0] = np_Stime.astype(object).year
                GD_Stime[1, 0] = np_Stime.astype(object).month
                GD_Stime[2, 0] = np_Stime.astype(object).day
                GD_Stime[3, 0] = np_Stime.astype(object).hour
                GD_Stime[4, 0] = np_Stime.astype(object).minute
                GD_Stime[5, 0] = (np_Stime.astype(object).second +
                                  (np_Stime.astype(object).microsecond * 1e-6))
                #  JD, GD values into pre-initialized vectors
                orbit_gd[:, t] = GD_Stime
                orbit_jd[:, t] = Gregorian2JD(orbit_gd[:, t])

            # Add Units from pint
            orbit_r[:, t] = orbit_r[:, t] * u.km
            orbit_v[:, t] = orbit_v[:, t] * (u.km / u.s)

    return orbit_r, orbit_v, orbit_jd, orbit_gd
