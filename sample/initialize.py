# Modules Utilized
import numpy as np
from sample.orbital_analyses import u  # Pint Units

# Script to be used for setting initial state as well as component data
# Also needs to handle user analysis/subsystem/component creation

# TODO: Find a good way to have users input all of this data
# TODO: Find a good way to represent data in system (Classes?)
# TODO: Determine default values for as many points as possible

# %% System Filepaths
# Set the path to faux-sat folder in the variable faux_path (Include faux-sat)
# Example: faux_path = r'\Users\Harley\Documents\faux-sat'
# Be sure to include "r'\'" before first directory
faux_path = r'\Users\Harley\Documents\faux-sat'

# %% Initial Orbital Values
# TODO: Convert all to rad, vel in km and km/sec in the inertial frame
radius_i = np.array([[],      # i component of radius vector
                     [],      # j component of radius vector
                     []])     # k component of radius vector
radius_unit = ()  # u.km, u.m, u.ft, u.mi

velocity_i = np.array([[],    # i component of velocity vector
                       [],    # j component of velocity vector
                       []])   # k component of velocity vector
velocity_unit = ()  # u.km/u.sec, u.m/u.sec, u.ft/u.sec, u.mi/u.hr
# OR
keplerian_i = np.array([[],   # Semi-Major Axis (sma)
                        [],   # Eccentricity (ecc)
                        [],   # Inclination (inc)
                        [],   # Right Ascension of the Ascending Node (RAAN)
                        [],   # Argument of Periapsis (AOP)
                        []])  # True Anomaly (TA)
# AND
refrence_frame_i = ()       # Choose from a list of supported refrence frames

# %% Time
# Gregorian Date / Julian Date / Local Time + timezone
# TODO: Convert to JD
julian_date = ()
# OR
mod_julian = ()
# OR
gregorian_date_UTC = np.array([[],   # Year
                               [],   # Month
                               [],   # Day
                               [],   # Hour
                               [],   # Minute
                               []])  # Second w/ decimal
# OR
local_time = [np.array([[],   # Year
                        [],   # Month
                        [],   # Day
                        [],   # Hour
                        [],   # Minute
                        []]),  # Second w/ decimal
              ()]  # Timezone from list of timezones

# %% Spacecraft Attitude
# includes rate of change if not already included in quaternions
# TODO: Convert to quaternions
quaternion = np.array([[],
                       [],
                       [],
                       []])


# %% Propagator
# Include Perturbations here


# %% Ground Locations
# Both point locations and area locations


# %% Launch Vehicle
# For WAYYY down the line


# %% Environment Models
# Gravity / Atmosphere


# %% Component Data
# Split up by subsection and then by component type (make general input form?)


# %% Component Layout


# %% Additionally have the user do:
# Fill out state machine
# Fill out mission sequence
# Add any user defined systems/functions/frames/...
# General assumptions for missing data in model
