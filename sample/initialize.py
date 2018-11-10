# Modules Utilized
import numpy as np

# Script to be used for setting initial state as well as component data
# Also needs to handle user analysis/subsystem/component creation

# TODO: Find a good way to have users input all of this data
# TODO: Find a good way to represent data in system (Classes?)

# %% Initial Orbital Values
radius_i = np.array([[],      # i component of radius vector
                     [],      # j component of radius vector
                     []])     # k component of radius vector

velocity_i = np.array([[],    # i component of velocity vector
                       [],    # j component of velocity vector
                       []])   # k component of velocity vector
# OR
keplerian_i = np.array([[],   # Semi-Major Axis (sma)
                        [],   # Eccentricity (ecc)
                        [],   # Inclination (inc)
                        [],   # Right Ascension of the Ascending Node (RAAN)
                        [],   # Argument of Periapsis (AOP)
                        []])  # True Anomaly (TA)
# AND
refrence_frame_i = None       # Choose from a list of supported refrence frames

# %% Time
# Gregorian Date / Julian Date / Local Time + timezone


# %% Spacecraft Attitude


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
