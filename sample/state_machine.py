# User defined state determiner per step of simulation

# %% Vehicle State (Stored in Spacecraft class object from spacecraft.py)

# Definitely a dictionary file per step (within an array?)
# Determine preffered representations to convert to before setting state

# Time
# Step (in relation to all otehr simulation steps)
# Current operation running from mission_sequence.py
# Actions being preformed (Determine based on flags produced by state machine) (aka current spacecraft modes)
# From actions / state machine, component states
# Location (radius / velocity / acceleration & refrence frame)
# Attitude (euler angles / quaternions / rates of change)
# General Analyses data from functional_analyses.py


state = {
        'time': time_c,  # Julian Date
        'step': step_c,  # int
        'sequence' : sequence_c,  # str of function name in mission_sequence.py
        'action' : action_c,  # array of flags? or related functions
        'component' : component_c,  # array of boolish values?
        'location' : location_c,  # rad, vel, acc, refrence frame in tuple
        'attitude' : attitude_c,  # euler angles / quaternions + rates in array
        'sensor_data' : sensor_data_c,  # Data from all sensors in array
        'general_data' : general_data_c  # sun / GS contact etc. in array
        }
