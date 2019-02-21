# Modules Used
from sample.orbital_analyses import Transform_State as TS
from sample.orbital_analyses import Transform_Coordinate as TC

# Classes for subsystems consisting of methid functions defining different
# function/analysis modes

# Subsystem Definitions

# TODO: Make sure all info determined here is stored and accessable in classes for the state machine
# Function: Something that alters the state of the spacecraft (Spacecraft Mode)
# Analysis: Something that returns data about the current spacecraft state
# TODO: flag object should be tuple of all possible flags from satte machine, pull needed touple per subsection

# TODO: How to save intermediate values (Not in vehicle_state) to be accessable by all other analyses that need them?
# ^^ Save back into flags to be called later and values accessed then.


class ADC:
    """ DocString here """
    def __init__(self, flag, vehicle_state, **kwargs):
        self.flag = flag
        self.vehicle_state = vehicle_state

    # Subsystem Functions (Modes)
    # Pointing to New Direction (FlagArray)
    # Maintain Pointing Direction (FlagArray)
    # Detumble? (FlagBasic)
    # TODO: Account for non-reaction wheel based systems
    # 

    # Subsystem Analyses
    # Current Angular Momentum (FlagBasic)
    # Current Pointing Direction (FlagBasic)
    # Stored Wheel momentum? (FlagBasic)
    # 


class Astrodynamics:
    """ DocString here """
    def __init__(self, flag, vehicle_state, **kwargs):
        self.flag = flag
        self.vehicle_state = vehicle_state

    # Subsystem Functions (Modes)
    # Targeting (FlagArray) (Requires solvers and optimizers)
    # Achieve (FlagArray)
    # Optimize (Flag?) (Possibly a decorator on the command sequence)
    # 

    # Subsystem Analyses
    # TODO: Run through all functions and analyses, and the flag value should determine if rest of analysis is run (0 or 1)
    # TODO: If the value is 1, then any additional agruments of the flag will help determine the rest of the function's activity

    # State Converstions (FlagValue)
    # TODO: Consider moving LLA (spherical) to this section of transforms
    def convert_state(self):
        """ DocString here """
        # From State machine, always cartesian
        rad, vel, acc, state_current = self.vehicle_state['location']
        # What state representation to convert to
        state_desired = self.flag.input_value  # Should be a string value
        # Eithwe keplerian or perifocal so far
        # Set conversion from flag.imput_value and convert that state
        # TODO: Find cleaner way to do this (low)
        if state_desired == "keplerian":
            # needs acceleration handling
            new_kep = TS.Cartesian2Keplerian(rad, vel)
            state_returned = (new_kep)
        if state_desired == "circular":
            # TODO: Function does not yet exist
            # TODO: Needs acceleration handling
            state_returned = (azimuth, elevation)
        if state_desired == "cartesian":
            new_rad = rad
            new_vel = vel
            new_acc = acc
            state_returned = (new_rad, new_vel, new_acc)
        # Place new state values into flag.return_value
        self.flag.return_value = state_returned
        return self.flag

    # Coordinate Transforms (FlagValue)
    # TODO: Move perifocal transformation here
    def convert_coordinate(self):
        """ DocString Here """
        # from state machine, always cartesian and in inertial system (GCRF)
        rad, vel, acc, state_current = self.vehicle_state['location']
        # What coordinate system to convert to
        coordinate_desired = self.flag.input_value  # Should be a string value
        # String Values: ITRF, geocentric, perifocal,
        # Set conversion from flag.imput_value and convert that coordinate
        # TODO: Find cleaner way to do this (low)
        if coordinate_desired == "ITRF":
            
            coordinate_returned = (new_rad, new_vel, new_acc)
        if coordinate_desired == "perifocal":
            # Does not exist yet, needs acceleration handlng too
            new_rad, new_vel, new_acc = TS.Cartesian2Perifocal(rad, vel, acc)
            coordinate_returned = (new_rad, new_vel, new_acc)
        if coordinate_desired == "geocentric":
            
            coordinate_returned = (new_rad, new_vel, new_acc)
        if coordinatedesired == "GCRF":
            new_rad = rad
            new_vel = vel
            new_acc = acc
            coordinate_returned = (new_rad, new_vel, new_acc)
        # Place new state values into flag.return_value
        self.flag.return_value = coordinate_returned
        return self.flag
    #


class CDH:
    """ DocString here """
    def __init__(self, flag, vehicle_state, **kwargs):
        self.flag = flag
        self.vehicle_state = vehicle_state

    # Subsystem Functions (Modes)
    # 

    # Subsystem Analyses
    # Determine stored data (for downlink) (FlagBasic)
    # TODO: Represent all internal system messages being sent/recieved


class Communication:
    """ DocString here """
    def __init__(self, flag, vehicle_state, **kwargs):
        self.flag = flag
        self.vehicle_state = vehicle_state

    # Subsystem Functions (Modes)
    # Downlink data (FlagValue)
    # Uplink data (FlagValue)
    # Listen for comm (FlagBasic)

    # Subsystem Analyses
    # GS Contact Times (FlagBasic)
    # 


class Power:
    """ DocString here """
    def __init__(self, flag, vehicle_state, **kwargs):
        self.flag = flag
        self.vehicle_state = vehicle_state

    # Subsystem Functions (Modes)
    # Generate Power(FlagBasic)
    # 

    # Subsystem Analyses
    # Current Power(FlagBasic)
    # Power Draw (FlagBasic) (Which components using power)
    # Sun Contact Times (FlagBasic)
    # 


class Propulsion:
    """ DocString here """
    def __init__(self, flag, vehicle_state, **kwargs):
        self.flag = flag
        self.vehicle_state = vehicle_state

    # Subsystem Functions (Modes)
    # Impulsive Burn (FlagArray - delta-v vector or time of burn)
    # Begin Finite Burn (FlagArray)
    # End Finite Burn (FlagArray)
    # 

    # Subsystem Analyses
    # Propellant/Delta-V Remaining (FlagBasic)
    # 


class Structure:
    """ DocString here """
    def __init__(self, flag, vehicle_state, **kwargs):
        self.flag = flag
        self.vehicle_state = vehicle_state

    # Subsystem Functions (Modes)
    # 

    # Subsystem Analyses
    # 


class Thermal:
    """ DocString here """
    def __init__(self, flag, vehicle_state, **kwargs):
        self.flag = flag
        self.vehicle_state = vehicle_state

    # Subsystem Functions (Modes)
    # Achieve Temp (FlagValue)
    # Maintain Temp (FlagValue)
    # 

    # Subsystem Analyses
    # Current Temp (FlagBasic)
    # 


class UserDesined:
    """ DocString here """
    def __init__(self, flag, vehicle_state, **kwargs):
        self.flag = flag
        self.vehicle_state = vehicle_state

    # Subsystem Functions (Modes)
    # 

    # Subsystem Analyses
    # 


# Run each step to help determine the state of the spacecraft at the next step
# Make private methods
class General:
    """ DocString here """
    def __init__(self, flag, vehicle_state, **kwargs):
        self.flag = flag
        self.vehicle_state = vehicle_state

    # Subsystem Functions (Modes)
    # Report Data (FlagArray)
    # Plot Data (FlagArray)
    # Error Handling? (Contingency)
    # 

    # Subsystem Analyses
    # In Sun (boolean)
    # In Ground Contact (Point & Area) (boolean)
    # Attitude (array)
    # Always on sensor data (array)
    # 
