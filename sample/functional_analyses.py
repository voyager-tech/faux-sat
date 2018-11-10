# Classes for subsystems consisting of methid functions defining different
# function/analysis modes

# Subsystem Definitions

# TODO: Make sure all info determined here is stored and accessable in classes for the state machine
# Function: Something that alters the state of the spacecraft (Spacecraft Mode)
# Analysis: Something that returns data about the current spacecraft state


class ADC:

    def __init__(self, flag, vehicle_state, **kwargs):

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

    def __init__(self, flag, vehicle_state, **kwargs):

    # Subsystem Functions (Modes)
    # Targeting (FlagArray) (Requires solvers and optimizers)
    # Achieve (FlagArray)
    # Optimize (Flag?) (Possibly a decorator on the command sequence)
    # 

    # Subsystem Analyses
    # State Converstions (FlagValue)
    # Coordinate Transforms (FlagValue)
    #

class CDH:

    def __init__(self, flag, vehicle_state, **kwargs):

    # Subsystem Functions (Modes)
    # 

    # Subsystem Analyses
    # Determine stored data (for downlink) (FlagBasic)
    # TODO: Represent all internal system messages being sent/recieved

class Communication:

    def __init__(self, flag, vehicle_state, **kwargs):

    # Subsystem Functions (Modes)
    # Downlink data (FlagValue)
    # Uplink data (FlagValue)
    # Listen for comm (FlagBasic)

    # Subsystem Analyses
    # GS Contact Times (FlagBasic)
    # 


class Power:

    def __init__(self, flag, vehicle_state, **kwargs):

    # Subsystem Functions (Modes)
    # Generate Power(FlagBasic)
    # 

    # Subsystem Analyses
    # Current Power(FlagBasic)
    # Power Draw (FlagBasic) (Which components using power)
    # Sun Contact Times (FlagBasic)
    # 


class Propulsion:

    def __init__(self, flag, vehicle_state, **kwargs):

    # Subsystem Functions (Modes)
    # Impulsive Burn (FlagArray - delta-v vector or time of burn)
    # Begin Finite Burn (FlagArray)
    # End Finite Burn (FlagArray)
    # 

    # Subsystem Analyses
    # Propellant/Delta-V Remaining (FlagBasic)
    # 

class Structure:

    def __init__(self, flag, vehicle_state, **kwargs):

    # Subsystem Functions (Modes)
    # 

    # Subsystem Analyses
    # 

class Thermal:

    def __init__(self, flag, vehicle_state, **kwargs):

    # Subsystem Functions (Modes)
    # Achieve Temp (FlagValue)
    # Maintain Temp (FlagValue)
    # 

    # Subsystem Analyses
    # Current Temp (FlagBasic)
    # 

class UserDesined:

    def __init__(self, flag, vehicle_state, **kwargs):

    # Subsystem Functions (Modes)
    # 

    # Subsystem Analyses
    # 


# Run each step to help determine the state of the spacecraft at the next step
# Make private methods
class General:

    def __init__(self, flag, vehicle_state, **kwargs):

    # Subsystem Functions (Modes)
    # Report Data (FlagArray)
    # Plot Data (FlagArray)
    # Error Handling? (Contingency)
    # 

    # Subsystem Analyses
    # In Sun (boolean)
    # In Ground Contact (Point & Area) (boolean)
    # Attitude (array)
    # 
