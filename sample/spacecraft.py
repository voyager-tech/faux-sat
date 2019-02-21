# Utilized Modules
import sample.helpers
import sample.propagate
import sample.constants
from sample.orbital_analyses import u
import numpy as np

# Functional Definitions/Functions
# TYPE = constants.TYPE

# =============================================================================
# def array_converter(init):
#     if helpers.step == 0:
#         import inspect
#         paramater_value = inspect.getfullargspec(init)[0]
#         paramater_name = inspect.getfullargspec(init)[3]
#         print(paramater_value)
#         print(paramater_name)
#
#         def wrapper(self, *args):
#             for i in range(len(paramater_value) - 1):
#                 if type(paramater_value[i + 1]) != constants.TYPE:  # Handing for nonetype?
#                     paramater_value[i + 1] = np.asarray(paramater_value[i + 1])
#                 setattr(self, paramater_name[i + 1], paramater_value[i + 1])
#             init(self, *args)
#         return wrapper
# =============================================================================

# Subsystem Definitions


class ADCSubsystem:
    @array_converter
    def __init__(self, component=None, attitude=None, angular_velocity=None,
                 angular_acceleration=None, pointing_vector=None, **kwargs):

        self.component = component
        self.attitude = attitude
        self.angular_velocity = angular_velocity
        self.angular_acceleration = angular_acceleration
        self.pointing_vector = pointing_vector


class AstrodynamicsSubsystem:
    # @array_converter
    def __init__(self, keplerian=None, radius=None, velocity=None,
                 acceleration=None, **kwargs):
#        global TYPE
#        Size = (len(locals()) + len(kwargs) - 2)  # minus self & **kwargs
#        # self.Size = Size  # makes accessible by calling class.Size
#        if helpers.step == 0:
#            paramater_list = [keplerian, radius, velocity, acceleration]
#            for i in range(4):
#                if type(paramater_list[i]) != TYPE:
#                    paramater_list[i] = np.asarray(paramater_list[i])

        self.keplerian = keplerian
        self.radius = radius
        self.velocity = velocity
        self.acceleration = acceleration
        self.__dict__.update(kwargs)


class CDHSubsystem:

    def __init__(self, component=None, stored_data=None, **kwargs):


class CommunicationSubsystem:

    def __init__(self, component=None, to_downlink=None, pointing_vector=None,
                 **kwargs):


class PowerSubsystem:

    def __init__(self, component=None, stored_power=None, power_draw=None,
                 **kwargs):


class PropulsionSubsystem:

    def __init__(self, component=None, fuel_remaining=None, delta_v=None,
                 pointing_vector=None, **kwargs):


class StructureSubsystem:

    def __init__(self, component=None, **kwargs):


class ThermalSubsystem:

    def __init__(self, component=None, temperature=None, **kwargs):


class Time:

    def __init__(self, julian_date=None, gregorian_date=None, **kwargs):


class UserSubsystem:

    def __init__(self, **kwargs):


class Spacecraft:

    def __init__(self, ADCSubsystem=None, AstrodynamicsSubsystem=None,
                 CDHSubsystem=None, CommunicationSubsystem=None,
                 PropulsionSubsystem=None, StructureSubsystem=None,
                 ThermalSubsystem=None, Time=None, UserSubsystem=None,
                 **kwargs):
        self.ADC = ADC()
        self.Astrodynamics = Astrodynamics()
        self.CDH = CDH()
        self.Communication = Communication()
        self.Propulsion = Propulsion()
        self.Structure = Structure()
        self.Thermal = Thermal()
        self.Time = Time()
        self.UserSubsystem = UserSubsystem()


        


