# Utilized Modules
import numpy as np
from orbital_analyses import u

# Spaceflight Constants - Vallado, Fundamentals of Astrodynamics & Applications
GRAV_PARAMATER = {
                  'SUN': 1.32712440018e11,
                  'MERCURY': 2.2032e4,
                  'VENUS': 3.24859e5,
                  'EARTH': 3.986004418e5,
                  'MOON': 4.9048695e3,
                  'MARS': 4.282837e4,
                  'JUPITER': 1.26686534e8,
                  'SATURN': 3.7931187e7,
                  'URANUS': 5.793939e6,
                  'NEPTUNE': 6.836529e6,
                  'PLUTO': 8.71e2}  # km^3 / sec^2
RADIUS = {
          'SUN': 695000,
          'MERCURY': 2439.7,
          'VENUS': 6051.8,
          'EARTH': 6378.137,
          'MOON': 1738.1,
          'MARS': 3396.2,
          'JUPITER': 71492,
          'SATURN': 60268,
          'URANUS': 25559,
          'NEPTUNE': 24764,
          'PLUTO': 1195}  # km

#  EARTH_RADIUS = 6378.137 * u.km
#  EARTH_MU = 398600.4418 * ((u.km ** 3) / (u.sec ** 2))

# Earth specific constants
EARTH_ROTATION = 0.0000729211585530 * (u.rad / u.sec)
EARTH_ECC = 0.081819221456
EARTH_ECC2 = 0.006694385000
EARTH_FLATTENING = (1.0 / 298.257)
EARTH_J2 = 0.0010826267
EARTH_J3 = -0.0000025327
EARTH_J4 = -0.0000016196



