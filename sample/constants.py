# Utilized Modules
import numpy as np
from orbital_analyses import u

# Array type to check against
TYPE = type(np.array(0))

# Spaceflight Constants - Vallado, Fundamentals of Astrodynamics & Applications
EARTH_RADIUS = 6378.137 * u.km
EARTH_GRAV_PARAM = 398600.4418 * ((u.km ** 3) / (u.sec ** 2))
EARTH_ROTATION = 0.0000729211585530 * (u.rad / u.sec)

EARTH_ECC = 0.081819221456
EARTH_ECC2 = 0.006694385000
EARTH_FLATTENING = (1.0 / 298.257)
EARTH_J2 = 0.0010826267
EARTH_J3 = -0.0000025327
EARTH_J4 = -0.0000016196



