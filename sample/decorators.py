# Utilized Modules
import numpy as np


# Decorators to do with input validation
def gregorian_date_validation(func):
    def wrapper(*args, **kwargs):
        # If input is array_like then convert to np.ndarray
        gd = np.asarray(args[0])
        # Check to see if all input values are numbers
        try:
            np.all(np.isnan(gd))
        except TypeError:
            raise TypeError('Input array should be composed of only numbers')
        # Raise exceptions for incorrectly arranged array inputs
        # Handling for 1-D arrays
        try:
            gd.shape[1]
        except IndexError:
            if gd.shape[0] != 6:
                raise IndexError('Input array should be 6 entries long')
        else:
            # Handling for array that needs to be transposed to (n, 6)
            if gd.shape[1] != 6:
                raise Exception('Input array should be arranged (n, 6)')
            # Handling for (6, 6) array by checking against max possible values
            elif gd[0, 1] > 12 or gd[0, 2] > 31 or gd[0, 5] > 60:
                raise Exception('Input array should be arranged (n, 6)')
        # Run decorated function
        conversion = func(*args, **kwargs)
        return conversion
    return wrapper
