# Defines types of flags that can be instanciated by the state machine

# TODO: Add DocStrings to all classes
# TODO: Possible combine into one class with variable arguments


# Flag (basic - only represents if a function should be run)
class FlagBasic:
    """ DocString here """
    def __init__(self, boolean):
        self.bool = boolean


# Array Flag (refrences a single number as an attribute of the flag)
class FlagValue:

    def __init__(self, boolean, value):
        self.bool = boolean
        self.value = value


# Array Flag (refrences an array as an attribute of the flag)
class FlagArray:

    def __init__(self, boolean, array):
        self.bool = boolean
        self.array = array


# Object Flag (Refrences an object as an attribute of the flag)
class FlagObject:

    def __init__(self, boolean, objects):
        self.bool = boolean
        self.obj = objects
