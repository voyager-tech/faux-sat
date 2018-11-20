# Defines types of flags that can be instanciated by the state machine


# Now can use any additional value with the boolean from just one class
class FlagValue:
    """ DocString here """
    def __init__(self, boolean, input_value=None, return_value=None, **kwargs):
        self.bool = boolean
        self.input_value = input_value
        self.return_value = return_value
