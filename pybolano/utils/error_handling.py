__all__ = []

############################################################

class InvalidTypeError(TypeError):
    def __init__(self, expected, got):
        """
        Raise an error for bad input types.
        
        Parameters
        ----------
        
        expected : type or list of type
            Expected object types.
            
        got : type
            Input type.
        """
        if not(isinstance(expected, list)):
            expected = [expected]
        
        msg = f"Expected \n ["
        msg += f"{expected[0]}"
        for obj in expected[1:]:
            msg += "] \n"
            msg += " or \n ["
            msg += f"{obj}"
        msg += f"], \n got \n [{got}] \n instead."
        super().__init__(msg)