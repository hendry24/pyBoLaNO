__all__ = []

class InvalidTypeError(TypeError):
    def __init__(self, expected, got):
        """
        Raise an error for bad input types.
        
        Parameters
        ----------
        
        expected : type or list of type
            Expected object types.
            
        got : type:
            Input type.
        """
        if not(isinstance(expected), list):
            expected = [expected]
            
        for i, obj in enumerate(expected): 
            if type(obj) != "type":
                expected[i] = type(obj)
        if type(got) != "type":
            got = type(got)
        
        msg = f"Expected ["
        msg += f"{expected[0]}"
        for obj in expected[1:]:
            msg += "or"
            msg += f"{expected}"
        msg += f"], got [{got}] instead."
        super().__init__(msg)