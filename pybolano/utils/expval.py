from sympy import (
    Add,
    Mul,
    Pow,
    Symbol,
    Number,
    latex
)
from sympy.physics.secondquant import (
    AnnihilateBoson,
    CreateBoson
)
from .operators import (
    is_ladder_contained
)

############################################################

__all__ = []

############################################################

class _expval(Symbol):
    """
    The expectation value object. Only for presentation of
    results in .core.Lindblad_ME, not for users to use.
    
    Parameters
    ----------
    
    q : sympy.Expr
        The quantity to be put inside the braket.
    """
    def __new__(cls, *args, **kwargs):
        """
        Assumes a polynomial in ladder operators. Each term in
        the input Add is properly enclosed in bra-kets
        """
        
        def _braket(q):
            return r"{\left\langle " + latex(q) + r" \right\rangle}"

        def _process(q):
            """
            Assumes q to be one term to enclose in braket. This function
            separates the scalars from the operators enclosed in bra-ket.
            """
            if not(is_ladder_contained(q)):
                return q, Number(1) # Not string, no bra-ket
            
            elif isinstance(q, (Pow, 
                                CreateBoson, 
                                AnnihilateBoson)):
                return Number(1), _braket(q)
            
            elif isinstance(q, Mul):
                scalars = []
                opers = []
                for arg in q.args:
                    if is_ladder_contained(arg):
                        opers.append(arg)
                    else:
                        scalars.append(arg)
                        
                return Mul(*scalars), _braket(Mul(*opers))
        ###
        
        try:
            q = args[0]
        except:
            q = kwargs.get("q")
        
        if q is None:
            return Number(0)
        q = q.expand()
            
        ###
        
        constructor = []
        if isinstance(q, Add):
            q_args = q.args
        else:
            q_args = [q]

        for qq in q_args:
            scal, braket = _process(qq)
            if isinstance(braket, str):
                braket = super().__new__(cls, braket)
            constructor.append(Mul(scal, braket))

        return Add(*constructor)
        
    def __init__(self, q):
        pass