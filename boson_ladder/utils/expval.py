from sympy import \
    Add, \
    Mul, \
    Pow, \
    Symbol, \
    Number, \
    sympify, \
    latex
from .operators import \
    is_ladder, \
    is_ladder_contained

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
        Assumes a polynomial in ladder operators. The 
        returned expectation values are thus a 
        """
        
        def _get_braket(q):
            if is_ladder_contained(q):
                return r"{\left\langle " + latex(q) + r" \right\rangle}"
            else:
                return latex(q)

        def _process(q):
            if not(isinstance(q, (Mul, Pow))):
                return Number(1), _get_braket(q)
                        
            if q.has(Pow):
                q_args = []
                for qq in q.args:
                    if isinstance(qq, Pow):
                        if qq.args[1] > 0:
                            q_args.extend([qq.args[0]]*qq.args[1])
                        else:
                            q_args.append(qq) 
                            # Scalars in the denominator. 
                    else:
                        q_args.append(qq)
            else:
                q_args = q.args
            
            scalars = []
            opers = []
            for arg in q_args:
                if is_ladder(arg):
                    opers.append(arg)
                else:
                    scalars.append(arg)
                    
            return Mul(*scalars), _get_braket(Mul(*opers))
        ###
        
        try:
            q = args[0]
        except:
            q = kwargs.get("q")
        
        if q is None:
            return Number(0)
        q = sympify(q).expand()
            
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
            constructor.append(Mul(*[scal, braket]))

        return Add(*constructor)
        
    def __init__(self, q):
        pass