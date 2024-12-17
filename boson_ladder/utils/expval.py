from sympy import \
    Add, \
    Mul, \
    Pow, \
    Expr, \
    Number, \
    sympify, \
    latex
    
from .operators import \
    is_ladder, \
    is_ladder_contained

class _expval(Expr):
    """
    The expectation value object. Not accessed by
    the user; only for presentation of results. 
    """
    def __new__(cls, *args, **kwargs):
        """
        Only accepts the operator inside the braket. 
        No scalars, no additions. Can also be unity
        for easier time dealing with sums with scalars.
        """

        try:
            q = args[0]
        except:
            q = kwargs.get("q")

        q = sympify(q)
        
        return super().__new__(cls, q)
    
    def __init__(self, q):
        self.q = sympify(q)
        self.is_oper = not(self.q.is_number)
            
    def _get_latex(self):
        if not(self.is_oper):
            return latex(1)
        
        bra = r"{\left\langle "
        # {} to ensure correct exponentiation. 
        # Space is necessary.
        oper = latex(self.q)
        ket = r" \right\rangle}"
        return bra+oper+ket
    
    def __mul__(self, other):
        # Skip 1 explicitly during multiplication
        if not(self.is_oper):
            return other
        if isinstance(other, _expval) and not(other.is_oper):
            return self
        return Mul(self, other)

    def _eval_expand_multinomial(self, **hints):
        # This controls how the term expands in multinomial operations
        if not(self.is_oper):
            return 1
        return self
    
    def _latex(self, printer):
        return self._get_latex()
    
    def _repr_latex(self, printer):
        return self._get_latex()
    
    def __str__(self):
        return self._repr_latex(None)

    def __repr__(self):
        return self._repr_latex(None)
    
def _expval_sum(q):
    """
    Get the expectation value for a 
    scalar-weighted sum of expectation values.
    No a*(b+c) form allowed.
    """
    
    def _treat_Mul(_q):
        oper = []
        scalar = []
        for _arg in _q.args:
            if (isinstance(_arg, Pow) \
                    and is_ladder(_arg.args[0])) \
                or is_ladder(_arg):
                oper.append(_arg)
            else:
                scalar.append(_arg)
                
        return Mul(*scalar) * _expval(Mul(*oper))

    out = Number(0)
    
    if not(is_ladder_contained(q)):
        return q
    
    elif isinstance(q, (Pow)) or \
        is_ladder(q):
        return _expval(q)
    
    elif isinstance(q, Add):
        for qq in q.args:
            if not(is_ladder_contained(qq)):
                out += qq
            
            elif isinstance(qq, Mul):
                out += _treat_Mul(qq)
                
            else:
                out += _expval(qq)
    
    else: # Mul
        out += _treat_Mul(q)
            
    return out