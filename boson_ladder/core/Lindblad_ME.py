from sympy import \
    Add
from sympy.physics.secondquant import \
    Dagger
from .do_commutator import \
    do_commutator
from ..utils.normal_ordering import \
    normal_order as NO
from ..utils.expectation_values import \
    expval
    
__all__ = ["dissipator_trace"]

def dissipator_trace(O, A, normal_order=True):
    """
    tr(D(O) * A) where D(O) is the Lindblad dissipator 
    and A is a polynomial in the ladder operators.
    """
    O = O.expand()
    if isinstance(O, Add):
        O = [arg for arg in O.args]
    else:
        O = [O]
    
    comm = do_commutator
    
    out = 0
    for k, O_k in enumerate(O):
        Od_k = Dagger(O_k)
        out += Od_k * comm(A, O_k)
        out += comm(Od_k, A) * O_k
        for O_l in O[k+1:]:
            Od_l = Dagger(O_l)
            out += Od_k * comm(A, O_l)
            out += comm(Od_k, A) * O_l
            out += Od_l * comm(A, O_k)
            out += comm(Od_l, A) * O_k
    
    if normal_order:
        out = NO(out)
    return expval(out)
