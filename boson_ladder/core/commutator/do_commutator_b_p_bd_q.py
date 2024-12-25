from sympy import (
    Number,
    Pow,
    binomial,
    FallingFactorial
)
from sympy.physics.secondquant import (
    AnnihilateBoson,
    CreateBoson,
    AnnihilateFermion,
    CreateFermion
)
from ...utils.operators import (
    get_ladder_attr
)
from ...utils.error_handling import (
    InvalidTypeError
)

############################################################

__all__ = []

############################################################

def _do_commutator_b_p_bd_q(b_p, bd_q):
    """
    [b**p, bd**q] where b is an AnnihateBoson
    object and bd is a CreateBoson object.
    """
    
    ###
    
    # Shortcuts
    
    if b_p.has(CreateFermion, AnnihilateFermion) or\
        bd_q.has(CreateFermion, AnnihilateFermion):
        raise TypeError("Fermionic ladder operators are not accepted.")
    
    if b_p.has(CreateBoson):
        raise InvalidTypeError(AnnihilateBoson, type(b_p))
    if bd_q.has(AnnihilateBoson):
        raise InvalidTypeError(CreateBoson, type(bd_q))
    
    if not(b_p.has(AnnihilateBoson)) \
        or not(bd_q.has(CreateBoson)):
        return Number(0)
    
    sub_b_p, exp_b_p = get_ladder_attr(b_p)
    sub_bd_q, exp_bd_q = get_ladder_attr(bd_q)
    
    if sub_b_p != sub_bd_q:
        return Number(0)
    
    ### 
    
    if isinstance(b_p, AnnihilateBoson):
        b = b_p
        p = 1
    elif isinstance(b_p, Pow) \
        and b_p.has(AnnihilateBoson):
        b, p = b_p.args
    else:
        raise InvalidTypeError([AnnihilateBoson, Pow],
                               type(b_p))
    
    if isinstance(bd_q, CreateBoson):
        bd = bd_q
        q = 1
    elif isinstance(bd_q, Pow) \
        and bd_q.has(CreateBoson):
        bd, q = bd_q.args
    else:
        raise InvalidTypeError([CreateBoson, Pow],
                               type(bd_q))
    
    out = Number(0)
    for k in range(max(0, p-q), p):
        out += binomial(p, k) * FallingFactorial(q, p-k) \
                * bd**(q-p+k) * b**k 
    return out