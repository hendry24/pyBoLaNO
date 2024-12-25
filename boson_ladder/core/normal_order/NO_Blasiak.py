from sympy import (
    Mul,
    Pow,
    Number,
    factorial,
    FallingFactorial,
    binomial
)
from sympy.physics.secondquant import (
    AnnihilateBoson,
    CreateBoson
)
from ..commutator.do_commutator_b_p_bd_q import (
    _do_commutator_b_p_bd_q
)
from ...utils.operators import (
    is_ladder,
    get_ladder_attr,
    is_ladder_contained
)
from ...utils.error_handling import (
    InvalidTypeError
)
from .NO_find_and_swap import (
    _NO_find_and_swap
)

############################################################

__all__ = []    

############################################################

def _NO_Blasiak(q):
    """
    Normal ordering with the explicit formula derived in 
    Chapter 4 of Blasiak's PhD Dissertation, available at
    https://arxiv.org/abs/quant-ph/0507206.
    
    q is assumed to be a ``boson string'' per Blasiak, i.e.
    a monomial in the bosonic ladder operators. Refer to
    Eqs. (4.2), (4.10), (4.34). 
    
    Input is assumed to contain single subscript.
    """
        
    if not(is_ladder_contained(q)) \
        or is_ladder(q) \
        or isinstance(q, Pow):
        return q
            
    if isinstance(q, Mul):
        q_args = q.args
        if len(q_args) == 2 \
            and q_args[0].has(AnnihilateBoson) \
            and q_args[1].has(CreateBoson):
            return _do_commutator_b_p_bd_q(*q_args) \
                    + q_args[1]*q_args[0]
    else:
        raise InvalidTypeError([AnnihilateBoson, 
                                CreateBoson,
                                Pow,
                                Mul],
                               type(q))
        
    ###
    
    r = [] if q_args[0].has(CreateBoson) \
        else [Number(0)]    # in case monomial starts with b
    s = []
    for qq in q_args:
        sub, exp = get_ladder_attr(qq)
        if qq.has(CreateBoson):
            r.insert(0, exp)
        else:
            s.insert(0, exp)
    if len(r) != len(s):    # monomial ends with bd
        s.insert(0, Number(0))
    
    # To make indexing easier, we pad r and s 
    # with r_0 and s_0,
    # which do not exist in Blasiak's formulation.
    r.insert(0, Number(0))
    s.insert(0, Number(0))
        
    # Excess
    d = [] # d_0 is, however, used in Eq. (4.10).
    sum_val = Number(0)
    for r_m,s_m in zip(r,s):
        sum_val += (r_m-s_m)
        d.append(sum_val)

    ###
    
    def _S_rs(s,d,k):
        """
        Generalized Stirling number, Eq. (4.10). We use
        d instead of r since d is already calculated 
        before this function is called.
        """
        sum_val = Number(0)
        for j in range(k+1):
            prod_val = Number(1)
            for m in range(1, len(s)):
                prod_val *= FallingFactorial(d[m-1]+j, s[m])
            sum_val += binomial(k, j) * Number(-1)**(k-j) * prod_val
        return 1/factorial(k) * sum_val
    
    ###

    b = list(q.find(AnnihilateBoson))[0]
    bd = list(q.find(CreateBoson))[0]
    
    if d[-1] >= 0:
        R,S,D = r,s,d
        k_lst = range(s[1], sum(s)+1)
    else:
        k_lst = range(r[-1], sum(r)+1)
        
        """
        Somehow using the original expression
        in Eq. (4.34) does not work. However,
        it does work when we utilize the symmetry
        property stated in Eq. (4.37).
        """
        R = [Number(0)] + list(reversed(s[1:]))
        S = [Number(0)] + list(reversed(r[1:])) 
        D = []
        sum_val = Number(0)
        for r_m,s_m in zip(R,S):
            sum_val += (r_m-s_m)
            D.append(sum_val)
        
    out = Number(0)
    for k in k_lst:
        out += _S_rs(S,D,k) * bd**k * b**k
        
    if d[-1] >= 0:
        out = (bd**d[-1] * out).expand()
    else:
        out = (out * b**(-d[-1])).expand()
        
    return out