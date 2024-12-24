from sympy import (
    Mul,
    Number,
    factorial,
    FallingFactorial,
    binomial
)
from sympy.physics.secondquant import (
    AnnihilateBoson,
    CreateBoson
)
from ...utils.operators import (
    get_ladder_attr
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
    Eqs. (4.2), (4.10), (4.34)
    
    The shortcuts in the main function should allow no
    AnnihilateBoson, CreateBoson, or Pow to get here.
    As such, both b and bd must be in q and we can use 
    the .find method to get them. 
    """
    
    b = q.find(AnnihilateBoson)
    bd = q.find(CreateBoson)
    
    if isinstance(q, Mul):
        q_args = q.args
    else:
        q_args = [q]
    
    r = [] if q_args[0].has(AnnihilateBoson) \
        else [Number(0)]    # in case monomial starts with b
    s = []
    for i, qq in enumerate(q_args):
        sub, exp = get_ladder_attr(qq)
        if qq.has(CreateBoson):
            r.insert(0, exp)
        else:
            s.insert(0, exp)
    if len(r) != len(s):    # monomial ends with bd
        s.insert(0, Number(0))
    
    # To make indexing easier, we pad with r_0 and s_0,
    # which do not exist in Blasiak's formulation.
    r.insert(0, Number(0))
    s.insert(0, Number(0))
        
    # Excess
    d = [Number(0)] # d_0 is, however, used in Eq. (4.10).
    sum_val = 0
    for r_m,s_m in zip(r,s):
        sum_val += (r_m-s_m)
        d.append(sum_val)
    
    ###
    
    def _get_S_rsk(r, s, d, k):
        sum_val = Number(0)
        for j in range(k):
            prod_val = Number(1)
            for m in range(1, len(s)):
                prod_val *= FallingFactorial(d[m-1]+j, s[m])
            sum_val += binomial(k, j) * Number(-1)**(k-j) * prod_val
        return 1/factorial(k) * sum_val
    
    ###
    
    if d[-1] >= 0:
        k_lst = range(s[1], sum(s)+1)
    else:
        k_lst = range(r[-1], sum(r)+1)
        
    out = Number(0)
    for k in k_lst:
        out += _get_S_rsk(r,s,d,k) * bd**k * b**k
    
    if d[-1] >= 0:
        out = (bd**d[-1] * out).expand()
    else:
        out = (out * b**(-d[-1])).expand()
        
    return out