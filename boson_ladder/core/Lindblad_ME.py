from sympy import \
    Add, \
    I, \
    Number, \
    Derivative, \
    Symbol,\
    Equality
from sympy.physics.secondquant import \
    Dagger
from .do_commutator import \
    do_commutator
from .normal_order import \
    normal_ordering
from ..utils.expval import \
    _expval_sum
    
__all__ = ["Hamiltonian_trace",
           "dissipator_trace",
           "moment_evo"]

def Hamiltonian_trace(H, A, normal_order=True, _braket = True):
    """
    tr(-1j*[H,.]A) where -1j*[H,.] is the Hamiltonian
    superoperator. 
    """
    H = H.expand()
    if isinstance(H, Add):
        H = [arg for arg in H.args]
    else:
        H = [H]
        
    comm = do_commutator
    out = Number(0)
    
    for H_k in H:
        out += -I*comm(A, H_k)
    
    if normal_order:
        out = normal_ordering(out)
    
    return _expval_sum(out) if _braket else out 

def dissipator_trace(O, A, normal_order=True, _braket=True):
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
    
    out = Number(0)
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
        out = normal_ordering(out)
        
    out = (out/Number(2)).expand()
        
    return _expval_sum(out) if _braket else out

def moment_evo(H, D, A, normal_order = True):
    t = Symbol("t")
    RHS = Hamiltonian_trace(H, A, _braket=False)
    
    if not(isinstance(D, list)):
        RHS += dissipator_trace(D, A)
    else:
        for k, D_k in enumerate(D):
            if not(isinstance(D_k, list)):
                RHS += dissipator_trace(D_k, A, _braket=False)
            else:
                RHS += D_k[0]*dissipator_trace(D_k[1], A, _braket=False)
            
    if normal_order:
        RHS = normal_ordering(RHS)
        
    return Equality(Derivative(_expval_sum(A), t),
                    _expval_sum(RHS))