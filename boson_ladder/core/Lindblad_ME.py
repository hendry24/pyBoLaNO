from sympy import \
    Add, \
    I, \
    Number, \
    Derivative, \
    Symbol,\
    Equality
from sympy.physics.secondquant import \
    Dagger
from .commutator import \
    do_commutator
from .normal_order import \
    normal_ordering
from ..utils.expval import \
    _expval
    
__all__ = ["Hamiltonian_trace",
           "dissipator_trace",
           "LME_expval_evo"]

def Hamiltonian_trace(H, A, normal_order=True):
    """
    `tr(-i*[H,rho]A) = <-i[A,H]>` where `rho` 
    is the density matrix and `[.,.]` is the commutator.
    
    Parameters
    ----------
    
    H : sympy.Expr
        The Hamiltonian.
        
    A : sympy.Expr
        The operator to use in the trace.
        
    normal_order : bool, default: True
        Whether to normal-order the result.
    
    Returns
    -------
    
    out : sympy.Expr
        The Hamiltonian trace, which appears when the Lindblad
        master equation is used to calculate the evolution of 
        some expectation value.
    
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
    
    return _expval(out)

def dissipator_trace(O, A, normal_order=True):
    """
    `tr(D(O)[rho] * A)` where `rho` is the density matrix.
    
    Parameters
    ----------
    
    O : sympy.Expr
        The operator making up the Liouvillian superoperator 
        in the Lindblad form, also known as the Lindblad 
        dissipator, defined as
            
            `D(O)[rho] = O*rho*Od - 0.5*{Od*O, rho}`
        
        where `Od` is the Hermitian conjugate, `rho` is
        the system's density matrix, and {.,.} is the
        anticommutator.
        
    A : sympy.Expr
        The operator to use in the trace. 
    
    normal_order : bool, default: True
        Whether to normal-order the result.
    
    Returns
    -------
    
    out : sympy.Expr
        The dissipator trace, which appears when the Lindblad
        master equation is used to calculate the evolution of 
        some expectation value.
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
        
    return _expval(out)

def LME_expval_evo(H, D, A, normal_order = True):
    """
    Calculate the evolution of the expectation value
    of `A`, of the system described by the Lindblad master
    equation (LME):

        `d/dt expval(A) = Hamiltonian_trace(H, A) + sum_k D_k[0] * dissipator_trace(D_k[1], A)`
    
    for `D_k` in `D`.
    
    Parameters
    ----------
    
    H : sympy.Expr
        The Hamiltonian.
        
    D : list
        The Lindblad dissipators, specified as a nested list
        of lists of two elements. The first element is the
        multiplying scalar, which can be a `sympy.Expr`. The 
        second element is the operator defining the Lindblad
        dissipator.
        
    A : sympy.Expr
        The operator to calculate the expectation value evolution
        of.
        
    normal_order : bool, default: True
        Whether to normal order the result.
    
    Returns
    -------
    
    out : sympy.Equality
        The evolution equation.
    """
    t = Symbol("t")
    RHS = Hamiltonian_trace(H, A)
    
    for D_k in D:
        RHS += D_k[0]*dissipator_trace(D_k[1], A)
            
    if normal_order:
        RHS = normal_ordering(RHS)
    
    return Equality(Derivative(_expval(A), t),
                    RHS)