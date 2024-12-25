from sympy import (
    I,
    Number,
    Derivative,
    Symbol,
    Equality
)
from sympy.physics.secondquant import (
    Dagger
)
from .commutators import (
    NO_commutator
)
from .normal_ordering import (
    normal_ordering
)
from ..utils.expval import (
    _expval
)

############################################################

__all__ = ["Hamiltonian_trace",
           "dissipator_trace",
           "LME_expval_evo"]

############################################################

def Hamiltonian_trace(H, A):
    """
    `tr([H,rho]A) = <[A,H]>` where `rho` 
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
    A = A.expand()
    
    out = NO_commutator(A, H)
    
    return _expval(out)

############################################################

def dissipator_trace(O, A, P = None):
    """
    `tr(D(O, P)[rho] * A)` where `rho` is the density matrix.
    
    Parameters
    ----------
    
    O : sympy.Expr
        The operator making up the Liouvillian superoperator 
        in the Lindblad form, also known as the Lindblad 
        dissipator, defined as
            
            `D(O,P)[rho] = O*rho*Pd - 0.5*{Pd*O, rho}`
        
        where `Pd` is the Hermitian conjugate of P (another
        argument of this function), `rho` is
        the system's density matrix, and {.,.} is the
        anticommutator. 
        
    A : sympy.Expr
        The operator to use in the trace. 
        
    P : sympy.Expr, default: None
        The other operator making the dissipator. If not
        specified, then `P=O`. 
    
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
    if P is None:
        P = O
    else:
        P = P.expand()
    A = A.expand()
    
    comm = NO_commutator
    
    Pd = Dagger(P)
    
    out = comm(Pd, A)*O / Number(2)
    out += Pd*comm(A, O) / Number(2)
    
    out = normal_ordering(out.expand())
        
    return _expval(out)

############################################################

def LME_expval_evo(H, D, A, hbar_is_one=True):
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
        of lists of two or three elements. The first element is the
        multiplying scalar, which can be a `sympy.Expr`. The 
        second element is the operator defining the Lindblad
        dissipator. Optionally, the third element is another operator
        defining the dissipator alongside the second element.
        
    A : sympy.Expr
        The operator to calculate the expectation value evolution
        of.
        
    normal_order : bool, default: True
        Whether to normal order the result.
        
    hbar_is_one : bool, default: True
        Whether hbar is omitted in the Hamiltonian trace.
    
    Returns
    -------
    
    out : sympy.Equality
        The evolution equation.
    """
    
    RHS = Hamiltonian_trace(H, A)
    RHS *= -I if hbar_is_one else -I/Symbol(r"hbar")
                                    # Using sympy.physics.quantum.hbar 
                                    # seems to be meddlesome since it
                                    # is not a Number. 
    RHS = RHS.expand()
    
    for D_k in D:
        if len(D_k) == 2:
            D_k = D_k + [None]
        RHS += (D_k[0]*dissipator_trace(O = D_k[1], 
                                        A = A, 
                                        P = D_k[2])).expand()
    
    return Equality(Derivative(_expval(A), Symbol(r"t")),
                    RHS)