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

def Hamiltonian_trace(H, A, normal_order=True, hbar_is_one=True):
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
        
    hbar_is_one : bool, default: True
        Whether hbar is omitted. 
    
    Returns
    -------
    
    out : sympy.Expr
        The Hamiltonian trace, which appears when the Lindblad
        master equation is used to calculate the evolution of 
        some expectation value.
    
    """
    
    scal = -I if hbar_is_one else -I/Symbol(r"hbar")
    out = scal* do_commutator(A, H)
    
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
    
    comm = do_commutator
    
    Od = Dagger(O)
    out = comm(Od, A)*O - Od*comm(A, O)
    out = (out/Number(2)).expand()
        
    if normal_order:
        out = normal_ordering(out)
        
    out = (out/Number(2)).expand()
        
    return _expval(out)

def LME_expval_evo(H, D, A, normal_order = True, hbar_is_one=True):
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
        
    hbar_is_one : bool, default: True
        Whether hbar is omitted in the Hamiltonian trace.
    
    Returns
    -------
    
    out : sympy.Equality
        The evolution equation.
    """
    t = Symbol("t")
    RHS = Hamiltonian_trace(H, A,
                            normal_order=normal_order,
                            hbar_is_one = hbar_is_one)
    
    for D_k in D:
        RHS += D_k[0]*dissipator_trace(D_k[1], A, 
                                       normal_order=normal_order)
            
    if normal_order:
        RHS = normal_ordering(RHS)
    
    return Equality(Derivative(_expval(A), t),
                    RHS)