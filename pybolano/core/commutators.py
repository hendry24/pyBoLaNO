from sympy.physics.secondquant import (
    Commutator
)
from .normal_ordering import (
    normal_ordering
)

############################################################

__all__ = ["NO_commutator",
           "expand_comm_AB_C",
           "expand_comm_A_BC",
           "expand_comm_AB_CD"]

############################################################

def expand_comm_AB_C(A,B,C):
    """
    [AB,C] = A[B,C] + [A,C]B
    """
    return A*Commutator(B,C) + Commutator(A,C)*B

def expand_comm_A_BC(A,B,C):
    """
    [A,BC] = [A,B]C + B[A,C]
    """
    return Commutator(A,B)*C + B*Commutator(A,C)   

def expand_comm_AB_CD(A,B,C,D):
    return A*Commutator(B,C)*D + Commutator(A,C)*B*D \
            + C*A*Commutator(B,D) + C*Commutator(A,D)*B 
            
############################################################
            
def NO_commutator(A, B):
    """
    Return the normal-ordered equivalent of  commutator 
    of two arbitrary polynomials of bosonic ladder operators.
    
    Parameters
    ----------
    
    A : sympy.Expr
        Operator in the left-hand slot of the commutator bracket.
        
    B : sympy.Expr
        Operator in the right-hand slot of the commutator bracket.
        
    normal_order : bool, default: True
        Whether to normal-order the result. If `True`, then the
        algorithm uses Blasiak's formulae (see Eqs. (4.2) 
        and (4.34) in https://arxiv.org/abs/quant-ph/0507206)
        to evaluate the normal-ordered form of AB and BA, then
        take the difference. If `False', then the algorithm
        recursively uses [AB,CD] = A[B,C]D + [A,C]BD + CA[B,D] + C[A,D]B 
        to simplify the commutator  
        
    Returns
    ------- 
    
    out : sympy.Expr
        Normal-oredered commutator between A and B.
    """
    
    return normal_ordering(Commutator(A,B).doit().expand())