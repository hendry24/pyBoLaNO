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

def _break_comm(A, B):
    """
    To avoid SymPy evaluating the commutatotr to a Kronecker
    delta.
    """
    return (A*B-B*A).expand()

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
    """
    [AB,CD] = A[B,C]D + [A,C]BD + CA[B,D] + C[A,D]B
    """
    return A*Commutator(B,C)*D + Commutator(A,C)*B*D \
            + C*A*Commutator(B,D) + C*Commutator(A,D)*B 
            
############################################################
            
def NO_commutator(A, B):
    """
    Return the normal-ordered equivalent of the commutator 
    of two arbitrary polynomials of bosonic ladder operators.
    
    Parameters
    ----------
    
    A : sympy.Expr
        Operator in the left-hand slot of the commutator bracket.
        
    B : sympy.Expr
        Operator in the right-hand slot of the commutator bracket.
        
    Returns
    ------- 
    
    out : sympy.Expr
        Normal-oredered commutator between A and B.
    """
    
    return normal_ordering(_break_comm(A, B))