from sympy import \
    Mul, \
    Integer, \
    srepr
from sympy.physics.secondquant import \
    Commutator
    
__all__ = ["expand_AB_C",
           "expand_A_BC"]

def isolate_bracket(comm):
    """
    Isolate the bracket from the left and right factors.
    """
    
    if isinstance(comm, Commutator):
        left_factor = Integer(1)
        right_factor = Integer(1)
        comm = comm
    
    else:
        comm_idx = 0
        for q in comm.args:
            if "Commutator" in srepr(q):
                break
            else:
                comm_idx += 1

        left_factor = Mul(*comm.args[:comm_idx])
        right_factor = Mul(*comm.args[comm_idx+1:])
        comm = comm.args[comm_idx]     # this gets assigned last.
            
    return left_factor, comm, right_factor

def expand_AB_C(A,B,C):
    """
    [AB,C] = A[B,C] + [A,C]B
    """
    return A*Commutator(B,C) + Commutator(A,C)*B

def expand_A_BC(A,B,C):
    """
    [A,BC] = [A,B]C + B[A,C]
    """
    return Commutator(A,B)*C + B*Commutator(A,C)