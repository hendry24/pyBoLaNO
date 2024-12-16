from sympy import \
    Mul, \
    Integer, \
    srepr, \
    KroneckerDelta
from sympy.physics.secondquant import \
    Commutator
from .operators import \
    is_ladder
    
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

def _treat_Kron(q):
    """
    q is Kronecker Delta or Mul.
    """
    if isinstance(q, KroneckerDelta):
        return 1 if (q.args[0]==q.args[1]) else 0
    elif isinstance(q, Mul):
        if "KroneckerDelta" not in srepr(q):
            return q
        for arg in q.args:
            out = []
            if isinstance(arg, KroneckerDelta):
                if q.args[0]==q.args[1]:
                    continue
                else:
                    return Integer(0)
            out.append(arg)
        return Mul(*out)
    else: 
        msg = "Expected [KroneckerDelta] or [Mul], "
        msg += f"got [{type(q)}] instead."
        raise ValueError(msg) 