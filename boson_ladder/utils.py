from sympy import srepr, Integer, Mul
from sympy.physics.secondquant import Commutator, CreateBoson, AnnihilateBoson, B, Dagger

def is_laddder(q):
    return isinstance(q, (AnnihilateBoson, CreateBoson))

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

def ops():
    b = B(0)
    return b, Dagger(b)