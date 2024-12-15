from sympy import \
    Mul, \
    Integer, \
    srepr
from sympy.physics.secondquant import \
    AnnihilateBoson, \
    CreateBoson, \
    Commutator, \
    B, \
    Dagger
    
__all__ = ["ops",
           "is_ladder",
           "is_ladder_contained"]

def ops():
    b = B(0)
    return b, Dagger(b)

def is_ladder(q):
    return isinstance(q, (AnnihilateBoson, CreateBoson))

def is_ladder_contained(q):
    return ("AnnihilateBoson" in srepr(q) or "CreateBoson" in srepr(q))

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