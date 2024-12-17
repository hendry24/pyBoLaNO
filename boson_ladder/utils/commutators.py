from sympy import \
    Mul, \
    Number, \
    srepr, \
    KroneckerDelta
from sympy.physics.secondquant import \
    CreateBoson, \
    AnnihilateBoson, \
    Commutator
from .error_handling import InvalidTypeError
    
__all__ = ["expand_AB_C",
           "expand_A_BC"]

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

def _isolate_bracket(comm):
    """
    Isolate the commutator bracket from the left and right factors.
    Must only be called if a commutator bracket is present.
    
    Parameters
    ----------
    
    comm : sympy.physics.secondquant.Commutator or Mul containing a Commutator
    
    """
    
    if isinstance(comm, Commutator):
        left_factor = Number(1)
        right_factor = Number(1)
        comm = comm
    
    elif isinstance(comm, Mul):
        if "Commutator" not in srepr(comm):
            msg = "Input is Mul but does not contain Commutator."
            raise ValueError(msg)
        
        comm_idx = 0
        for arg in comm.args:
            if isinstance(arg, Commutator):
                """
                Power of a commutator should not occur.
                """
                break
            else:
                comm_idx += 1

        left_factor = Mul(*comm.args[:comm_idx])
        right_factor = Mul(*comm.args[comm_idx+1:])
        comm = comm.args[comm_idx]     # gets assigned last due to the variable assignment.
        
    else:
        raise InvalidTypeError([Commutator, Mul], comm)
            
    return left_factor, comm, right_factor

def _treat_Kron(q):
    """
    Due to the use of `sympy.Symbol` in the bosonic ladder objects, 
    commutator `[b_i,bd_j]` will be a Kronecker delta since SymPy
    cannot tell if the symbols `i` and `j` here are different. This
    function finishes the job.
    
    Parameters
    ----------
    
    q : sympy.functions.special.KroneckerDelta or Mul containing KroneckerDelta
        Object containing the Kronecker delta object.
    """
    
    if isinstance(q, (Number, CreateBoson, AnnihilateBoson)):
        return q
    
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
                    return Number(0)
            out.append(arg)
        return Mul(*out)
    
    else: 
        raise InvalidTypeError([KroneckerDelta, Mul], 
                               type(q))