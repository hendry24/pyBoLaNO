from sympy import \
    Add, \
    Mul, \
    Pow, \
    Number, \
    KroneckerDelta
from sympy.physics.secondquant import \
    CreateBoson, \
    AnnihilateBoson, \
    Commutator
from .error_handling import InvalidTypeError
    
__all__ = []

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
        if not(comm.has(Commutator)):
            msg = "Input is Mul but does not contain Commutator."
            raise ValueError(msg)
        
        comm_idx = 0
        for arg in comm.args:
            if isinstance(arg, Commutator):
                """
                Power of a commutator should not occur.
                """
                break
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
    
    q : sympy.Expr
        Object containing the Kronecker delta object.
    """
    
    if isinstance(q, Add):
        q = [arg for arg in q.args]
    else:
        q = [q]

    out = []
    for qq in q: 
        if isinstance(qq, (Number, CreateBoson, AnnihilateBoson, Pow)):
            """
            A power of KroneckerDelta should not be possible, so
            we can assume that this is the case when qq is a power
            of a ladder operator.
            """
            out.append(qq)
        
        elif isinstance(qq, KroneckerDelta):
            out.append(1 if (qq.args[0]==q.args[1]) 
                       else 0)
                
        elif isinstance(qq, Mul):
            if not qq.has(KroneckerDelta):
                out.append(qq)
            else:
                _out = []
                for arg in qq.args:
                    if isinstance(arg, KroneckerDelta):
                        if arg.args[0]==arg.args[1]:
                            continue
                        else:
                            _out = [Number(0)]
                            break
                    _out.append(arg)
                out.append(Mul(*_out))
        else: 
            raise InvalidTypeError([KroneckerDelta, Pow, Mul, Add], 
                                    type(qq))
            
    return Add(*out)