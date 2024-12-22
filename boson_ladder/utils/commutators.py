from sympy import \
    Add, \
    Mul, \
    Pow, \
    Number, \
    factorial, \
    KroneckerDelta
from sympy.physics.secondquant import \
    CreateBoson, \
    AnnihilateBoson, \
    Commutator
from .operators import \
    is_ladder_contained, \
    get_ladder_attr
from .error_handling import InvalidTypeError
    
__all__ = []

def _do_commutator_b_p_bd_q(b_p, bd_q):
    """
    [b**p, bd**q] where b is an AnnihateBoson
    object and bd is a CreateBoson object.
    """

    def _comb(a, b):
        return factorial(a) / (factorial(b)*factorial(a-b))
    
    ###
    
    # Shortcuts
    
    if b_p.has(CreateBoson):
        raise InvalidTypeError(AnnihilateBoson, type(b_p))
    if bd_q.has(AnnihilateBoson):
        raise InvalidTypeError(CreateBoson, type(bd_q))
    
    if not(b_p.has(AnnihilateBoson)) \
        or not(bd_q.has(CreateBoson)):
        return Number(0)
    
    sub_b_p, exp_b_p = get_ladder_attr(b_p)
    sub_bd_q, exp_bd_q = get_ladder_attr(bd_q)
    
    if sub_b_p != sub_bd_q:
        return Number(0)
    
    ### 
    
    if isinstance(b_p, AnnihilateBoson):
        b = b_p
        p = 1
    elif isinstance(b_p, Pow) \
        and b_p.has(AnnihilateBoson):
        b, p = b_p.args
    else:
        raise InvalidTypeError([AnnihilateBoson, Pow],
                               type(b_p))
    
    if isinstance(bd_q, CreateBoson):
        bd = bd_q
        q = 1
    elif isinstance(bd_q, Pow) \
        and bd_q.has(CreateBoson):
        bd, q = bd_q.args
    else:
        raise InvalidTypeError([CreateBoson, Pow],
                               type(bd_q))
    
    out = Number(0)
    for k in range(max(0, p-q), p):
        out += _comb(p, k) * _comb(q, p-k) * factorial(p-k)\
                * bd**(q-p+k) * b**k 
    return out

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