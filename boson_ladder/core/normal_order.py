from multiprocessing import \
    Pool
from sympy import \
    Add, \
    Pow, \
    Mul, \
    Number, \
    Symbol
from sympy.physics.secondquant import \
    AnnihilateBoson, \
    CreateBoson
from ..utils.operators import \
    is_ladder, \
    is_ladder_contained, \
    get_ladder_attr, \
    separate_mul_by_sub
from ..utils.commutators import \
    _do_commutator_b_p_bd_q
from ..utils.error_handling import \
    InvalidTypeError
from ..utils.multiprocessing import \
    mp_config

__all__ = ["normal_ordering"]

def _NO_single_sub(q):
    """
    Input is a Mul object containing ladder operators with 
    the same subscript. The input to this function should be
    the outputs of _separate_mul_by_sub, which should contain
    no scalars. 
    """
    
    if not(is_ladder_contained(q)):
        return q
    
    if isinstance(q, (Pow, 
                      CreateBoson, 
                      AnnihilateBoson)):
        return q

    def _NO_single_Add_term(qq):
        """
        In each step, iterate through the arguments qq of the
        output Add. Find all (b**p, bd**q) sequence and commute them. 
        Let the commutation output be an Add object which we will
        multiply together and return as the output.
        
        e.g. We want (b**2*bd**2 * b**2*bd**3) to become
        (2+4*bd*b+bd**2*b**2)(6*bd+6*bd**2*b*bd**3*b**2) 
        whose expansion is returned.
        
        stop_flag tells the recursion to stop when the normal
        ordering is done.
        """
        
        if not(is_ladder_contained(qq)) or \
            is_ladder(qq):
            return qq
                        
        out_Mul_args = []
        i = 0
        # qq must be Mul at this point
        while True:
            if i > (len(qq.args)-1):
                break
             
            qqq = qq.args[i]
            if qqq.has(AnnihilateBoson):
                if i == (len(qq.args)-1):
                    out_Mul_args.append(qqq)
                    break
                
                qqq_next = qq.args[i+1]
                # If the expression is not normal ordered yet,
                # next to a Pow of b we must have a Pow of bd since
                # there is only one subscript.
                
                out_Mul_args.append(_do_commutator_b_p_bd_q(qqq, 
                                                           qqq_next)
                                    + qqq_next*qqq)
                i += 2
            else:
                out_Mul_args.append(qqq)
                i += 1
            
        return Mul(*out_Mul_args).expand()

    out_Add_args = []
    def _recursion(qq):
        res = _NO_single_Add_term(qq)
        if isinstance(res, Add):
            for item in res.args:
                _recursion(item)
        else:
            out_Add_args.append(res)
            
    _recursion(q)
    
    return Add(*out_Add_args)

def _NO_task_per_q_arg(qq):
    _out_Mul_args = []
    for qq_single_sub in separate_mul_by_sub(qq):
        # NOTE: Not sure if its worthy to multiprocess here
        # except when the number of subscripts are large.
        _out_Mul_args.append(_NO_single_sub(qq_single_sub))
    return Mul(*_out_Mul_args).expand()

def _final_swap(q):
    if not(isinstance(q, Mul)):
        return q
    else:
        collect_scalar = []
        collect_b = []
        collect_bd = []
        for qq in q.args:
            # factor
            if qq.has(AnnihilateBoson):
                collect_b.append(qq)
            elif qq.has(CreateBoson):
                collect_bd.append(qq)
            else:
                collect_scalar.append(qq)
        return Mul(*(collect_scalar+collect_bd+collect_b))
        
def normal_ordering(q):
    """
    Normal order the operator q: all creation operators are
    written to the left of all annihilation operators within
    a single term. 
    
    Parameters
    ----------
    
    q : sympy.Expr
        The operator to normal order.
        
    Returns
    -------
    
    q_NO : sympy.Expr
        q, normal-ordered.
    """
    
    q = q.expand()
    
    if (not(q.has(CreateBoson)) \
            and not(q.has(AnnihilateBoson))) \
        or isinstance(q, (Pow, 
                          CreateBoson, 
                          AnnihilateBoson)):
        return q
    
    elif isinstance(q, Add):
        q_args = [qq for qq in q.args]
                    
    elif isinstance(q, Mul):
        if not(is_ladder_contained(q)):
            return q
        q_args = [q]
    
    else:
        raise InvalidTypeError([Pow,
                                CreateBoson,
                                AnnihilateBoson,
                                Add,
                                Mul],
                                type(q))
        
    use_mp = mp_config["enable"] \
                and (len(q_args) >= mp_config["min_num_args"])
                    
    if use_mp:
        with Pool(mp_config["num_cpus"]) as pool:
            _out = Add(*pool.map(_NO_task_per_q_arg, q_args))
    else:
        _out = Add(*[_NO_task_per_q_arg(qq) for qq in q_args])
    
    """
    At this point, the normal ordering is not done since there are
    probably the terms are written something like 
        bd_0**p*b_0**q * bd_1**r*b_1**s
    The last step is simply to swap the argument order to get the
    nice-looking output with the same subscript order between
    the creation and the annihilation operators.
    """
    
    if not(isinstance(_out, Add)):
        return _out

    if use_mp:
        with Pool(mp_config["num_cpus"]) as pool:
            return Add(*pool.map(_final_swap, _out.args))
    else:
        return Add(*[_final_swap(q) for q in _out.args])