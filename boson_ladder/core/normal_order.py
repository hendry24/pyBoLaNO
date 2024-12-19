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
    is_ladder_contained
from .do_commutator import \
    _do_commutator_b_p_bd_q
from ..utils.error_handling import \
    InvalidTypeError

__all__ = ["normal_ordering"]

def _NO_single_sub(q):
    """
    Input is a Mul object containing ladder operators with 
    the same subscript. The input to this function should be
    the outputs of _separate_mul_by_sub, which should contain
    no scalars. If a scalar is found, the input is returned as
    is.
    """

    if q.has((Number, Symbol)):
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
    
    # return Add(
    #         *[Mul(
    #             *[_NO_single_sub(qq_single_sub)
    #                    for qq_single_sub in _separate_mul_by_sub(qq)
    #                    ]
    #             ).expand()
    #             for qq in q_args]
    #         )
    
    # Equivalent to,
    #
    # q_NO_Add_args = []
    # for qq in q_args:
    #     qq_NO_Mul_args = []
    #     for qq_single_sub in _separate_mul_by_sub(qq):
    #         qq_NO_Mul_args.append(_NO_single_sub(qq_single_sub))
    #     q_NO_Add_args.append(Mul(*qq_NO_Mul_args).expand())
    # return Add(*q_NO_Add_args)