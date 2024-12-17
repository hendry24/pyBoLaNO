from sympy import \
    Add, \
    Pow, \
    Mul
from sympy.physics.secondquant import \
    AnnihilateBoson, \
    CreateBoson
from ..utils.operators import \
    _flatten_pow
from ..utils.error_handling import \
    InvalidTypeError

__all__ = ["normal_ordering"]

def _NO_one_step(q_args):
    """
    Normal order the first [b,bd] sequence found in q_args 
    and return the result as a list of two args.
    
    q_args = q.args. Here q is a single operator (Mul or Pow)
    to be normal ordered, and q_args must be made free of Pow
    by laying out the power in terms of multiplication instead.
    `get_args_no_pow` does this.
    
    Leading scalar does not matter.
    """
    
    for i, op in enumerate(q_args):
        
        # Find bd to the right of a b.
        if isinstance(op, CreateBoson):
            if i == 0 \
                or not(isinstance(q_args[i-1], AnnihilateBoson)):
                continue
        ###
        
            b_sub = q_args[i-1].args[0]
            bd_sub = op.args[0]
            args_left = q_args[:i-1] # to the left of the substituted b*bd.
            args_right = q_args[i+1:] # to the right of the substituted b*bd.
            
            b = AnnihilateBoson(b_sub)
            bd = CreateBoson(bd_sub)
            
            try:
                q_NO_mul_args = [args_left + [bd,b] + args_right]
            except:
                """
                In some cases (probably only the first recursion stack), args_left
                and args_right may not be lists, so tuples are needed. Anyway, it 
                seems to only happen when we normal order b*bd, so this may only
                rarely happen.
                """
                q_NO_mul_args = [args_left + (bd,b) + args_right]
            
            if b_sub == bd_sub:
                q_NO_mul_args.append(args_left + args_right)
                            
            """
            We return the quantities as arguments as they will be used
            for the recursion, so there is no need to call _flatten_pow
            repeatedly.
            """
                
            return q_NO_mul_args, False
        
    # If nothing is found, it means the NO is finished and we can return
    # a True stop flag. We use a boolean flag since both possible outputs
    # are lists of different
    return q_args, True
        
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
    
    if isinstance(q, (Pow, CreateBoson, AnnihilateBoson))\
        or (not(q.has(CreateBoson))
            and not(q.has(AnnihilateBoson))):
        return q
    
    elif isinstance(q, Add):
        q_args = [_flatten_pow(arg) 
                  for arg in q.args]
    
    elif isinstance(q, Mul):
        q_args = [_flatten_pow(q)]
    
    else:
        raise InvalidTypeError([Pow,
                                CreateBoson,
                                AnnihilateBoson,
                                Add,
                                Mul],
                                type(q))
    
    # The original operator is Add(*[Mul(*x) for x in q_args])
    # since get_args_no_pow returns a list of factors.
    
    ###
    
    q_args_NO = []
    
    def _NO_recursive(arg):
        """
        Recursively get the arguments of the sum of the
        normal-ordered operators and put them in q_args_NO
        such that the output of the main function is 
        Add(*q_args_NO).
        """
        
        res, stop_flag = _NO_one_step(arg)
            
        if stop_flag:
            """
            If the stop flag is True, then out is fully
            normal ordered. NO_one_step does
            not return a list of addends of the two
            terms arising from the substitution, but
            a list of factors that is input. We can then
            make the Mul object and append it to q_args_NO. 
            """  
            q_args_NO.append(Mul(*res))
        else:
            for arg in res:
                """
                For each addend in the resulting list of NO_one_step
                containing the two terms arising from the substitution,
                we do another step of normal ordering. arg here is a 
                list of the factors of each addend that is free of any
                Pow since it is the output of NO_one_step. 
                """
                _NO_recursive(arg)
    
    for arg in q_args:
        _NO_recursive(arg)
        
    return Add(*q_args_NO)
