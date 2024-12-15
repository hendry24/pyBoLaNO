from sympy import srepr, Mul, Add, Pow
from sympy.physics.secondquant import CreateBoson, AnnihilateBoson

from .utils import *

__all__ = ["normal_order"]

b, bd = ops()
assert isinstance(b, AnnihilateBoson) or isinstance(bd, CreateBoson)

def get_args_no_pow(q):
    """
    Layout any power expressions in q.args and 
    return a list of arguments without power.
    """
    
    if isinstance(q, Pow):
        return (q.args[0] for _ in range(q.args[1]))
    
    assert isinstance(q, Mul)
    
    if "Pow" in srepr(q):
        args_long = []
        for arg in q.args:
            if isinstance(arg, Pow):
                for _ in range(arg.args[1]):
                    args_long.append(arg.args[0])
            else:
                args_long.append(arg)
    else:
        args_long = q.args
        
    return args_long

def NO_one_step(q_args):
    """
    Normal order the first [b,bd] sequence found in q_args 
    and return the result as a list of two args.
    
    q_args = q.args. Here q is a single operator (Mul or Pow)
    to be normal ordered, and q_args must be made free of Pow
    by laying out the power in terms of multiplication instead.
    `get_args_no_pow` does this.
    
    Leading scalar does not matter.
    """
    
    find_bd = False 
    # We skip the algorithm for bd that are already
    # leftmost. The algorithm starts looking for bd
    # to swap after it meets the first b.
    for i, op in enumerate(q_args):
        if isinstance(op, AnnihilateBoson):
            find_bd = True
            continue
        
        if not(find_bd):
            continue
        
        # The first bd found must be right next
        # to a `b`.
        if isinstance(op, CreateBoson):
            args_left = q_args[:i-1] # To the left of the substituted b*bd.
            args_right = q_args[i+1:]
            
            try:           
                q_NO_mul_args = [args_left + args_right,   
                                    # Integer(1) is not necessary, since Mul(*()) = Integer(1)
                                args_left + [bd,b] + args_right]
                # where q_NO = Add(Mul(*q_NO_mul_args[0]), 
                #                  Mul(*q_NO_mul_args[1]))
                """
                We return the quantities as arguments as they will be used
                for the recursion, so there is no need to call get_args_no_pow
                repeatedly.
                """
            except:
                """
                In some cases (probably only the first recursion stack), args_left
                and args_right may not be lists, so tuples are needed. Anyway, it 
                seems to only happen when we normal order b*bd, so this may only
                rarely happen.
                """
                q_NO_mul_args = [args_left + args_right,   
                                args_left + (bd,b) + args_right]
                
            return q_NO_mul_args, False
        
    # If nothing is found, it means the NO is finished and we can return
    # a True stop flag. We use a boolean flag since both possible outputs
    # are lists of different depths.
    return q_args, True
    
def normal_order(q):
    """
    Normal order the operator q.
    
    Parameters
    ----------
    
    q : sympy.Mul or sympy.Add containing sympy.physics.secondquant.AnnihilateBoson and sympy.physics.secondquant.CreateBoson
        The operator to normal order.
        
    Returns
    -------
    
    q_NO : sympy.Mul or sympy.Add containing sympy.physics.secondquant.AnnihilateBoson and sympy.physics.secondquant.CreateBoson
        q, normal ordered.
    """
    if isinstance(q, (Pow, CreateBoson, AnnihilateBoson))\
        or (not("CreateBoson" in srepr(q))
            and not("AnnihilateBoson" in srepr(q))):
        return q
    elif isinstance(q, Add):
        q_args = [get_args_no_pow(arg) for arg in q.args]
    elif isinstance(q, Mul):
        q_args = [get_args_no_pow(q)]
    
    # The original operator is Add(*[Mul(*x) for x in q_args])
    # since get_args_no_pow returns a list of factors.
    
    ###
    
    q_args_NO = []
    
    def NO_recursive(arg):
        """
        Recursively get the arguments of the sum of the
        normal-ordered operators and put them in q_args_NO
        such that the output of the main function is 
        Add(*q_args_NO).
        """
        
        out, stop_flag = NO_one_step(arg)
            
        if stop_flag:
            """
            If the stop flag is True, then out is fully
            normal ordered. NO_one_step does
            not return a list of addends of the two
            terms arising from the substitution, but
            a list of factors that is input. We can then
            make the Mul object and append it to q_args_NO. 
            """  
            q_args_NO.append(Mul(*out))
        else:
            for arg in out:
                """
                For each addend in the resulting list of NO_one_step
                containing the two terms arising from the substitution,
                we do another step of normal ordering. arg here is a 
                list of the factors of each addend that is free of any
                Pow since it is the output of NO_one_step. 
                """
                NO_recursive(arg)
    
    for arg in q_args:
        NO_recursive(arg)
        
    return Add(*q_args_NO)
