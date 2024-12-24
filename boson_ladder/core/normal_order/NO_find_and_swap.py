from sympy import (
    Pow,
    Mul,
    Add,
)
from sympy.physics.secondquant import (
    CreateBoson,
    AnnihilateBoson
)
from ..commutator.do_commutator import (
    _do_commutator_b_p_bd_q
)
from ...utils.operators import (
    is_ladder,
    is_ladder_contained
)

############################################################

__all__ = []

############################################################

def _NO_find_and_swap(q):
    """
    Input is a Mul object containing ladder operators with 
    the same subscript. The input to this function should be
    the outputs of _separate_mul_by_sub, which should contain
    no scalars. 
    """

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

############################################################
