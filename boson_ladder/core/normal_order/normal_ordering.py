from sympy import (
    Mul,
    Add,
    Pow,
    Number
)
from sympy.physics.secondquant import (
    AnnihilateBoson,
    CreateBoson,
)
from multiprocessing import (
    Pool
)
from functools import (
    partial
)
from .NO_Blasiak import (
    _NO_Blasiak
)
from .NO_find_and_swap import (
    _NO_find_and_swap
)
from ...utils.operators import (
    separate_mul_by_sub,
    is_ladder_contained
)
from ...utils.multiprocessing import (
    mp_config
)
from ...utils.error_handling import (
    InvalidTypeError
)

############################################################

__all__ = ["normal_ordering"]    

############################################################

def _NO_preprocess(q):
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
    return q_args
    
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

def _NO_input_addend(qq, NO_single_sub):
    qq_Mul_args = separate_mul_by_sub(qq)
    
    out = Number(1) if is_ladder_contained(qq_Mul_args[0]) \
            else qq_Mul_args.pop(0)
    
    for qq_single_sub in qq_Mul_args:
        out *= NO_single_sub(qq_single_sub)
    
    return out.expand()

############################################################

def normal_ordering(q, method="Blasiak"):
    """
    Normal order the operator q: all creation operators are
    written to the left of all annihilation operators within
    a single term. 
    
    Parameters
    ----------
    
    q : sympy.Expr
        The operator to normal order.
        
    method : str, default: "Blasiak"
        Normal ordering method.
        
        (1) "find-and-swap" or "fns" or 1
        The conventional find-and-swap method. The algorithm
        recursively looks for `(b**p*bd**q)` sequences and 
        expands them using the commutation relations. Avoid
        this method for complex expressions as the recursion
        tree can be very large.
        
        (2) "Blasiak" or "2"
        Explicit formula derived by Blasiak (see 
        https://arxiv.org/abs/quant-ph/0507206. This is the
        faster method since it does not involve recursions.

    Returns
    -------
    
    q_NO : sympy.Expr
        q, normal-ordered.
        
    """
    
    # Shortcuts
    
    if not(is_ladder_contained(q)) \
        or isinstance(q, (Pow,
                          CreateBoson,
                          AnnihilateBoson)):
        return q
        
    # NOTE: checking if all bd are already to the
    # left of b may be too computationally expensive
    # to implement here.
    
    ###
    
    match method.lower():
        case "find-and-swap" | "fns" | 1:
            NO_single_sub = _NO_find_and_swap
        case default:
            NO_single_sub = _NO_Blasiak
    ###
    
    q_args = _NO_preprocess(q)
    use_mp = (mp_config["enable"] 
                and (len(q_args) >= mp_config["min_num_args"]))
    
    ###
    
    if use_mp:
        with Pool(mp_config["num_cpus"]) as pool:
            _out = Add(*pool.map(partial(_NO_input_addend, 
                                         NO_single_sub=NO_single_sub), 
                                 q_args))
    else:
        _out = Add(*[_NO_input_addend(qq, NO_single_sub) 
                     for qq in q_args])
    
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
            return Add(*pool.map(_final_swap, 
                                 _out.args))
    else:
        return Add(*[_final_swap(q) for q in _out.args])
