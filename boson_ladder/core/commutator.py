from sympy import \
    Number, \
    Mul, \
    Pow, \
    Add, \
    factorial
from sympy.physics.secondquant import \
    CreateBoson, \
    AnnihilateBoson, \
    Commutator
from ..utils.operators import \
    is_ladder, \
    is_ladder_contained, \
    get_ladder_attr
from ..utils.commutators import \
    _isolate_bracket, \
    _treat_Kron
from ..utils.error_handling import \
    InvalidTypeError
from .normal_order import \
    normal_ordering

__all__ = ["do_commutator",
           "expand_AB_C",
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

def _do_commutator_b_p_bd_q(b_p, bd_q):
    """
    [b**p, bd**q] where b is an AnnihateBoson
    object and bd is a CreateBoson object.
    """

    def _comb(a, b):
        return factorial(a) / (factorial(b)*factorial(a-b))
    
    ###
    
    # Shortcuts
    
    sub_b_p, exp_b_p = get_ladder_attr(b_p)
    sub_bd_q, exp_bd_q = get_ladder_attr(bd_q)
    
    if not(is_ladder_contained(b_p)) or \
        not(is_ladder_contained(bd_q)) or \
        sub_b_p != sub_bd_q:
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
    

def _eval_sole_comm(comm):
    """
    Get A, B, C for the expansion and expand accordingly. 
    
    The eligible entry contains either a power or a multiplication
    between different powers of the two operators. We can find the
    leading contribution by checking the first three letters of
    the srepr, then accordingly separate the quantity with respect 
    to which we expand the commutator.
    
    Find the first occurence of operator in the bracket with
    # respect to which the commutator is expanded.
    
    The commutator to this input may be flipped from what the user
    inputs, but there is no problem. This function only cares about
    what its input is and will do its job accordingly, i.e. to expand
    with respect to the leftmost eligible operator. 
    """
    
    if not(isinstance(comm, Commutator)):
        raise InvalidTypeError(Commutator, type(comm))
    
    comm_1, comm_2 = comm.args
    
    """
    == ZERO CASES ==
    
    We check if the commutator is zero to speed up the evaluation. 
    Some of these cases should never occur since SymPy automatically 
    evaluates simple commutators.
    
    NOTE: The order in which the conditions are written are optimized.
        
    >   Conditions 1 and 2 applies regardless of the subscripts. If there
        are only either one of the two operators, then the commutator is 
        zero.
        
    NOTE: May want to add a condition to check if, in the case that conditions
        1 and 2 are false, the subscripts are such that the commutator is zero.
        This is not implemented since it can potentially be more computationally
        expensive than letting going through the recursions. 
      
    >   Condition 3 and 4 are obvious. If at least one slot does not contain
        any ladder operators, then the commutator is zero. This gets put last
        since SymPy should automatically evaluate commutators in this form,
        which means that it never gets to this function.

    """
        
    b_in_comm = comm_1.has(AnnihilateBoson) or comm_2.has(AnnihilateBoson)
    bd_in_comm = comm_1.has(CreateBoson) or comm_2.has(CreateBoson)
    
    if (b_in_comm and not(bd_in_comm)) \
        \
        or (not(b_in_comm) and bd_in_comm) \
        \
        or not(is_ladder_contained(comm_1)) \
        \
        or not(is_ladder_contained(comm_2)):
    
        return Number(0)
    
    ##########
    
    if not(is_ladder(comm_1)):
        # skip the check if the entry is a single operator
        if isinstance(comm_1, Mul):
            A = comm_1.args[0]
            B = Mul(*comm_1.args[1:])
        elif isinstance(comm_1, Pow):
            if not(comm_1.args[1].is_number):
                msg = "The exponent of the ladder operator "
                msg += "must be a Number."
                raise ValueError(msg)
            A = comm_1.args[0]
            B = Pow(A, comm_1.args[1]-1)
        else:
            raise InvalidTypeError([Mul, Pow], type(comm_1))
        C = comm_2
        
        return expand_AB_C(A, B, C)
        
    elif not(is_ladder(comm_2)):
        if isinstance(comm_2, Mul):
            B = comm_2.args[0]
            C = Mul(*comm_2.args[1:])
        elif isinstance(comm_2, Pow):
            if not(comm_2.args[1].is_number):
                msg = "The exponent of the ladder operator "
                msg += "must be a Number."
                raise ValueError(msg)
            B = comm_2.args[0]
            C = Pow(B, comm_2.args[1]-1)
        else:
            raise InvalidTypeError([Mul, Pow], type(comm_2))
        A = comm_1
        
        return expand_A_BC(A, B, C)
        
    else:   
        return comm
        # though both entries cannot both be ladder operators since SymPy
        # would evaluate comm upon construction. 
    
def _expand_addend(q):
    """
    Utilizing the properties
        [AB,C] = A[B,C] + [A,C]B
        [A,BC] = [A,B]C + B[A,C]
    expand the commutator into a sum of simpler commutators.
    Each call does this once, starting with the leftmost 
    eligible operator, e.g. A in [AB,C] and B in [A,BC].
    
    The expanded commutators are returned as a list of single
    commutators for the recursion used by the main function.
    
    NOTE: The output not being a list when there are no brackets left 
    is necessary to quit the recursion in the main function.
    """

    # We start with a commutator in the bracket form, if
    # not directly solved by SymPy. This means that comm
    # must be input without calling .doit first.
    
    if not(q.has(Commutator)):
        return _treat_Kron(q), True
                                # stop flag.
    
    left_factor, comm, right_factor = _isolate_bracket(q)
    # At this point, comm should purely be a Commutator object.
                
    comm = _eval_sole_comm(comm)
    # Should generally return an Add. 
    comm = (left_factor*comm*right_factor).expand()
                                            # lays it flat.
    
    if isinstance(comm, Add):
        return [arg for arg in comm.args], False
    elif not(comm.has(Commutator)):
        return _treat_Kron(comm), True
    else:
        return [comm], False
        
###################################################
                
def do_commutator(A, B, normal_order = True):
    """
    Calculate the commutator of two arbitrary
    polynomials of bosonic ladder operators.
    
    Parameters
    ----------
    
    A : sympy.Expr
        Operator in the left-hand slot of the commutator bracket.
        
    B : sympy.Expr
        Operator in the right-hand slot of the commutator bracket.
        
    normal_order : bool, default: True
        Whether to normal-order the result. If `True`, `normal_ordering`
        is called before the result is returned.
        
    Returns
    ------- 
    
    out : sympy.Expr
        Commutator between A and B, optionally normal ordered.
    """
    
    comm = Commutator(A, B)
    """
    When the Commutator object is initialized, any sum in the
    input will automatically result in a sum of Commutators.
    """
    
    if not(comm.has(Commutator)):
        return _treat_Kron(comm)
    
    if isinstance(comm, Add):
        recursion_input = comm.args
    else:
        recursion_input = [comm]
    
    single_comms = []
    def expand_recursive(comm):
        res, stop_flag = _expand_addend(comm)
        if stop_flag:
            single_comms.append(res)
        else:
            for item in res:
                expand_recursive(item)
                
    for item in recursion_input:
        expand_recursive(item)
        
    out = Add(*single_comms)
    
    if normal_order:
        out = normal_ordering(out)
    
    return out