from sympy import \
    srepr, \
    Integer, \
    Mul, \
    Pow, \
    Add
from sympy.physics.secondquant import \
    Commutator
from ..utils.operators import \
    is_ladder
from ..utils.commutators import \
    expand_A_BC, \
    expand_AB_C, \
    _isolate_bracket, \
    _treat_Kron
from ..utils.normal_ordering import \
    normal_order as NO

__all__ = ["do_commutator"]

def get_ABC_and_expand(comm):
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
    
    def comm_entry_err_msg(entry):
        msg = "Commutator entry is neither Mul or Pow. Here is the srepr: \n"
        msg += srepr(entry)
        raise ValueError(msg)
    
    comm_1, comm_2 = comm.args
    
    # First, we need to check if the commutator is zero, since 
    # the following algorithm does not take into account cases 
    # like [b, b**9], which is treated as an expandable expression
    # due to the Pow in its srepr. This leads to faulty evaluation.

    concatenated_srepr = srepr(comm_1) + srepr(comm_2)
    b_in_comm = "AnnihilateBoson" in concatenated_srepr
    bd_in_comm = "CreateBoson" in concatenated_srepr
    if (b_in_comm and not(bd_in_comm)) or (not(b_in_comm) and bd_in_comm):
        return Integer(0)
    # NOTE: Might want to change this for many body cases.
    
    if not(is_ladder(comm_1)):
        # skip the check if the entry is a single operator
        operation = srepr(comm_1)[:3]
        if operation == "Mul":
            A = comm_1.args[0]
            B = Mul(*comm_1.args[1:])
        elif operation == "Pow":
            A = comm_1.args[0]
            B = Pow(A, comm_1.args[1]-1)
        else:
            comm_entry_err_msg(comm_1)
        C = comm_2
        
        return expand_AB_C(A, B, C)
        
    elif not(is_ladder(comm_2)):
        operation = srepr(comm_2)[:3]
        if operation == "Mul":
            B = comm_2.args[0]
            C = Mul(*comm_2.args[1:])
        elif operation == "Pow":
            B = comm_2.args[0]
            C = Pow(B, comm_2.args[1]-1)
        else:
            comm_entry_err_msg(comm_2)
        A = comm_1
        
        return expand_A_BC(A, B, C)
        
    else:   
        return comm
        # though both entries cannot both be ladder operators since SymPy
        # would evaluate comm upon construction. 
    
def expand_single_comm(comm):
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

    if "Commutator" not in srepr(comm):
        return _treat_Kron(comm) 
                # in case a Symbol is used as the input to the ladder operators.
    
    if "Add" in srepr(comm):
        """
        When Commutator is initialized, any sum in the
        input will automatically result in a sum of Commutators.
        This ensures that the function only deals with one bracket
        at a time. 
        """
        msg = "Input accepts only one commutator term, "
        msg += "not a sum of them"
        raise ValueError(msg)
    
    left_factor, comm, right_factor = _isolate_bracket(comm)
    # At this point, comm should purely be a commutator object.
                
    comm = get_ABC_and_expand(comm)
    comm = (left_factor*comm*right_factor).expand()
            # Should generally return a sum of commutators.
    
    if not("Commutator" in srepr(comm)) \
        or len(comm.args) == 1:
        return _treat_Kron(comm)
    else:
        return [term for term in comm.args]
                
def do_commutator(A, B, normal_order = True):
    comm = Commutator(A, B)
    
    single_comms = []
    
    def expand_recursive(comm):
        expanded_comm = expand_single_comm(comm)
        if isinstance(expanded_comm, list):
            for item in expanded_comm:
                expand_recursive(item)
        else:
            single_comms.append(expanded_comm)

    expand_recursive(comm)
    
    out = Add(*single_comms)
    
    if normal_order:
        out = NO(out)
    
    return out