from ast import Expr

from sympy import sympify

from pybolano.core.normal_ordering import normal_ordering

############################################################

__all__ = ["NO_commutator", "expand_comm_AB_C", "expand_comm_A_BC", "expand_comm_AB_CD"]

############################################################


def _break_comm(A: Expr, B: Expr) -> Expr:
    """
    To avoid SymPy evaluating the commutatotr to a Kronecker
    delta.
    """
    
    return sympify(A * B - B * A).expand()

############################################################


def NO_commutator(A: Expr, B: Expr) -> Expr:
    """
    Return the normal-ordered equivalent of the commutator
    of two arbitrary polynomials of bosonic ladder operators.

    Parameters
    ----------

    A : sympy.Expr
        Operator in the left-hand slot of the commutator bracket.

    B : sympy.Expr
        Operator in the right-hand slot of the commutator bracket.

    Returns
    -------

    out : sympy.Expr
        Normal-oredered commutator between A and B.
    """

    return normal_ordering(_break_comm(A, B))
