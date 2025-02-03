from random import randrange
from typing import TypeGuard, Union

from sympy import Expr, Mul, Number, Pow, Symbol, latex
from sympy.physics.secondquant import AnnihilateBoson, CreateBoson, SqOperator

from pybolano.utils.error_handling import InvalidTypeError

############################################################

__all__ = [
    "ops",
    "is_ladder",
    "is_ladder_contained",
    "get_ladder_attr",
    "separate_mul_by_sub",
    "random_ladder",
]

############################################################


def ops(k: str | Symbol | None = None) -> tuple[AnnihilateBoson, CreateBoson]:
    """
    SymPy's bosonic ladder operator objects.

    Parameters
    ----------

    k : scalar or `sympy.Number` or `sympy.Symbol`, default: None
        Subscript of the boson ladder objects, used
        do differentiate the ladder operators for
        different subsystems in a multipartite system.
        Anything other than a sympy.Symbol is converted to
        one using `sympy.Symbol(latex(k))`. If this fails, an
        error is raised.

    Returns
    -------

    b : AnnihilateBoson
        Boson annihilation object.

    bd : CreateBoson
        Boson creation object.

    """

    if k is None:
        k = Symbol("")
    elif isinstance(k, Symbol):
        pass
    else:
        try:
            k = Symbol(latex(k))
        except:
            raise ValueError("Invalid k.")

    b = AnnihilateBoson(k)
    bd = CreateBoson(k)
    return b, bd


############################################################

def is_ladder(q: object) -> TypeGuard[Union[AnnihilateBoson, CreateBoson]]:
    """
    Check if the input object is a ladder operator.

    Parameters
    ----------

    q : object
        Object to check.

    Returns
    -------

    out : bool
        `True` if q is a ladder operator. `False` otherwise.
    """
    return isinstance(q, (CreateBoson, AnnihilateBoson))


############################################################


def is_ladder_contained(q: Expr) -> bool:
    """
    Check if at least one ladder operator is contained within the
    input object.

    Parameters
    ----------

    q : sympy.Expr
        Object to check.

    Returns
    -------

    out : bool
        `True` if a ladder operator is contained. `False` otherwise.
    """
    return q.has(AnnihilateBoson, CreateBoson)


############################################################


def get_ladder_attr(q: Expr) -> tuple[Symbol, Number]:
    """
    Return the index and exponent of the ladder operator.

    Parameters
    ----------

    q : sympy.Expr
        Either a ladder operator or its power.

    Returns
    -------

    sub : sympy.Symbol or sympy.Number
        Subscript (index) of the ladder operator.

    exp : sympy.Number
        Exponent of the expression.
    """
    if is_ladder(q):
        sub = q.args[0]
        exp = Number(1)
    elif isinstance(q, Pow):
        sub = q.args[0].args[0]
        exp = q.args[1]
    else:
        raise InvalidTypeError([CreateBoson, AnnihilateBoson, Pow], type(q))

    return sub, exp


############################################################


def separate_mul_by_sub(q: Expr) -> list[Expr]:
    """
    Separate a Mul object by the ladder operator subscript.

    Parameters
    ----------

    q : sympy.Expr
        Input quantity.

    Returns
    -------

    out : list
        A list containing the arguments of `q` separated
        by subscript. Scalars are put in one group as the
        first entry.
    """
    if isinstance(q, (Number, Symbol, CreateBoson, AnnihilateBoson, Pow)):
        return [q]
    elif isinstance(q, Mul):
        out = {}
        for qq in q.args:
            if not (is_ladder_contained(qq)):
                if Number not in out:
                    out[Number] = []
                out[Number].append(qq)
            else:
                sub, exp = get_ladder_attr(qq)
                if sub not in out:
                    out[sub] = []
                out[sub].append(qq)
        return [Mul(*args) for args in list(out.values())]
    else:
        raise InvalidTypeError(
            [Number, Symbol, CreateBoson, AnnihilateBoson, Pow, Mul], type(q)
        )


############################################################

def random_ladder(n_ladder: int, 
                  k: str | Expr | list[str, Expr] = ""
                  ) -> Expr:
    if not(isinstance(k, list)):
        k = [k]
    return Mul(*[ops(k[randrange(len(k))])[randrange(2)] for _ in range(n_ladder)])