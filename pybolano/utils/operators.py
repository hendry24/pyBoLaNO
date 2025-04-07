import random
from typing import TypeGuard

from sympy import Expr, Mul, Add, Number, Pow, Basic, Symbol, latex, conjugate, sympify
from sympy.physics.quantum import Dagger
from sympy.physics.quantum.boson import BosonOp

from pybolano.utils.error_handling import InvalidTypeError

############################################################

__all__ = [
    "BosonicAnnihilationOp",
    "BosonicCreationOp",
    "ops",
    "is_ladder",
    "is_ladder_contained",
    "dagger",
    "get_ladder_attr",
    "separate_mul_by_sub",
    "random_ladder",
]

############################################################

def _handle_sub(s):
    if s is None:
        s = ""
    if not(isinstance(s, Expr)):
        s = Symbol(str(s))
    return s

class pybolanoOp(Expr):
    
    symb = NotImplemented
    is_commutative = False
    
    def __new__(cls, sub = None):
        if sub is None:
            sub = ""
        if not(isinstance(sub, Basic)):
            sub = Symbol(str(sub))
        return Basic.__new__(cls, sub)
    
    @property
    def sub(self):
        return self.args[0]
    
    def __str__(self):
        return r"{%s_{%s}}" % (self.symb, latex(self.sub))
    
    def __repr__(self):
        return str(self)
    
    def _latex(self, printer):
        return str(self)
    
class BosonicAnnihilationOp(pybolanoOp):
    symb = r"b"
    
class BosonicCreationOp(pybolanoOp):
    symb = r"b^{\dagger}"

###

def ops(sub: str | Symbol | None = None) -> tuple[BosonicAnnihilationOp, BosonicCreationOp]:
    """
    Bosonic ladder operator objects.

    Parameters
    ----------

    sub : scalar or `sympy.Number` or `sympy.Symbol`, default: None
        Subscript of the boson ladder objects, used
        do differentiate the ladder operators for
        different subsystems in a multipartite system.
        Anything other than a sympy.Symbol is converted to
        one using `sympy.Symbol(latex(k))`. If this fails, an
        error is raised.

    Returns
    -------

    b : BosonicAnnihilationOp
        Boson annihilation object.

    bd : BosonicCreationOp
        Boson creation object.

    """
    
    return BosonicAnnihilationOp(sub=sub), BosonicCreationOp(sub=sub)

############################################################

def is_ladder(q: object) -> TypeGuard[pybolanoOp]:
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
    return isinstance(q, pybolanoOp)

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
    return q.has(pybolanoOp)

############################################################

def dagger(q : Expr) -> Expr:

    def Op_dag(op : pybolanoOp) -> pybolanoOp:
        if isinstance(op, BosonicAnnihilationOp):
            return BosonicCreationOp(sub = op.sub)
        elif isinstance(op, BosonicCreationOp):
            return BosonicAnnihilationOp(sub = op.sub)
    
    q = sympify(q).expand()
    if not(isinstance(q, Add)):
        q_args  = [q]
    else:
        q_args = q.args
    
    q_dag = 0        
    for qq in q_args:
        qq_dag = 1
        if not(is_ladder_contained(qq)):
            qq_dag *= conjugate(qq)
        elif is_ladder(qq):
            qq_dag *= Op_dag(qq)
        elif isinstance(qq, Pow):
            qq_dag *= Op_dag(qq.args[0]) ** qq.args[1]
        elif isinstance(qq, Mul):
            qq_dag = 1
            for qqq in reversed(qq.args):
                qq_dag *= dagger(qqq)
        else:
            raise ValueError("Input is not a polynomial in the ladder operators.")
        q_dag += qq_dag
    
    return q_dag

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
        sub = q.sub
        exp = Number(1)
    elif isinstance(q, Pow):
        if is_ladder_contained(q):
            sub = q.args[0].sub
            exp = q.args[1]
        else:
            raise ValueError("q is Pow but does not contain pybolanoOp.")
    else:
        raise InvalidTypeError([pybolanoOp, Pow], type(q))

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
    
    if isinstance(q, (Number, Symbol, Pow, pybolanoOp)):
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
            [Number, Symbol, pybolanoOp, Pow, Mul], type(q)
        )

############################################################

def random_ladder(n_ladder: int, 
                  k: str | Expr | list[str, Expr] = ""
                  ) -> Expr:
    """
    Generate a random ladder operator monomial containing a specified
    number of ladder operators whose subscripts are randomly chosen from a list
    of input subscripts.
    
    Parameters
    ----------
    
    n_ladder : int
        Total number of ladder operators in the monomial.
        
    k : str or sympy.Expr or list
        List of subscripts from which the function randomly samples and assign to
        a given ladder operator in the monomial.
    """
    if not(isinstance(k, list)):
        k = [k]
    return Mul(*[ops(random.choice(k))[random.randrange(2)] for _ in range(n_ladder)])

############################################################

def to_sympy_physics(q : Expr) -> Expr:
    """
    Convert expressions such that this package's ladder operators is replaced by the
    corresponding `BosonOp` objects in `sympy.physics.quantum`. 
    
    Parameters
    ----------
    
    q : Expr
        Expression to convert. Does not work with expectation value objects, which are
        visual-only. 
        
    Returns
    -------
    
    out : Expr
        Equivalent expression compatible with `sympy.physics.quantum`.
    """
    def _treat_factor(q):
        if not(is_ladder_contained(q)):
            return q
        else:
            k, pow = get_ladder_attr(q)
            b_sympy = BosonOp(r"b_{%s}" % k)
            if q.has(BosonicAnnihilationOp):
                return b_sympy**pow
            else:
                return Dagger(b_sympy)**pow

    q = sympify(q).expand()
    
    if not(is_ladder_contained(q)):
        return q
    elif isinstance(q, Add):
        return Add(*[to_sympy_physics(qq) for qq in q.args])
    elif isinstance(q, Mul):
        return Mul(*[_treat_factor(qq) for qq in q.args])
    else:
        return _treat_factor(q)