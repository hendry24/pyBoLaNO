from sympy import \
    Symbol, \
    Add, \
    Mul, \
    Pow, \
    Number, \
    latex
from sympy.physics.secondquant import \
    CreateBoson, \
    AnnihilateBoson
from .error_handling import \
    InvalidTypeError
    
__all__ = ["ops",
           "is_ladder",
           "is_ladder_contained",
           "get_ladder_attr",
           "separate_mul_by_sub"]

def ops(k=None):
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

def is_ladder(q):
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

def is_ladder_contained(q):
    """
    Check if at least one ladder operator is contained within the
    input object.
    
    Parameters
    ----------
    
    q : object
        Object to check.
        
    Returns
    -------

    out : bool
        `True` if a ladder operator is contained. `False` otherwise.
    """
    return q.has(AnnihilateBoson, CreateBoson)
            
def get_ladder_attr(q):
    """
    Return the index and exponent of the ladder
    operator.
    """
    if is_ladder(q):
        sub = q.args[0]
        exp =  Number(1)
    elif isinstance(q, Pow):
        sub = q.args[0].args[0]
        exp =  q.args[1]
    else:
        raise InvalidTypeError([CreateBoson,
                                AnnihilateBoson,
                                Pow],
                               type(q))
    
    return sub, exp

def separate_mul_by_sub(q):
    if isinstance(q, (Number,
                      Symbol,
                      CreateBoson,
                      AnnihilateBoson,
                      Pow)):
        return [q]
    elif isinstance(q, Mul):
        out = {}
        for qq in q.args:
            if not(is_ladder_contained(qq)):
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
        raise ValueError
