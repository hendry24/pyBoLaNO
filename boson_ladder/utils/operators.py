from sympy import \
    Symbol, \
    Add, \
    Mul, \
    Pow, \
    Number, \
    srepr, \
    latex
from sympy.physics.secondquant import \
    CreateBoson, \
    AnnihilateBoson
from .error_handling import \
    InvalidTypeError
    
__all__ = ["ops",
           "is_ladder",
           "is_ladder_contained",
           "separate_by_subscript"]

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
            
def _flatten_pow(q):
    """
    Layout any power expressions in q.args and 
    return a list of arguments without power, i.e.
    q = Mul(*output).
    
    Bad expressions such as a Power object with 
    negative power and non-Number object are
    kept as is.
    
    """
    
    if isinstance(q, Add):
        raise ValueError("q is not supposed to be Add.")
    
    if isinstance(q, Pow):
        if q.args[1] is not Number \
            or q.args[1] < 2: # bad power expressions that may raise error in the program.
            return [q]
        return [q.args[0] for _ in range(q.args[1])]
    
    if not(isinstance(q, Mul)):
        return [q]
    
    if "Pow" in srepr(q):
        args_long = []
        for arg in q.args:
            if isinstance(arg, Pow):
                if q.args[1] is not Number\
                    or q.args[1] < 2:
                    args_long.append(arg)
                else:
                    for _ in range(arg.args[1]):
                        args_long.append(arg.args[0])
            else:
                args_long.append(arg)
    else:
        args_long = q.args
        
    return args_long

def separate_by_subscript(q):
    """
    Separate an operator product according to the subscript.
    q is the output of flatten_pow. Returns a list, whose
    each element is a group of operators with the same
    index. Scalars are put into the first operator group.
    
    Parameters
    ----------
    
    q : sympy.Expr
        Operator product or a list of its arguments. 
        
    Returns
    -------
    
    out : list
        List of operator arguments separated by the subscript.
        Reconstruct the original operator with `sympy.Mul(*out)`.
    """
    
    if isinstance(q, (CreateBoson, AnnihilateBoson)):
        return q
    
    if isinstance(q, (Mul, Pow)):
        q = _flatten_pow(q)
    
    if not(isinstance(q, list)):
        raise InvalidTypeError([CreateBoson, 
                                AnnihilateBoson,
                                Mul,
                                Pow,
                                list],
                               q)
    
    sub_args = {}
    scalar = Number(1)
    for qq in q:
        if not(is_ladder(qq)):
            scalar *= qq
        
        sub = qq.args[0]
        if sub not in sub_args:
            sub_args[sub] = []
        
        sub_args[sub].append(qq)
        
    sub_args[list(sub_args.keys())[0]].append(scalar)
    
    return list(sub_args.values())
