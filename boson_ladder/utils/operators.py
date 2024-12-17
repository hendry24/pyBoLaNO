from sympy import \
    Symbol, \
    srepr
from sympy.physics.secondquant import \
    CreateBoson, \
    AnnihilateBoson
    
__all__ = ["ops",
           "is_ladder",
           "is_ladder_contained"]

def ops(k = 0):
    """
    SymPy's bosonic ladder operator objects.
    
    Parameters
    ----------
    
    k : scalar or `sympy.Integer` or `sympy.Symbol`, default: 0
        Subscript of the boson ladder objects, used 
        do differentiate the ladder operators for 
        different subsystems in a multipartite system.
        If a Python scalar or a `sympy.Integer` is input, 
        it is turned into the corresponding `sympy.Symbol`
        which is the only one to behave nicely as k in
        the current version of SymPy (1.1.13).
    
    Returns
    -------
    
    b : AnnihilateBoson
        Boson annihilation object.
        
    bd : CreateBoson
        Boson creation object.
    
    """
    
    if not(isinstance(k, Symbol)):
        k = Symbol(r"%s" % (k))
    
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
    return ("CreateBoson" in srepr(q)) \
            or ("AnnihilateBoson" in srepr(q))