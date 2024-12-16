from sympy import \
    Expr
    
from sympy import *
from sympy.physics.secondquant import *
from boson_ladder import *

__all__ = ["ladder",
           "ladders"]
    
class Ladder(Expr):
    def __new__(cls, *args, **kwargs):
        return super(Ladder, cls).__new__(cls)
    
    def __init__(self, sub = None, dagger = False, symbol = "a"):  
        super().__init__()
        
        s = symbol
        if sub is not None:
            s += r"_{%s}"%(sub)
        if dagger:
            s += r"^{\dagger}"
        self.s = s
        self.sympy = Symbol(s)
        
        self.symbol = symbol
        self.sub = sub
        self.dagger = dagger
        
    def dag(self):
        return Ladder(sub=self.sub,
                      dag = not(self.dagger), 
                      symbol=self.symbol)
        
    def _latex(self, printable):
        return r"{%s}"%(printable.doprint(self.sympy))
    
    def __str__(self):
        return self.sympy.__str__()
    
    def __repr__(self):
        return self.sympy.__repr__()
        
    def __add__(self, x):
        if is_ladder(x):
            x = x.sympy
        return Add(self.sympy, x)
    
    def __sub__(self, x):
        return self.__add__(-x)
    
    def __mul__(self, x):
        if is_ladder(x):
            x = x.sympy
        return Mul(self.sympy, x)
    
    def __pow__(self, x):
        if is_ladder(x):
            x = x.sympy
        return Pow(self.sympy, x)
    
def ladders(sub = None, symbol = "a"):
    a = Ladder(sub=sub, symbol=symbol, dagger=False)
    ad = Ladder(sub=sub, symbol=symbol, dagger=True)
    return a, ad