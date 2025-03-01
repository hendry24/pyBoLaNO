import pytest

from sympy import Add, Mul, Pow

from pybolano.utils.operators import ops
from pybolano.utils.expval import _expval

@pytest.mark.order(8)
def test_expval():
    b, bd = ops()
    
    x = _expval(b)
    print(x)
    
    try:
        x + x
        2*x
        x**2
        
        s = "_expval is for display only. "
        s += "Algebraic operations should not be allowed."
        raise AttributeError(s)
    except:
        pass
    
    #
    
    assert isinstance(_expval(b+bd), Add)
    assert isinstance(_expval(2*b), Mul)
    assert not(isinstance(_expval(b*bd), Mul))
    assert not(isinstance(_expval(b**2), Pow))