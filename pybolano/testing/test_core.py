import pytest
import random

from sympy import Symbol, I, symbols, E, Derivative

from pybolano.core.commutators import _break_comm
from pybolano.core import normal_ordering, NO, LME_expval_evo
from pybolano.utils.operators import ops, random_ladder
from pybolano.utils.expval import _expval

###

@pytest.mark.order(7)
def test_break_comm():
    b, bd = ops()
    for _ in range(100):
        qs = [1, 1]
        for i in range(2):
            for _ in range(random.randrange(1, 11)):
                if random.randrange(2):
                    qs[i] *= b
                if random.randrange(2):
                    qs[i] *= bd
        
        A, B = qs
        
        expect = (A*B-B*A).expand()
        
        assert _break_comm(A, B) == expect

@pytest.mark.order(8)
def test_NO():
    b_1, bd_1 = ops(1) 
    b_2, bd_2 = ops(2)
    x = Symbol("x")
    
    for _ in range(100):
        q = [b_1, bd_1][random.randrange(2)]**random.randrange(10) * x**random.randrange(5)
        assert q == normal_ordering(q)
        
        q = bd_1**random.randrange(10) * b_1**random.randrange(10) * x**random.randrange(5)
        assert q == normal_ordering(q)
    
    assert (bd_1*b_2) == normal_ordering(bd_1*b_2)
    assert (bd_2*b_1) == normal_ordering(b_1*bd_2)
    
    assert (1+bd_1*b_1) == normal_ordering(b_1*bd_1)
    
    assert (2+bd_1*b_1+bd_2*b_2) == normal_ordering(b_1*bd_1 + b_2*bd_2)
    
    # nonnegative excess
    assert (bd_1**2*b_1 + 2*bd_1) == normal_ordering(b_1*bd_1**2)
    
    # negative excess
    assert (2*b_1+bd_1*b_1**2) == normal_ordering(b_1**2 * bd_1)
    
    # alias
    for _ in range(100):
        q = random_ladder(10, ["a", "b", "c"])
        assert NO(q) == normal_ordering(q)
    
@pytest.mark.order(9)
def test_LME_expval_evo():
    b_1, bd_1 = ops("1")
    b_2, bd_2 = ops("2")
    
    #
    
    res = LME_expval_evo(H = bd_1*b_1,
                         D = [],
                         A = b_1)
    assert res.lhs == Derivative(_expval(b_1), Symbol("t"))
    
    got = res.rhs
    expect = _expval(-I * b_1)
    assert got == expect
    
    #
    
    got = LME_expval_evo(H = bd_1*b_1,
                         D = [],
                         A = b_2).rhs
    expect = _expval(0)
    assert got == expect
    
    # Europhys. Lett. 140 (3) (2022) 35001, Eq. [1--4].
    # Also available in the tutorial notebook as the last example.
    # Here we only check for the evolution of b_1.
        
    Delta, Omega, g, theta = symbols("Delta Omega g theta")
    gamma, Gamma, phi = symbols("gamma Gamma phi")

    H = Delta*(bd_1*b_1 + bd_2*b_2) \
        + Omega*(b_1+bd_1) \
        + g*(E**(I*theta)*bd_1*b_2 + E**(-I*theta)*bd_2*b_1)

    D = [[gamma, b_1],
        [gamma, b_2],
        [Gamma*E**(I*phi), b_2, b_1],
        [Gamma*E**(-I*phi), b_1, b_2]]

    A = b_1

    got = LME_expval_evo(H=H, D=D, A=A).rhs
    
    expect = _expval(
        - I * (Delta - I/2*gamma) * b_1
        - (I*g*E**(I*theta) + Gamma*E**(I*phi)/2) * b_2
        - I * Omega
    )
    
    assert got == expect