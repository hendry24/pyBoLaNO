import pytest
import random

from sympy import Symbol, I, Mul, Pow, sympify
from sympy.physics.quantum.boson import BosonOp

from pybolano.utils.operators import (BosonicAnnihilationOp,
                                      BosonicCreationOp,
                                      ops, 
                                      is_ladder,
                                      is_ladder_contained, 
                                      dagger,
                                      random_ladder, 
                                      get_ladder_attr, 
                                      separate_mul_by_sub,
                                      to_sympy_physics)

###

@pytest.mark.order(1)
def test_ops():
    b, bd = ops()
    assert isinstance(b, BosonicAnnihilationOp) and isinstance(bd, BosonicCreationOp)
    
    ops("iashdak")
    ops(r"\\gamma")
    ops(2)
    ops(Symbol("x"))

@pytest.mark.order(2)
def test_ladder_check():
    x = Symbol("x")
    b, bd = ops()
    
    assert all([not(is_ladder(I)),
                not(is_ladder(x)),
                is_ladder(b),
                is_ladder(bd),
                not(is_ladder(b+bd)),
                not(is_ladder(2*b)),
                not(is_ladder(b*bd)),
                not(is_ladder(b**2))])
    
    assert all([not(is_ladder_contained(I)),
                not(is_ladder_contained(x)),
                is_ladder_contained(b),
                is_ladder_contained(bd),
                is_ladder_contained(b+bd),
                is_ladder_contained(2*b),
                is_ladder_contained(b*bd),
                is_ladder_contained(b**2)])

@pytest.mark.order(3)
def test_random_ladder():
    assert not(is_ladder_contained(random_ladder(0)))
    assert is_ladder_contained(random_ladder(5))
    assert is_ladder_contained(random_ladder(5, Symbol("0")))
    assert is_ladder_contained(random_ladder(5, [Symbol("x"), "yeah", 0]))

@pytest.mark.order(4)
def test_dagger():
    for _ in range(100):
        q = random_ladder(10)
        assert q == dagger(dagger(q))

@pytest.mark.order(5)
def test_get_ladder_attr():
    for _ in range(10):
        exp = random.randrange(1, 10)
        sub = Symbol(f"{range(10)[random.randrange(10)]}")
        b, bd = ops(sub)
        assert (sub, exp) == get_ladder_attr(b**exp)

@pytest.mark.order(6)
def test_separate_by_sub():
    q = random_ladder(20, ["a", "b", "c", "d"])
    q_sep = separate_mul_by_sub(q)
    sub_checked = []
    for qq in q_sep:
        sub_lst = []
        
        if not(is_ladder_contained(qq)):
            continue
        
        elif is_ladder(qq) or isinstance(qq, Pow):
            sub = get_ladder_attr(qq)[0]
            if sub in sub_checked:
                s = "Same sub in different entries."
                raise ValueError(s)
            sub_lst.append(sub)
            sub_checked.append(sub)
        
        elif isinstance(qq, Mul):
            for i, qqq in enumerate(qq.args):
                sub = get_ladder_attr(qqq)[0]
                if i == 0:
                    if sub in sub_checked:
                        s = "Same sub in different entries."
                        raise ValueError(s)
                    sub_checked.append(sub)
                sub_lst.append(sub)
            
        else:
            s = "Expected scalar, SqOperator, Pow, or Mul."
            raise ValueError(s)
        
        assert len(set(sub_lst)) == 1
        
@pytest.mark.order(7)
def test_to_sympy_physics():
    assert to_sympy_physics(5) == sympify(5)
    assert to_sympy_physics(I) == I
    assert to_sympy_physics(Symbol("x")) == Symbol("x")

    b, bd = ops(1)
    bb, bbd = ops(2)
    
    b_s = BosonOp(r"b_{1}", True)
    bd_s = BosonOp(r"b_{1}", False)
    bb_s = BosonOp(r"b_{2}", True)
    bbd_s = BosonOp(r"b_{2}", False)
    
    assert to_sympy_physics(b) == b_s
    assert to_sympy_physics(bbd) == bbd_s
    assert to_sympy_physics(b+1) == b_s+1
    assert to_sympy_physics(b+bbd) == b_s + bbd_s
    assert to_sympy_physics(50*Symbol(r"\gamma")*bb)  == 50*Symbol(r"\gamma")*bb_s
    assert to_sympy_physics(bd**2) == bd_s**2
    assert to_sympy_physics(5*bbd**2*b + b*bbd*b) == (5*bbd_s**2*b_s + b_s*bbd_s*b_s)