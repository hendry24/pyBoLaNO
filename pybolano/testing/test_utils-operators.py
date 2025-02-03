import pytest
import random

from sympy import Symbol, I, Mul, Pow
from sympy.physics.secondquant import AnnihilateBoson, CreateBoson

from pybolano.utils.operators import (ops, 
                                      is_ladder, 
                                      is_ladder_contained, 
                                      random_ladder, 
                                      get_ladder_attr, 
                                      separate_mul_by_sub)

###

@pytest.mark.order(1)
def test_ops():
    b, bd = ops()
    assert isinstance(b, AnnihilateBoson) and isinstance(bd, CreateBoson)
    
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
def test_get_ladder_attr():
    for _ in range(10):
        exp = random.randrange(1, 10)
        sub = Symbol(f"{range(10)[random.randrange(10)]}")
        b, bd = ops(sub)
        assert (sub, exp) == get_ladder_attr(b**exp)

@pytest.mark.order(5)
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