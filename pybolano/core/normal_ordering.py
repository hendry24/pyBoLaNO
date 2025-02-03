from ast import Expr
from multiprocessing import Pool

from sympy import Add, FallingFactorial, Mul, Number, Pow, binomial, factorial, sympify
from sympy.physics.secondquant import (
    AnnihilateBoson,
    CreateBoson,
)

from pybolano.utils.error_handling import InvalidTypeError
from pybolano.utils.multiprocessing import mp_config
from pybolano.utils.operators import (
    get_ladder_attr,
    is_ladder,
    is_ladder_contained,
    separate_mul_by_sub,
)

############################################################

__all__ = ["normal_ordering", "NO"]

############################################################


def _NO_Blasiak(q):
    """
    Normal ordering with the explicit formula derived in
    Chapter 4 of Blasiak's PhD Dissertation, available at
    https://arxiv.org/abs/quant-ph/0507206.

    q is assumed to be a ``boson string'' per Blasiak, i.e.
    a monomial in the bosonic ladder operators. Refer to
    Eqs. (4.2), (4.10), (4.34).

    Input is assumed to contain single subscript.
    """

    if not (is_ladder_contained(q)) or is_ladder(q) or isinstance(q, Pow):
        return q

    if isinstance(q, Mul):
        q_args = q.args
    else:
        raise InvalidTypeError([AnnihilateBoson, CreateBoson, Pow, Mul], type(q))

    ###

    r = (
        [] if q_args[0].has(CreateBoson) else [Number(0)]
    )  # in case monomial starts with b
    s = []
    for qq in q_args:
        sub, exp = get_ladder_attr(qq)
        if qq.has(CreateBoson):
            r.insert(0, exp)
        else:
            s.insert(0, exp)
    if len(r) != len(s):  # monomial ends with bd
        s.insert(0, Number(0))

    # To make indexing easier, we pad r and s
    # with r_0 and s_0,
    # which do not exist in Blasiak's formulation.
    r.insert(0, Number(0))
    s.insert(0, Number(0))

    # Excess
    d = []  # d_0 is, however, used in Eq. (4.10).
    sum_val = Number(0)
    for r_m, s_m in zip(r, s):
        sum_val += r_m - s_m
        d.append(sum_val)

    ###

    def _S_rs(s, d, k):
        """
        Generalized Stirling number, Eq. (4.10). We use
        d instead of r since d is already calculated
        before this function is called.
        """
        sum_val = Number(0)
        for j in range(k + 1):
            prod_val = Number(1)
            for m in range(1, len(s)):
                prod_val *= FallingFactorial(d[m - 1] + j, s[m])
            sum_val += binomial(k, j) * Number(-1) ** (k - j) * prod_val
        return 1 / factorial(k) * sum_val

    ###

    b = list(q.find(AnnihilateBoson))[0]
    bd = list(q.find(CreateBoson))[0]

    if d[-1] >= 0:
        R, S, D = r, s, d
        k_lst = range(s[1], sum(s) + 1)
    else:
        k_lst = range(r[-1], sum(r) + 1)

        """
        Somehow using the original expression
        in Eq. (4.34) does not work. However,
        it does work when we utilize the symmetry
        property stated in Eq. (4.37).
        """
        R = [Number(0)] + list(reversed(s[1:]))
        S = [Number(0)] + list(reversed(r[1:]))
        D = []
        sum_val = Number(0)
        for r_m, s_m in zip(R, S):
            sum_val += r_m - s_m
            D.append(sum_val)

    out = Number(0)
    for k in k_lst:
        out += _S_rs(S, D, k) * bd**k * b**k

    if d[-1] >= 0:
        out = (bd ** d[-1] * out).expand()
    else:
        out = (out * b ** (-d[-1])).expand()

    return out


def _NO_preprocess(q):
    q = q.expand()

    if (not (q.has(CreateBoson)) and not (q.has(AnnihilateBoson))) or isinstance(
        q, (Pow, CreateBoson, AnnihilateBoson)
    ):
        return q

    elif isinstance(q, Add):
        q_args = [qq for qq in q.args]

    elif isinstance(q, Mul):
        if not (is_ladder_contained(q)):
            return q
        q_args = [q]

    else:
        raise InvalidTypeError([Pow, CreateBoson, AnnihilateBoson, Add, Mul], type(q))
    return q_args


def _final_swap(q):
    """
    Each input must be a summand. One subscript
    must have exactly one Pow of b and one Pow
    of bd.
    """
    if not (isinstance(q, Mul)):
        return q
    else:
        collect_scalar = []
        collect_b = {}
        collect_bd = {}
        for qq in q.args:
            if not (is_ladder_contained(qq)):
                collect_scalar.append(qq)
                continue

            sub = get_ladder_attr(qq)[0].name
            if (sub not in collect_b) or (sub not in collect_bd):
                collect_b[sub] = Number(1)
                collect_bd[sub] = Number(1)
                # to force the subscript order for
                # b and bd.

            if qq.has(AnnihilateBoson):
                collect_b[sub] = qq
            else:
                collect_bd[sub] = qq

        collect_b = {k: collect_b[k] for k in sorted(collect_b)}
        collect_bd = {k: collect_bd[k] for k in sorted(collect_bd)}

        return Mul(
            *(collect_scalar + list(collect_bd.values()) + list(collect_b.values()))
        )


def _NO_input_addend(qq):
    qq_Mul_args = separate_mul_by_sub(qq)

    out = Number(1) if is_ladder_contained(qq_Mul_args[0]) else qq_Mul_args.pop(0)

    for qq_single_sub in qq_Mul_args:
        out *= _NO_Blasiak(qq_single_sub)

    return out.expand()


############################################################


def normal_ordering(q: Expr) -> Expr:
    """
    Normal order the operator q: all creation operators are
    written to the left of all annihilation operators within
    a single term. This function uses Blasiak's formulae
    [see Eqs. (4.2), (4.10), (4.34), (4.37) in https://arxiv.org/abs/quant-ph/0507206 ].

    Parameters
    ----------

    q : sympy.Expr
        An expression containing SymPy's bosonic ladder operator objects.

    Returns
    -------

    q_NO : sympy.Expr
        Normal ordering of q.

    See Also
    --------

    ops : Get the bosonic ladder operator objects.
    """

    q = sympify(q)

    # Shortcuts

    if not (is_ladder_contained(q)) or isinstance(
        q, (Pow, CreateBoson, AnnihilateBoson)
    ):
        return q

    # NOTE: checking if all bd are already to the
    # left of b may be too computationally expensive
    # to implement here.

    ###

    q_args = _NO_preprocess(q)
    use_mp = mp_config["enable"] and (len(q_args) >= mp_config["min_num_args"])

    ###

    if use_mp:
        with Pool(mp_config["num_cpus"]) as pool:
            _out = Add(*pool.map(_NO_input_addend, q_args))
    else:
        _out = Add(*[_NO_input_addend(qq) for qq in q_args])

    """
    At this point, the normal ordering is not done since there are
    probably the terms are written something like 
        bd_0**p*b_0**q * bd_1**r*b_1**s
    The last step is simply to swap the argument order to get the
    nice-looking output with the same subscript order between
    the creation and the annihilation operators.
    """

    if not(isinstance(_out, Add)):
        return _final_swap(_out)

    if use_mp:
        with Pool(mp_config["num_cpus"]) as pool:
            return Add(*pool.map(_final_swap, _out.args))
    else:
        return Add(*[_final_swap(q) for q in _out.args])


def NO(q: Expr) -> Expr:
    """
    Alias for `normal_ordering`.
    """
    return normal_ordering(q)
