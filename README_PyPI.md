The core working principle of `pyBoLaNO` is simple&mdash;the package is based on the commutation relations $\left[\hat{b}_j , \hat{b}_k^\dagger\right]= 1 \mathrm{if} j=k,\ 0 \mathrm{otherwise}$ and $\left[\hat{b}_j,\hat{b}_k\right]=\left[\hat{b}_j^\dagger,\hat{b}_k^\dagger\right]=0$ of the bosonic creation $\hat{b}_j^\dagger$ and annihilation $\hat{b}_j^\dagger$ operators, where the subscript ($j$ here) indexes the bosonic mode. More precisely, we make use of the explicit formula for the normal ordering of any monomial in bosonic operators presented by Blasiak ([arXiv link for his PhD thesis](https://arxiv.org/abs/quant-ph/0507206) and the [journal article](https://doi.org/10.1063/1.1990120)).

#### > [`normal_ordering`](https://github.com/hendry24/pyBoLaNO/blob/main/pybolano/core/normal_ordering.py#L221)

allows the user to normal-order any polynomial of bosonic ladder operators. It separates each monomial in the input (most generally a polynomial) by the subscripts of the ladder operators. For each subscript, normal ordering is performed using Blasiak's formulae (see Eqs. (4.2), (4.10), (4.34), (4.37) of his thesis linked above). Lastly, the algorithm moves the operators with different indices (which commute) around to give a nice-looking output.

#### > [`NO_commutator`](https://github.com/hendry24/pyBoLaNO/blob/main/pybolano/core/commutators.py#L35)

allows the user to evaluate the any commutation relation of two polynomials of bosonic ladder operators.
It is just a shorthand to save you the time of typing `normal_ordering(A*B-B*A)`.

#### > [`LME_expval_evo`](https://github.com/hendry24/pyBoLaNO/blob/main/pybolano/core/Lindblad_ME.py#L115) 

allows the user to compute the normal-ordered expression for the expectation value evolution of a quantity represented by the operator $\hat{A}$ for a system described in the Lindblad master equation framework. The user simply needs to input: (1) the Hamiltonian $\hat{H}$; (2) the Lindblad dissipator operators $\hat{O}_j,\hat{P}_j$ as well as their nonnegative multiplier $\gamma_j$; and (3) the operator $\hat{A}$ to calculate the expectation value evolution of.

Inside `LME_expval_evo`, the function [`Hamiltonian_trace`](https://github.com/hendry24/boson_ladder/blob/main/boson_ladder/core/Lindblad_ME.py#L21) is called to evaluate the contribution from the Hamiltonian, while [`dissipator_trace`](https://github.com/hendry24/boson_ladder/blob/main/boson_ladder/core/Lindblad_ME.py#L64) is called to evaluate the contribution from each dissipator term indexed $j$ above. These functions are available for the user to call, as well.

---

### **A quick guide**

We provide a quick tutorial of this package, in the file `tutorial.ipynb` in the repository tree. Here is a quick [link](https://github.com/hendry24/pyBoLaNO/blob/main/tutorial.ipynb) that will take you there. The notebook includes examples of use alongside a more detailed explanation of the way the package works.

---
