## **`pyBoLaNO`: A `Python` Symbolic Package for Normal Ordering Involving Bosonic Ladder Operators**

<div align="center">
<img src="logo.svg" alt="boson-ladder-logo" style = "width:40%">
</div>

### **Hello, there!**

The `pyBoLaNO` package is _fully_ based on the `SymPy` symbolic package and serves as a complement for its bosonic ladder operator features that have yet to be implemented.

---

### **Who needs this?**

Are you tired of doing a **cumbersome bosonic ladder operator algebra with complicated commutation relations** and keep making mistakes in the page-long derivation but do not realize until you are at the end? We are.

---

### **We present you...**

At `ver. 1.0.2`, this package offers you useful functions to do your algebra, including:

-   **Normal-ordering** _any_ polynomial of bosonic ladder operators.
-   Normal-ordering _any_ **commutator** of two polynomials of bosonic ladder operators.
-   And this is what motivates us to do this: evaluating the normal ordering of _any_ 
    **expectation value evolution equation** for systems described in the 
    **Linbdlad Master Equation** framework using bosonic ladder operators.

What's more, it works for **multipartite systems**!

---

### **Installation**

You can install the package via `pip`:
```
pip install pybolano
```
or
```
pip install git+https://github.com/hendry24/pyBoLaNO
```

---

### **Let's be transparent**

The core working principle of `pyBoLaNO` is simple; the package is based on the commutation relations, i.e.,
$[\hat{b}_j,{\hat{b}_k}^\dagger]= 1~\textrm{if}j = k \quad~(\textrm{otherwise:}~[\hat{b}_j,{\hat{b}_k}^\dagger]= 0)$ and $\left[\hat{b}_j,\hat{b}_k\right]=\left[\hat{b}_j^\dagger,\hat{b}_k^\dagger\right]=0$ for the bosonic creation $\hat{b}_j^\dagger$ and annihilation $\hat{b}_j^\dagger$ operators, where the subscript ($j$ here) indexes the bosonic mode. More precisely, we make use of the explicit formula for the normal ordering of any monomial in bosonic operators presented by Blasiak (see arXiv links for [his PhD thesis](https://arxiv.org/abs/quant-ph/0507206) and a relevant [journal article](https://doi.org/10.1063/1.1990120)).

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

### **Cite us, please!**

Pitiful as it may be, researchers nowadays are valued based on their citation counts. If you find our package helpful for your work, feel free to acknowledge the use of `pyBoLaNO` in your publications. Here is a `bibtex` entry you can copy:

```
@misc{lim2025pybolano,
      title={pyBoLaNO: A Python symbolic package for normal ordering involving bosonic ladder operators}, 
      author={Hendry M. Lim and Donny Dwiputra and M. Shoufie Ukhtary and Ahmad R. T. Nugraha},
      year={2025},
      eprint={2501.01603},
      archivePrefix={arXiv},
      primaryClass={quant-ph},
      url={https://arxiv.org/abs/2501.01603},
}
```

If you are interested in how we handle the peer review from the original version to the latest version of the preprint manuscript, we also provide a [preprint folder](preprint) in this repository in addition to the [arXiv page](https://arxiv.org/abs/2501.01603).

---

### **Parting words**

This program is far from perfect and we would appreciate any critics and suggestions that you may have. In particular, we would appreciate it if you could inform us of any bugs you encounter while using this package. Feel free to reach out to us via email to [hendry01@ui.ac.id](mailto:hendry01@ui.ac.id).

Enjoy the package. \\( ﾟヮﾟ)/ 

\- The authors.