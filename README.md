## **`boson_ladder`: Python Symbolic Package for Bosonic Ladder Operators**

<div align="center">
<img src="https://i.ibb.co.com/JKSzGdM/boson-ladder-logo.png" alt="boson-ladder-logo" style = "width:60%">
</div>

### **Hello, there!**

The `boson_ladder` package is _fully_ based on the `SymPy` symbolic package and serves as a complement for its bosonic ladder operator features that have yet to be implemented.

---

### **Who needs this?**

Are you tired of doing a **cumbersome bosonic ladder operator algebra with complicated commutation relations** and keep making mistakes in the page-long derivation but do not realize until you are at the end? We feel the same.

---

### **We present you...**

At `ver. 1.0.0`, this package offers you useful functions to do your algebra, including:

-   **Normal-ordering** _any_ polynomial of bosonic ladder operators.
-   Evaluating _any_ **commutation** relations of two polynomials of bosonic ladder operators.
-   And this is what motivates us to do this: evaluating _any_ **evolution equation for the
    expectation value** for _any_ system described in the **Linbdlad Master Equation** framework.

What's more, it works for **multipartite systems**!

---

### **Let's be transparent**

The core working principle of `boson_ladder` is simple&mdash;the package is based on the commutation relations $\left[\hat{b}_j , \hat{b}_k^\dagger\right]= 1 \mathrm{if} j=k,\ 0 \mathrm{otherwise}$ and $\left[\hat{b}_j,\hat{b}_k\right]=\left[\hat{b}_j^\dagger,\hat{b}_k^\dagger\right]=0$ of the bosonic creation $\hat{b}_j^\dagger$ and annihilation $\hat{b}_j^\dagger$ operators, where the subscript ($j$ here) indexes the bosonic mode. More precisely, we make use of the closed form of the commutator $\left[\hat{b}_j^{p},\hat{b}_k^{\dagger q}\right]$ where $p,q$ are integers to quickly evaluate a given commutator.

#### > [`normal_ordering`](https://github.com/hendry24/boson_ladder/blob/main/boson_ladder/core/normal_order.py#L80)

allows the user to normal-order any polynomial of bosonic ladder operators. It works like how we humans would normal-order a given expression. Assuming a sum (the most general form), the algorithm separates each addend by its subscripts. For each factor of one subscript, it looks for $\left(\hat{b}_j,\hat{b}_j^\dagger\right)$ sequences in the expression and applies the commutation relations to swap their places. This is done until all the creation operators are to the left of all the annihilation operators. The factors are then multiplied, and the products for different terms are summed to give an almost-normal-ordered expression. Lastly, the algorithm moves the operators with different indices (which commute) around to give a nice-looking output.

#### > [`do_commutator`](https://github.com/hendry24/boson_ladder/blob/main/boson_ladder/core/do_commutator.py#L162)

allows the user to evaluate any commutation relation of two polynomials of bosonic ladder operators. It uses commutator properties to express complex commutators in terms of commutators of the form $\left[\hat{b}_j^{p},\hat{b}_k^{\dagger q}\right]$, which can then be straightforwardly evaluated.

#### > [`LME_expval_evo`](https://github.com/hendry24/boson_ladder/blob/main/boson_ladder/core/Lindblad_ME.py#L123) 

allows the user to compute the expectation value evolution of a quantity represented by the operator $\hat{A}$ for a system described in the Lindblad master equation framework. The user simply needs to input: (1) the Hamiltonian $\hat{H}$; (2) the Lindblad dissipator operators $\hat{O}_j$ as well as their nonnegative multiplier $\gamma_j$; and (3) the operator $\hat{A}$ to calculate the expectation value evolution of.

Inside `LME_expval_evo`, the function [`Hamiltonian_trace`](https://github.com/hendry24/boson_ladder/blob/main/boson_ladder/core/Lindblad_ME.py#L21) is called to evaluate the contribution from the Hamiltonian, while [`dissipator_trace`](https://github.com/hendry24/boson_ladder/blob/main/boson_ladder/core/Lindblad_ME.py#L64) is called to evaluate the contribution from each dissipator term indexed $j$ above. These functions are available for the user to call, as well.

---

### **A quick guide**

We provide a quick tutorial of this package, in the file `boson_ladder_tutorial.ipynb` in the repository tree. Here is a quick [link](https://github.com/hendry24/boson_ladder/blob/main/boson_ladder_tutorial.ipynb) that will take you there. The notebook includes examples of use alongside a more detailed explanation of the way the package works.

---

### **Cite us, please!**

Pitiful as it may be, researchers nowadays are valued based on their citation counts. If you find our package helpful for your work, feel free to acknowledge the use of `boson_ladder` in your publications. Here is a `bibtex` entry you can copy:

```
@article{hendlim2024,
    title = "boson_ladder: Python Symbolic Package for Bosonic Ladder Operators",
    year = 2024,
    author = "Lim, Hendry M. and Dwiputra, Donny and Ukhtary, M. Shoufie and Nugraha, A. R. T"
}
```

---

### **Parting words**

This program is far from perfect and we would appreciate any critics and suggestions that you may have. In particular, we would appreciate it if you could inform us of any bugs you encounter while using this package. Feel free to reach out to us via email to [hendry01@ui.ac.id](mailto:hendry01@ui.ac.id).

Enjoy the package. \\( ﾟヮﾟ)/ 

\- The authors.

---

---

### **CHANGELOG**

`ver 1.0.0`
-   First release.