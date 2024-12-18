## **`boson_ladder`: Python Symbolic Package for Bosonic Ladder Operators**

---

### **Hello, there.**

The `boson_ladder` package is _fully_ based on the `SymPy` symbolic package, and serves as a complement for its bosonic ladder operator features that is not yet implemented.

### **Who needs this?**

Are you tired of doing a **cumbersome bosonic ladder operator algebra with complicated commutation relations** and keep on making mistakes in the page-long derivation but did not realize until you are at the end? We feel the same.

### **We present you...**

At `ver. 1.0.0`, this package offers you useful functions to do your algebra, including:

-   Evaluating _any_ **commutation** relations of two polynomials of bosonic ladder operators.
-   **Normal-ordering** _any_ polynomial of bosonic ladder operators.
-   And this is what motivates us to do this: evaluating _any_ **evolution equation for the
    expectation value** for _any_ system described in the **Linbdlad Master Equation** framework.

What's more, it works for **multipartite systems**!

### **Let's be transparent**

The core working principle of `boson_ladder` is simple. The package is based on the commutation relations $[b_j,b_k^\dagger]=\delta_{jk}$ and $[b_j,b_k]=[b_j^\dagger,b_k^\dagger]=0$ of the bosonic creation $b_j^\dagger$ and annihilation $b_j^\dagger$ operators, where the subscript ($j$ here) indexes the bosonic mode. 

#### > `do_commutator` 

allows the user to evaluate any commutation relation of two polynomials of bosonic ladder operators, based on the identities 
$$
[AB,C] = A[B,C] + [A,C]B
$$ 
$$
[A,BC] = [A,B]C + B[A,C]
$$
To evaluate a commutator, `boson_ladder` applies these identites recursively to expand the commutator into a sum of simpler commutators. This is done until the commutators can be automatically evaluated by SymPy.

#### > `normal_ordering` 

allows the user to normal order any polynomial of bosonic ladder operators. It works like how we humans would normal-order the operator: it looks for a creation operator to the right of an annihilation operator, then apply the commutation relation to swap their places. This is done recursively until all the creation operators are
to the left of all the annihilation operators. 

#### > `LME_expval_evo` 

allows the user to compute the equation for the evolution of the expectation value of a quantity represented by the operator $A$, i.e. 
$$
\frac{\mathrm{d}\left\langle A\right\rangle}{\mathrm{d}t} = \mathrm{tr}\left(-i\left[H,\rho\right]A+\sum_j\gamma_j \mathcal{D}\left(O_j\right)\left[\rho\right]A\right)
$$
where $H$ and $O_k$ can be any polynomials in the ladder operators. The following identities are used to compute the RHS (you can easily derive these):

\begin{align}
\mathrm{tr}\left(-i\left[H,\rho\right]A\right) &= -i\left[A,H\right]
\notag\nonumber
\\
\mathrm{tr}\left(\mathcal{D}\left(\sum_k O_k\right)\left[\rho\right]A\right) &= \sum_k \left\langle O_k^\dagger \left[A, O_k\right]\right\rangle + \left\langle \left[O_k^\dagger, A\right] O_k \right\rangle 
\notag\nonumber
\\ & \quad
+\sum_{l>k}\left\langle O_k^\dagger \left[A,O_l\right]\right\rangle + \left\langle \left[ O_k^\dagger,A\right]O_l\right\rangle + \left\langle O_l^\dagger \left[A,O_k\right]\right\rangle + \left\langle \left[ O_l^\dagger,A\right]O_k\right\rangle
\notag\nonumber
\end{align}

The user simply needs to input: (1) the Hamiltonian $H$; (2) the Lindblad dissipator (formally the Liouvillian superoperator in Lindblad form) operators $O_j$ and their multiplying nonnegative scalar $\gamma_j$; and (3) the operator $A$ to calculate the expectation value evolution of&mdash;the function will do the rest! 

Inside `LME_expval_evo`, the function `Hamiltonian_trace` is called to evaluate the contribution from the Hamiltonian, while `dissipator_trace` is called to evaluate the contribution from each dissipator term indexed $j$ above. These functions are available for the user to call, as well.

### **A quick guide**

We provide a quick tutorial of this package, in the file `boson_ladder_tutorial.ipynb` in the repository tree. Here is a quick ![link.](https://github.com/hendry24/boson_ladder/blob/main/boson_ladder_tutorial.ipynb)

### **Cite us, please!**

Pitiful as it may be, researchers nowadays are valued based on their citation counts. If you find our package helpful, feel free to acknowledge the use of `boson_ladder` in your publications. Here is a `bibtex` entry you can copy:

```
@article{hendlim2024,
    title = "boson_ladder: Python Symbolic Package for Bosonic Ladder Operators",
    year = 2024,
    author = "Lim, Hendry M. and Dwiputra, Donny and Ukhtary, M. Shoufie and Nugraha, A. R. T"
}
```

Enjoy the package. \\( ﾟヮﾟ)/

### CHANGELOG

`ver 1.0.0`
    - First release.