`v1.1.0`
-   Add: `pybolanoOp`, `BosonicCreationOp`, and `BosonicAnnihilationOp` to replace
SymPy's implementations (avoids future hassle with the possible merge conflict between `sympy.physics.secondquant` and `sympy.physics.quantum`).
-   Add: `dagger` returns the complex conjugate of the input, replacing Dagger in
`sympy.physics`.
-   Add: unit test for `dagger`.

`v1.0.2`
-   Add: unit tests.
-   Add: type annotations; change import to absolute (PR [#1](https://github.com/hendry24/pyBoLaNO/pull/1) by [@](https://github.com/Matt-Ord)).
-   Add: sympify in `_break_comm` and `_expval` to handle input better.
-   Update: Algebraic operations are prohibited for `_expval`.
-   Update: `random_ladder` now can also shuffle subscripts. 
-   Bugfix: no `_final_swap` call for non-`Add` result in `normal_ordering`.

`v1.0.1`
-   Fix PyPI deployment.

`v1.0.0`
-   First release.