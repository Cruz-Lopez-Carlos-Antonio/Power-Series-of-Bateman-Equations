---
layout: default
title: Numerical Validation
math: true
---

## Purpose

The goal of this section is to document how the analytical solutions for \(n(t)\) and \(C(t)\) are validated against a high-precision numerical reference obtained with a fourth–order Runge–Kutta (RK4) method.

---

## RK4 Reference Solver

The script `RK4_reference_mpmath.py` integrates the NPKE system

$$
\frac{dn}{dt}=\frac{\rho(t)-\beta}{\Lambda}\,n(t)+\lambda\,C(t)+q,
\qquad
\frac{dC}{dt}=\frac{\beta}{\Lambda}\,n(t)-\lambda\,C(t),
$$

with \(\rho(t)=at+b\), using a sufficiently small time step and a working precision of 32 decimal digits.

The initial conditions are

$$
n(0) = \frac{q\,\Lambda}{|\rho(0)|}, \qquad
C(0) = \frac{\beta}{\lambda\Lambda}\,n(0),
$$

corresponding to a stationary state at \(t=0\) for the initial reactivity value.

---

## Error Measures

To compare the analytical and numerical solutions, absolute percentage errors (APE) are computed, for example:

$$
\operatorname{APE}_n(t_k)
= 100\,\frac{\bigl|n_{\text{RK4}}(t_k)-n_{\text{analytic}}(t_k)\bigr|}
               {\bigl|n_{\text{RK4}}(t_k)\bigr|},
$$

and analogously for the precursor concentration \(C(t)\).

Tables similar to those reported in the manuscript can be reproduced by:

1. Running `RK4_reference_mpmath.py` to generate the reference solution.
2. Running `Neutron_density_SciPyNumPy.py` and `C_precursor_SciPyNumPy.py` for the same time grid.
3. Computing the APE values and exporting them to a CSV or LaTeX table.

Further details on the numerical tests (ramp parameters, time grids and tolerance values) can be found in the article.
