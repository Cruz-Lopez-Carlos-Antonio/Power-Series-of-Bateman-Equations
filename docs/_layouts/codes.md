---
layout: default
title: Python Codes
---

## Overview of the Python Scripts

The repository contains seven main Python 3 scripts.  
They implement the analytical solution for $$n(t)$$ and the associated delayed neutron precursor concentration $$C(t)$$, both developed with the Modified Integration Method proposed in the submited paper. A RK4 reference solution, as well as the computational implementation of the Zhang et al. (2008) and the Palma et al. (2010) solutions. 

---

### 1. `Neutron_density_SciPyNumPy.py`
<div style="padding:8px; border-left:4px solid #3c6e71; margin-bottom:10px;">
  <a href="https://github.com/Cruz-Lopez-Carlos-Antonio/Ramp_analytical_solution/blob/main/Neutron_density_SciPyNumPy.py" 
     target="_blank" style="font-size:16px; color:#22577a; font-weight:bold;">
     👉 Click here to view the code in a new tab
  </a>
</div>

This script implements the analytical solution of the neutron density $$n(t)$$ using **SciPy** and **NumPy**.  
The core of the implementation is based on the integral representation shown in the *Equations* page, where the involved integrals are evaluated using `scipy.integrate.quad`.

The script also contains a linear system, derived from the initial conditions $$n(0)$$ and $$\dot n(0)$$, to determine the constants $$K_1$$ and $$K_2$$ via a least-squares procedure (`numpy.linalg.lstsq`), including column-wise normalization for numerical stability.

---

### 2. `Neutron_density_mpmath.py`
[Link to access to the code](https://github.com/Cruz-Lopez-Carlos-Antonio/Ramp_analytical_solution/blob/main/Neutron_density_mpmath.py)
This script provides a **high-precision** version of the analytical solution for $$n(t)$$, implemented with **mpmath**.  
It includes:

- Multiprecision evaluation of the integrals $$I_1,\dots,I_6$$ described in the manuscript,
- A robust $$2\times 2$$ linear solver with row/column scaling and Tikhonov regularization,
- Control of the working precision through `mp.mp.dps`.

This implementation is used as a benchmark to assess conditioning effects and to validate the double-precision results obtained with SciPy/NumPy.

---

### 3. `C_precursor_SciPyNumPy.py`

This script computes the delayed neutron precursor concentration $$C(t)$$ using the convolution formula

$$
C(t)
= C(0)\,e^{-\lambda t}
+ \frac{\beta}{\Lambda}\,e^{-\lambda t}
\int_0^t e^{\lambda\tau}\,n(\tau)\,d\tau.
$$

The integral is evaluated numerically using `scipy.integrate.cumulative_trapezoid`.  
The function $$n(t)$$ is imported from `Neutron_density_SciPyNumPy.py`, and the script returns a vectorized approximation of $$C(t)$$ over a prescribed time grid.

---

### 4. `RK4_reference_mpmath.py`

This script implements a **fourth–order Runge–Kutta (RK4)** solver in 32-digit precision (via mpmath) for the NPKE system.  

It solves simultaneously for \(n(t)\) and \(C(t)\) using a fine time step, and the resulting numerical solution is used as a high-accuracy reference to validate the analytical formulations and their numerical implementation.

Parameters, time step, and integration interval can be adjusted to reproduce the tables and figures reported in the manuscript.
