---
layout: default
title: Overview
math: true
---
# Overview

The present repository documents the computational implementation of a new analytical solution of the following system of differential equations: 

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\begin{matrix}\dfrac{dn\left(t\right)}{dt}&=&\dfrac{\rho\left(t\right)-\beta}{\Lambda}n\left(t\right)+\lambda C\left(t\right)+q,\\\dfrac{dC\left(t\right)}{dt}&=&\dfrac{\beta}{\Lambda}n\left(t\right)-\lambda C\left(t\right).\\\end{matrix}\, 
$$
</div>

known as the Neutron Point Kinetic Equations (NPKE). Such analytical solution corresponds to the particular case of a linear ramp reactivity, $$\rho(t)=at+b$$ and was obtained using the Modified Integral Method (MIM). 

## Python Codes and Implementations

<div style="margin: 1.5rem 0; padding: 1.25rem 1.5rem; background:#fafafa; border-radius:10px; border:1px solid #e0e0e0;">
  <h3 style="margin-top:0; margin-bottom:0.75rem; font-size:1.1rem;">
    Python 3 implementations included in this repository:
  </h3>

  <ul style="margin:0; padding-left:1.2rem; list-style-type:disc;">
    <li style="margin:0.4rem 0;">
      <strong>Analytical neutron density \( n(t) \)</strong><br/>
      <span>Implementation using <code>SciPy</code>/<code>NumPy</code> libraries.</span>
    </li>
    <li style="margin:0.4rem 0;">
      <strong>High-precision analytical neutron density \( n(t) \)</strong><br/>
      <span>Implementation using <code>mpmath</code> with <em>32-digit precision</em> by default.</span>
    </li>
    <li style="margin:0.4rem 0;">
      <strong>Delayed neutron precursor concentration \( C(t) \)</strong><br/>
      <span>Analytical evaluation consistent with the neutron density solvers.</span>
    </li>
    <li style="margin:0.4rem 0;">
      <strong>RK4 reference solver</strong><br/>
      <span>A high-precision classical Runge–Kutta 4 (RK4) scheme used as benchmark.</span>
    </li>
    <li style="margin:0.4rem 0;">
      <strong>Analytical benchmark solution by Zhang et al. (2008)</strong><br/>
      <span>Closed-form formulation for comparison and validation.</span>
    </li>
    <li style="margin:0.4rem 0;">
      <strong>Analytical benchmark solution by Palma et al. (2010)</strong><br/>
      <span>Alternative analytical solution for cross-checking the proposed method.</span>
    </li>
    <li style="margin:0.4rem 0;">
  <strong>Comparative execution-time benchmark</strong><br/>
  <span>
    A dedicated performance routine that evaluates and compares the mean execution time of all methods.  
    The implementation makes use of <code>time</code>,to quantify the computational cost of each approach.
  </span>
</li>
  </ul>
</div>

## The related article

These results accompany the manuscript:  

<div style="display:inline-block; padding:0.6em 0.9em; margin:0.4em 0;
            background-color:#f5f5f5; border:1px solid #dddddd; border-radius:6px;"
     markdown="1">
***A New Analytical Solution to the Neutron Point Kinetics Equations for Linear Ramp Reactivity: $$\rho(t)=a t + b$$. A Theoretical and Computational Framework***
</div>

that was recently submitted to the journal *Computer Physics Communications*.

# Map of the site 
You can use the navigation bar above, or the following buttons, to explore:

<div style="display:grid; grid-template-columns:repeat(auto-fit,minmax(220px,1fr)); gap:0.9rem; margin:1rem 0;">

  <a href="{{ '/equations.html' | relative_url }}" style="text-decoration:none; color:inherit;">
    <div style="background:#fafafa; border:1px solid #e0e0e0; border-radius:10px; padding:0.75rem 0.9rem; cursor:pointer;">
      <strong>Equations</strong><br/>
      <span style="font-size:0.95rem; color:#555;">
        Main analytical expressions and key formulae used in the article.
      </span>
    </div>
  </a>

  <a href="{{ '/codes.html' | relative_url }}" style="text-decoration:none; color:inherit;">
    <div style="background:#fafafa; border:1px solid #e0e0e0; border-radius:10px; padding:0.75rem 0.9rem; cursor:pointer;">
      <strong>Codes</strong><br/>
      <span style="font-size:0.95rem; color:#555;">
        Summary of the Python scripts, interfaces, and numerical settings.
      </span>
    </div>
  </a>

  <a href="{{ '/validation.html' | relative_url }}" style="text-decoration:none; color:inherit;">
    <div style="background:#fafafa; border:1px solid #e0e0e0; border-radius:10px; padding:0.75rem 0.9rem; cursor:pointer;">
      <strong>Validation</strong><br/>
      <span style="font-size:0.95rem; color:#555;">
        Comparison against the high-precision RK4 benchmark and analytical references.
      </span>
    </div>
  </a>

  <a href="{{ '/about.html' | relative_url }}" style="text-decoration:none; color:inherit;">
    <div style="background:#fafafa; border:1px solid #e0e0e0; border-radius:10px; padding:0.75rem 0.9rem; cursor:pointer;">
      <strong>About</strong><br/>
      <span style="font-size:0.95rem; color:#555;">
        Authorship, affiliations, and financial support acknowledgements.
      </span>
    </div>
  </a>

</div>


<div style="border-left: 4px solid #f39c12; padding: 0.7em 1em; background: #fff7e6;">
<b style="color:#c0392b;">⚠️  Important:</b><br>
For more details of the derivation of the analytical solution, its computational implementation as well as a numerical analysis, please see the manuscript mentioned above.
</div>
