---
layout: default
title: Main Equations
math: true
---
# Overview of the methodology
The present section contains the main equations used in the present work. Please, see the full manuscript for a more detailed description of them.

## Neutron Point Kinetics Equations
The Neutron Point Kinetic Equations with a single group of delayed neutron precursors are given as:

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\begin{matrix}\dfrac{dn\left(t\right)}{dt}&=&\frac{\rho\left(t\right)-\beta}{\Lambda}n\left(t\right)+\lambda C\left(t\right)+q,\\\dfrac{dC\left(t\right)}{dt}&=&\frac{\beta}{\Lambda}n\left(t\right)-\lambda C\left(t\right).\\\end{matrix}\ .  
$$
</div>

Here $$\Lambda$$ is the prompt generation time, $$\beta$$ is the total delayed neutron fraction and $$\lambda$$. In the present case, the ramp reactivity, $$\rho(t)$$ has the following linear form:

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\rho(t) = a t + b, \quad a>0.
$$
</div>

---

## Modified Integration Method (MIM) 

The proposed solutions were developed using a more efficient, formal and advanced Modified Integration Method (MIM), originally developed by Smets (1957). Such procedure consists of assuming that analytical solutions have the following form:

$$
n\left(t\right)=\int_{\Omega}{\bar{n}\left(s\right)e^{st}ds},\ \ C_i\left(t\right)=\int_{\Omega}{{\bar{C}}_i\left(s\right)e^{st}ds}$$

## Integral Representation of $$n(t)$$

$$
\begin{aligned}
n(t)=\;&
K_{1}\,e^{-\lambda t}
\int_0^\infty e^{-y^{2}/2+\mathcal{E}_{3}(t)\,y}\, y^{\lambda\beta/a}\,dy
\\[4pt]
&+K_{2}\,e^{-\lambda t}
\int_0^\infty e^{-y^{2}/2-\mathcal{E}_{3}(t)\,y}\,y^{\lambda\beta/a}\,dy
\\[4pt]
&+B_{p,1}
\int_0^\infty 
\exp\!\left[
-\frac{1}{a}\left(\frac{\Lambda}{2}s^{2}+(\beta-b-at)s\right)
\right]
\,(s+\lambda)^{\lambda\beta/a}\,ds.
\end{aligned}
$$

For completeness, we recall that the auxiliary function $$\mathcal{E}_3(t)$$ arises directly from the analytical reduction of the point kinetics equation with a linear reactivity. It is defined as

$$
\mathcal{E}_3(t)
=
\frac{1}{\sqrt{2a\Lambda}}
\left[
at + b - \beta + \frac{\Lambda\lambda}{2}
\right],
$$

which is obtained after completing the square in the exponent of the transformed equation for $n(t)$. The function $B_{p,1}$ appearing in the third integral term of $$n(t)$$ follows analogously from the same reduction procedure and depends linearly on the parameters $$a$$, $$b$$, $$\beta$$, $$\Lambda$$, $$\lambda$$ and $$q$$, as detailed in the associated manuscript.

---

## Delayed Neutron Precursor Concentration

Once $$n(t)$$ is known, the concentration of delayed neutron precursors is obtained from

$$
C(t)
= C(0)\,e^{-\lambda t}
+ \frac{\beta}{\Lambda}\,e^{-\lambda t}
\int_0^t e^{\lambda\tau}\,n(\tau)\,d\tau.
$$

The numerical implementation of this integral is described in the **Codes** section.

---

## RK4 Reference System

For validation purposes, a high-precision RK4 solver is applied directly to the NPKE system:

$$
\frac{dn}{dt}=\frac{\rho(t)-\beta}{\Lambda}\,n(t)+\lambda\,C(t)+q,
\qquad
\frac{dC}{dt}=\frac{\beta}{\Lambda}\,n(t)-\lambda\,C(t),
$$

with the same ramp reactivity $$\rho(t)=at+b$$.  

Further details and numerical parameters are provided in the **Validation** page.
