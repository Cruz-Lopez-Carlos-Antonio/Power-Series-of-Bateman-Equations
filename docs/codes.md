---
layout: default
title: MATLAB Codes
math: true
---

## Overview of the MATLAB Scripts

The repository contains five main MATLAB scripts. Two of them implement the analytical power-series solution of the Bateman equations using the Mittag-Leffler function, another two contain the analytical solution based on the Laplace method developed previously [1, 2], and the final one implements the solution developed by Velhinho et al. [3].

The first two codes are based on the power-series solution and differ fundamentally in their computational strategy: one uses a cache file where the non-negative Diophantine solutions $k_0+k_1+\dots+k_{p-1}=m$ are precomputed, while the other computes them on the fly during execution.

The two codes implementing the Laplace transform use the same procedure, but with different digit precisions: one utilizes the Variable-Precision Arithmetic (VPA) approach, and the other uses standard MATLAB precision.

<hr style="border: none; border-top: 3px solid #3b5998; margin: 1.5rem 0;">

## ⚠️ Prerequisite: Mittag-Leffler Evaluation Code

Please note that both of the proposed codes, based on the power-series solution require the `ml.m` file to be located in the same directory for proper execution. This file contains the highly robust numerical implementation of the Mittag-Leffler function.

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
E_{\alpha,\beta}(z)=
\sum_{k=0}^{\infty}
\frac{z^k}{\Gamma(\alpha k+\beta)},
\qquad
\Re(\alpha)>0,
\quad z,\beta\in\mathbb{C}.
$$
</div>
  
For the case of the derivatives, the following formulation is used in the work:

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
E_{\alpha,\beta}^{(m)}(z)
=m!\,E_{\alpha,\beta+\alpha m}^{\,m+1}(z).
$$
</div>


where $E_{\alpha,\beta}^{\gamma}(z)$ is the Prabhakar function, which is defined as:

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
E_{\alpha,\beta}^{\,\gamma}(z)
=\sum_{k=0}^{\infty}
\frac{\Gamma(\gamma+k)}
{\Gamma(\gamma)\,k!\,\Gamma(\alpha k+\beta)}
\,z^k.
$$
</div>

Such module was developed by Roberto Garrappa. You can download it directly from the MATLAB Central File Exchange:

👉 <a href="[https://www.mathworks.com/matlabcentral/fileexchange/48154-the-mittag-leffler-function](https://www.mathworks.com/matlabcentral/fileexchange/48154-the-mittag-leffler-function)" target="_blank" rel="noopener noreferrer">The Mittag-Leffler function by Roberto Garrappa</a>

Disclaimer: The `ml.m` script is the intellectual property of its author, Dr. Roberto Garrappa. When downloading and utilizing this code, please ensure you carefully read, adhere to, and respect the specific licensing terms and conditions established by Dr. Garrappa.

<hr style="border: none; border-top: 3px solid #3b5998; margin: 1.5rem 0;">

## 1. Power_series_mittag.m

<div style="padding:8px; border-left:4px solid #3c6e71; margin-bottom:10px; background-color:#f9f9f9;">
  <a href="[https://github.com/Cruz-Lopez-Carlos-Antonio/Power-Series-of-Bateman-Equations/blob/main/Power_series_mittag.m](https://github.com/Cruz-Lopez-Carlos-Antonio/Power-Series-of-Bateman-Equations/blob/main/Power_series_mittag.m)" 
     target="_blank" style="font-size:16px; color:#22577a; font-weight:bold; text-decoration:none;">
     👉 Click here to view the code in a new tab
  </a>
</div>

This script implements the optimized analytical solution using **MATLAB**.  
This version uses the precomputed cache file to load the non-negative Diophantine solutions, avoiding redundant combinatorial calculations. 

The script evaluates the following analytical expression:

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\begin{equation}
\begin{aligned}
x_n(t)&=x_1(0)
\left(\prod_{k=1}^{n-1}\lambda_k
\right)
\frac{1}{a_0}\sum_{m=0}^{\infty}
(-1)^m\sum_{\substack{
k_0+k_1+\cdots+k_{n-2}=m \\
k_0,k_1,\ldots,k_{n-2}\geq 0
}}\frac{1}{k_0!\,k_1!\cdots k_{n-2}!}
\\[6pt]
&\quad\times
\left(\prod_{i=0}^{n-2}\left(\frac{a_{n-i}}{a_0}
\right)^{k_i}\right)
t^{m+\beta-1}E_{1,\beta}^{(m)}
\left(-\frac{a_1}{a_0}t\right),
\end{aligned}
\end{equation}
$$
</div>

where the coefficients $a_i$ depends on the lambda parameters as follows:

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
a_0=1, \ \ a_k =
\sum_{i_1=1}^{\,n-k+1}
\sum_{i_2=i_1+1}^{\,n-k+2}
\cdots \sum_{i_k=i_{k-1}+1}^{\,n}
\prod_{j=1}^{k}\lambda_{i_j}.
$$
</div>

Finally, this expression can be written in a more compact vectorial form, which is very convenient for computational implementation:

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
x_n(t)=
\frac{x_1(0)}{a_0}
\left(\prod_{j=1}^{n-1}\lambda_j\right)
\sum_{\mathbf{k}\in\mathbb{N}_0^{\,n-1}}
(-1)^{|\mathbf{k}|}
\frac{\mathbf{c}^{\mathbf{k}}}{\mathbf{k}!}
\,t^{\gamma(\mathbf{k})}
E_{1,\beta(\mathbf{k})}^{(|\mathbf{k}|)}
\bigl(-\mu t\bigr).
$$
</div>

**Inputs:**  
It receives physical and temporal parameters such as the half-lives vector (`half_lives`), the initial concentration (`x10`), and the evaluation time grid (`Time_vector`), alongside the preloaded `DioCache`.

<div style="background:#f4f4f4; border:1px solid #ddd; border-left:4px solid #4a90e2; border-radius:4px; padding:10px; margin-bottom:15px; overflow-x:auto;">
<pre style="margin: 0; background: transparent; border: none; font-family: monospace; color: #333;">
% Load precomputed Diophantine solutions
load('dio_cache_p14_m5.mat','DioCache');

% Example of physical parameters
half_lives = [2,2,3,3,3,4,4,4,4];
lambda = log(2) ./ half_lives;
Time_vector = [1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100];
x10 = 6.023e23;
</pre>
</div>

The default data belongs to the scheme proposed by Dreher, which considers the following linear chain:

$$X_1(t) \to X_2(t) \to X_3(t) \to \cdots \to X_9(t),$$

where:
$$\lambda_1 = \lambda_2, \quad \lambda_3 = \lambda_4 = \lambda_5, \quad \text{and} \quad \lambda_6 = \cdots = \lambda_9,$$

That is, a linear chain of nine nuclides, some of which share identical decay constants, resulting in only three distinct eigenvalues. Specifically, the assigned half-lives are $2$, $3$, and $4$ seconds (defining the decay constants as $\lambda_i = \ln(2)/T_{1/2,i}$), with multiplicities of 2, 3, and 4, respectively.

**Outputs:**  
The script outputs the computed concentration $X_n(t)$ corresponding to the final isotope in the decay chain and automatically exports these results into a text file named `Bateman_superposition_results_optimized.txt` for further plotting or numerical analysis.

<hr style="border: none; border-top: 3px solid #3b5998; margin: 1.5rem 0;">

## 2. Power_series_mittag_autonomus.m

<div style="padding:8px; border-left:4px solid #3c6e71; margin-bottom:10px; background-color:#f9f9f9;">
  <a href="[https://github.com/Cruz-Lopez-Carlos-Antonio/Power-Series-of-Bateman-Equations/blob/main/Power_series_mittag_autonomus.m](https://github.com/Cruz-Lopez-Carlos-Antonio/Power-Series-of-Bateman-Equations/blob/main/Power_series_mittag_autonomus.m)" 
     target="_blank" style="font-size:16px; color:#22577a; font-weight:bold; text-decoration:none;">
     👉 Click here to view the code in a new tab
  </a>
</div>

This script provides a standalone, or **autonomous**, version of the analytical solution given above.  
Unlike the optimized version, this script does *not* require an external cache file. Instead, it computes the combinations "on the fly" by generating the weak compositions for the Diophantine sum:

$$
k_0 + k_1 + \dots + k_{p-1} = m, \qquad k_i \ge 0
$$

This combinatorial generation is achieved using a "stars and bars" approach directly within the code. While mathematically identical to the cached version, it provides high flexibility when external `.mat` files are not desired or when the user requires studying linear chains with more than 15 nuclides. 

**Inputs & Outputs:**  
It receives the exact same inputs (`half_lives`, `x10`, `Time_vector`) but omits the dictionary loading procedure. It yields the exact same concentration array $X_n(t)$ and outputs a text file named `Bateman_superposition_results_optimized_autonomous.txt`.

<hr style="border: none; border-top: 3px solid #3b5998; margin: 1.5rem 0;">

## 3. Bateman_reference_double_precision.m

<div style="padding:8px; border-left:4px solid #3c6e71; margin-bottom:10px; background-color:#f9f9f9;">
  <a href="[https://github.com/Cruz-Lopez-Carlos-Antonio/Power-Series-of-Bateman-Equations/blob/main/Bateman_reference_double_precision.m](https://github.com/Cruz-Lopez-Carlos-Antonio/Power-Series-of-Bateman-Equations/blob/main/Bateman_reference_double_precision.m)" 
     target="_blank" style="font-size:16px; color:#22577a; font-weight:bold; text-decoration:none;">
     👉 Click here to view the code in a new tab
  </a>
</div>

This script solves the generalized Bateman equations by implementing the closed-form analytical solution for linear chains with repeated decay constants. This exact mathematical formulation was developed and simplified in the works by Cruz-López et al. [1, 2]. 

The implementation evaluates the following analytical expression, which avoids nested summations and simplifies the combinatorial structure using Frobenius-type Diophantine equations (Cauchy products):

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
X_N(t) = \frac{X_1(0)}{\lambda_N} \left( \prod_{k=1}^n \lambda_k^{\mu_k+1} \right) \sum_{i=1}^n \frac{\exp(-\lambda_i t)}{\prod_{\substack{j=1 \\ j \neq i}}^n (\lambda_j - \lambda_i)^{\mu_j+1}} \sum_{\ell=0}^{\mu_i} \frac{t^\ell}{\ell!} \chi_{i,\mu_i-\ell}
$$
</div>

where the coefficient $\chi_{i,r}$ is defined as:

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\chi_{i,r} = \sum_{h_1+\dots+h_{i-1}+h_{i+1}+\dots+h_n=r} \prod_{\substack{k=1 \\ k \neq i}}^n \binom{h_k+\mu_k}{\mu_k} \frac{1}{(\lambda_i-\lambda_k)^{h_k}}
$$
</div>

**Variables Description:**
*   $X_N(t)$: Concentration of the $N$-th (final) nuclide in the decay chain at evaluation time $t$.
*   $X_1(0)$: Initial concentration (number of atoms) of the first member of the chain.
*   $n$: Number of distinct decay constants in the linear chain.
*   $\lambda_i, \lambda_k$: The grouped distinct decay constants.
*   $\mu_k$: Multiplicity parameter, meaning a specific decay constant appears $\mu_k+1$ times in the linear chain.
*   $h_k$: Non-negative integers that satisfy the Diophantine sum condition for the restricted partitions.

**Inputs & Initial Conditions:**  
This script uses standard arrays for the physical parameters and the initial concentration of the first nuclide.

<div style="background:#f4f4f4; border:1px solid #ddd; border-left:4px solid #4a90e2; border-radius:4px; padding:10px; margin-bottom:15px; overflow-x:auto;">
<pre style="margin: 0; background: transparent; border: none; font-family: monospace; color: #333;">
% Physical parameters
Half_lifes = [2, 2, 3, 3, 3, 4];
D = log(2) ./ Half_lifes;

% Evaluation time grid
Time_vector = [0.001, 0.002, 0.003, 0.004, 0.005, ...
               0.006, 0.007, 0.008, 0.009, 0.010];

% Initial number of atoms in the first member of the chain
X0 = 6.023e23;
</pre>
</div>

**Numerical Note:**  
It receives the exact same initial parameters (`half_lives`, `x10`, `Time_vector`). However, it is important to note that this specific code intentionally utilizes standard IEEE MATLAB double precision. Consequently, for very small evaluation times, the output may exhibit a loss of numerical accuracy due to severe cancellation among large terms.

<hr style="border: none; border-top: 3px solid #3b5998; margin: 1.5rem 0;">

## References
1. Cruz-López, C.-A., Espinosa-Paredes, G., & François, J.-L. (2024). General solution of Bateman equations using Cauchy products and the Theory of Divided Differences. *Annals of Nuclear Energy, 207*, 110729 (p. 7).
2. Cruz-López, C.-A., Jornet, M., Espinosa-Paredes, G., & François, J.-L. (2026). On the generalized summation of series with rational coefficients. *Computer Physics Communications, 325*, 110198 (p. 4).
3. Velhinho, J., Fonseca, E., & Serôdio, R. (2023). General solutions to decay chain equations. *Computer Physics Communications, 283*, 108582.
