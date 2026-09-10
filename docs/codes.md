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

$$
E_{\alpha,\beta}^{(m)}(z)
=m!\,E_{\alpha,\beta+\alpha m}^{\,m+1}(z).
$$
</div>

Such module was developed by Roberto Garrappa. You can download it directly from the MATLAB Central File Exchange:

👉 <a href="https://www.mathworks.com/matlabcentral/fileexchange/48154-the-mittag-leffler-function" target="_blank" rel="noopener noreferrer">The Mittag-Leffler function by Roberto Garrappa</a>

Disclaimer: The `ml.m` script is the intellectual property of its author, Dr. Roberto Garrappa. When downloading and utilizing this code, please ensure you carefully read, adhere to, and respect the specific licensing terms and conditions established by Dr. Garrappa.

<hr style="border: none; border-top: 3px solid #3b5998; margin: 1.5rem 0;">

## 1. Power_series_mittag.m

<div style="padding:8px; border-left:4px solid #3c6e71; margin-bottom:10px; background-color:#f9f9f9;">
  <a href="https://github.com/Cruz-Lopez-Carlos-Antonio/Power-Series-of-Bateman-Equations/blob/main/Power_series_mittag.m" 
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

which can be written in a more compact vectorial form, highly convenient for computational implementation:

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

```matlab
% Load precomputed Diophantine solutions
load('dio_cache_p14_m5.mat','DioCache');

% Example of physical parameters
half_lives = [2,2,3,3,3,4,4,4,4];
lambda = log(2) ./ half_lives;
Time_vector = [1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100];
x10 = 6.023e23;
