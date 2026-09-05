---
layout: default
title: MATLAB Codes
math: true
---

## Overview of the MATLAB Scripts

The repository contains five main MATLAB scripts. 
They implement the analytical power-series solution of the Bateman equations using the Mittag-Leffler function. Specifically, two of these scripts are dedicated to the proposed analytical solution, differing fundamentally in their computational strategy: one uses a cache file where the non-negative Diophantine solutions $k_0+k_1+...+k_{p-1}=m$ are precomputed, while the other computes them on the fly during execution.

## ⚠️ Prerequisite: Mittag-Leffler Evaluation Code

Please note that both of the proposed codes described, based on the power-series solution require the `ml.m` file to be located in the same directory for proper execution. This file contains the highly robust numerical implementation of the Mittag-Leffler function developed by Roberto Garrappa. You can download it directly from the MATLAB Central File Exchange:

👉 <a href="https://www.mathworks.com/matlabcentral/fileexchange/48154-the-mittag-leffler-function" target="_blank" rel="noopener noreferrer">The Mittag-Leffler function by Roberto Garrappa</a>

Disclaimer: The `ml.m` script is the intellectual property of its author, Dr. Roberto Garrappa. When downloading and utilizing this code, please ensure you carefully read, adhere to, and respect the specific licensing terms and conditions established by Dr. Garrappa.
---

# 1. Power_series_mittag.m

<div style="padding:8px; border-left:4px solid #3c6e71; margin-bottom:10px; background-color:#f9f9f9;">
  <a href="https://github.com/Cruz-Lopez-Carlos-Antonio/Power-Series-of-Bateman-Equations/blob/main/Power_series_mittag.m" 
     target="_blank" style="font-size:16px; color:#22577a; font-weight:bold; text-decoration:none;">
     👉 Click here to view the code in a new tab
  </a>
</div>

This script implements the optimized analytical solution using **MATLAB** and relies on Roberto Garrappa's algorithm for the Mittag-Leffler function evaluation.  
This version uses the precomputed cache file to load the non-negative Diophantine solutions, avoiding redundant combinatorial calculations. 

**Inputs:**  
It receives physical and temporal parameters such as the half-lives vector (`half_lives`), the initial concentration (`x10`), and the evaluation time grid (`Time_vector`), alongside the preloaded `DioCache`.

<div style="background:#f4f4f4; border:1px solid #ddd; border-left:4px solid #4a90e2; border-radius:4px; padding:10px; margin-bottom:15px; overflow-x:auto;">
<pre style="margin: 0; background: transparent; border: none;"><code>% Load precomputed Diophantine solutions
load('dio_cache_p14_m5.mat','DioCache');

% Example of physical parameters
half_lives = [2,2,3,3,3,4,4,4,4];
lambda = log(2) ./ half_lives;
Time_vector = [1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100];
x10 = 6.023e23;</code></pre>
</div>

The default data belongs to the scheme proposed by Dreher, which considers the following linear chain:

$$X_1(t) \to X_2(t) \to X_3(t) \to \cdots \to X_9(t),$$

where:
$$\lambda_1 = \lambda_2, \quad \lambda_3 = \lambda_4 = \lambda_5, \quad \text{and} \quad \lambda_6 = \cdots = \lambda_9,$$

That is, a linear chain of nine nuclides, some of which share identical decay constants, resulting in only three distinct eigenvalues. Specifically, the assigned half-lives are $2$, $3$, and $4$ seconds (defining the decay constants as $\lambda_i = \ln(2)/T_{1/2,i}$), with multiplicities of 2, 3, and 4, respectively.

**Outputs:**  
The script outputs the computed concentration $X_n(t)$ corresponding to the final isotope in the decay chain and automatically exports these results into a text file named `Bateman_superposition_results_optimized.txt` for further plotting or numerical analysis.

---

## 2. `Power_series_mittag_autonomus.m`
<div style="padding:8px; border-left:4px solid #3c6e71; margin-bottom:10px; background-color:#f9f9f9;">
  <a href="https://github.com/Cruz-Lopez-Carlos-Antonio/Power-Series-of-Bateman-Equations/blob/main/Power_series_mittag_autonomus.m" 
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
