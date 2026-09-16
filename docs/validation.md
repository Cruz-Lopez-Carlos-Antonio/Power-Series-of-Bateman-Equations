---
layout: default
title: Numerical Validation
math: true
---

## Purpose

The goal of this section is to document how the proposed analytical solutions are validated against exact reference standards to ensure their accuracy and reliability.

---

## Case 1. Dreher's scheme with short time values

To verify the proposed power-series solution, the scheme reported by Dreher is employed as a benchmark for systems with repeated roots. This test problem consists of a linear chain of nine nuclides:

$$
X_1(t) \to X_2(t) \to X_3(t) \to \cdots \to X_9(t)
$$

where the decay constants satisfy the following conditions:

$$
\lambda_1 = \lambda_2, \quad \lambda_3 = \lambda_4 = \lambda_5, \quad \text{and} \quad \lambda_6 = \cdots = \lambda_9
$$

Specifically, the assigned half-lives are 2, 3, and 4 seconds, with multiplicities of 2, 3, and 4, respectively. Following Dreher's original setup, the initial concentration for the first member of the chain is established as $x_1(0) = 6.023 \times 10^{23}$ atoms, with all subsequent nuclides being initially absent. 

The results obtained from the Mittag-Leffler function-based formulation are compared against an exact analytical reference previously derived via the Laplace transform method, confirming excellent agreement without detailing the asymptotic decay behavior extensively here[cite: 1].

👉 <a href="https://github.com/Cruz-Lopez-Carlos-Antonio/Power-Series-of-Bateman-Equations/blob/main/Verification/Complete_tables_verification.txt" target="_blank" rel="noopener noreferrer">Complete verification tables for Dreher's scheme</a>

### Results for $X_1(t)$

The following figure illustrates the time evolution of the nuclide concentrations.

*(Nota para ti: compila tu Figura_4.txt a formato PNG/JPG, súbela a tu repositorio y reemplaza el enlace de abajo)*
![Time evolution of the nuclide concentrations](ruta/a/tu_figura_4.png)

Below is a condensed view of the verification results for the first nuclide, $X_1(t)$, at selected time steps. 

| $t$ (s) | Reference Solution (Laplace) | Series Solution (Mittag-Leffler) | APE (%) |
| :--- | :--- | :--- | :--- |
| 10 | 18821875000000001428776.46 | 1.8821874999999985e+22 | 8.7286e-14 |
| 20 | 588183593750000112850.53 | 5.8818359374999906e+20 | 1.7900e-13 |
| 50 | 17949938774108896.41 | 1.7949938774108814e+16 | 4.5909e-13 |
| 100 | 534949862.19 | 5.3494986218538117e+08 | 9.0781e-13 |

<details>
<summary>👉 Click to view full data</summary>

<br>

| $t$ (s) | Reference Solution (Laplace) | Series Solution (Mittag-Leffler) | APE (%) |
| :--- | :--- | :--- | :--- |
| 1 | 425890414308657561831457.99 | 4.2589041430865752e+23 | 9.8221e-15 |
| 2 | 301149999999999994925185.48 | 3.0114999999999992e+23 | 2.4880e-14 |
| 3 | 212945207154328785854028.04 | 2.1294520715432873e+23 | 2.6229e-14 |
| 4 | 150575000000000000954497.48 | 1.5057499999999994e+23 | 4.0481e-14 |
| 5 | 106472603577164395396163.54 | 1.0647260357716435e+23 | 4.2636e-14 |
| 6 | 75287500000000002223201.11 | 7.5287499999999963e+22 | 5.2098e-14 |
| 7 | 53236301788582198932656.53 | 5.3236301788582165e+22 | 6.3740e-14 |
| 8 | 37643750000000001984576.74 | 3.7643749999999973e+22 | 7.6997e-14 |
| 9 | 26618150894291100083615.65 | 2.6618150894291078e+22 | 8.2964e-14 |
| 10 | 18821875000000001428776.46 | 1.8821874999999985e+22 | 8.7286e-14 |
| 20 | 588183593750000112850.53 | 5.8818359374999906e+20 | 1.7900e-13 |
| 30 | 18380737304687505657.87 | 1.8380737304687456e+19 | 2.7016e-13 |
| 40 | 574398040771484618.41 | 5.7439804077148256e+17 | 3.5836e-13 |
| 50 | 17949938774108896.41 | 1.7949938774108814e+16 | 4.5909e-13 |
| 60 | 560935586690903.08 | 5.6093558669090000e+14 | 5.4868e-13 |
| 70 | 17529237084090.72 | 1.7529237084090611e+13 | 6.4014e-13 |
| 80 | 547788658877.84 | 5.4778865887783118e+11 | 7.2727e-13 |
| 90 | 17118395589.93 | 1.7118395589932211e+10 | 8.1700e-13 |
| 100 | 534949862.19 | 5.3494986218538117e+08 | 9.0781e-13 |

</details>
