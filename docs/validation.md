---
layout: default
title: Numerical Validation
math: true
---

## Purpose

The goal of this section is to document how the proposed analytical solutions are validated against exact reference standards to ensure their accuracy and reliability.

---

## Case 1. Dreher's scheme with short time values

To verify the proposed power-series solution, the scheme reported by Dreher is employed as a benchmark for systems with repeated roots[cite: 1]. This test problem consists of a linear chain of nine nuclides[cite: 1]:

$$
X_1(t) \to X_2(t) \to X_3(t) \to \cdots \to X_9(t)
$$

where the decay constants satisfy the following conditions[cite: 1]:

$$
\lambda_1 = \lambda_2, \quad \lambda_3 = \lambda_4 = \lambda_5, \quad \text{and} \quad \lambda_6 = \cdots = \lambda_9
$$

Specifically, the assigned half-lives are 2, 3, and 4 seconds, with multiplicities of 2, 3, and 4, respectively[cite: 1]. Following Dreher's original setup, the initial concentration for the first member of the chain is established as $x_1(0) = 6.023 \times 10^{23}$ atoms, with all subsequent nuclides being initially absent[cite: 1]. 

The results obtained from the Mittag-Leffler function-based formulation are compared against an exact analytical reference previously derived via the Laplace transform method, confirming excellent agreement without detailing the asymptotic decay behavior extensively here[cite: 1].

👉 <a href="https://github.com/Cruz-Lopez-Carlos-Antonio/Power-Series-of-Bateman-Equations/blob/main/Verification/Complete_tables_verification.txt" target="_blank" rel="noopener noreferrer">Complete verification tables for Dreher's scheme</a>
