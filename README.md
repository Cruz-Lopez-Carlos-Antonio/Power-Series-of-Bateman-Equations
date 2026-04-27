# A Novel Power-series Solution of the Bateman Equations using the Mittag-Leffler Function

## Overview of the Repository

The present repository contains the **MATLAB** codes associated with a novel power-series formulation for the Bateman equations based on the two-parameter Mittag-Leffler function $E_{\alpha,\beta}(x)$ and recent advances in fractional calculus. 

The Bateman equations describe the evolution of nuclide concentrations under radioactive decay and transmutation. Although their classical analytical solution is well known, its direct numerical evaluation can suffer from catastrophic cancellation in certain regimes, particularly when the product $\lambda t$ is small or when two decay constants are very close in value. This new representation departs significantly from classical expansions, exhibiting a more compact structure that avoids deeply nested summations and facilitates its practical implementation.

The repository also includes a set of computational strategies (such as recursive evaluation, time discretization, and a superposition-based stepping scheme) to ensure numerical stability over extended time intervals. These optimized strategies significantly lower the execution time, achieving reductions of over 90% with respect to the original, non-optimized version.

## Theoretical Framework

The core mathematical foundation of this computational implementation is established by the analytical solution:

$$\begin{aligned} x_n(t) &= x_1(0) \left( \prod_{k=1}^{n-1}\lambda_k \right) \sum_{m=0}^{\infty} (-1)^m \sum_{\substack{k_0+k_1+\cdots+k_{n-2}=m \\ k_0,k_1,\dots,k_{n-2}\geq 0}} \frac{1}{k_0!\,k_1!\cdots k_{n-2}!} \\[6pt] &\quad\times \left(\prod_{i=0}^{n-2} a_{n-i}^{k_i}\right) t^{\beta-1} E_{n,\beta-nm}^{(m)} \left( -a_nt^n \right) \end{aligned}$$

with the theoretical coefficients defined as:

$$a_0=1, \ \ a_k = \sum_{i_1=1}^{\,n-k+1} \sum_{i_2=i_1+1}^{\,n-k+2} \cdots \sum_{i_k=i_{k-1}+1}^{\,n} \prod_{j=1}^{k}\lambda_{i_j}$$

## Computational Implementation

For programmatic evaluation and improved computational performance, the analytical equation is structured into the following optimized form:

$$x_n(t) = \frac{x_1(0)}{a_0} \left( \prod_{j=1}^{n-1}\lambda_j \right) \sum_{\mathbf{k}\in\mathbb{N}_0^{\,n-1}} (-1)^{|\mathbf{k}|} \frac{\mathbf{c}^{\mathbf{k}}}{\mathbf{k}!} \,t^{\gamma(\mathbf{k})} E_{1,\beta(\mathbf{k})}^{(|\mathbf{k}|)} \bigl(-\mu t\bigr)$$

where $\mathbf{k}=(k_0,k_1,\dots,k_{n-2})\in\mathbb{N}_0^{\,n-1}$, with:

$$|\mathbf{k}|:=\sum_{i=0}^{n-2}k_i, \qquad \mathbf{k}!:=\prod_{i=0}^{n-2}k_i!, \qquad \mathbf{c}^{\mathbf{k}}:=\prod_{i=0}^{n-2}c_i^{k_i}$$

where $c_i=\frac{a_{n-i}}{a_0}$ for $i=0,\ldots,n-2$ and $\mu=\frac{a_1}{a_0}$. The parameters regulating the time exponents and the orders of the Mittag-Leffler function are defined by:

$$\beta(\mathbf{k}) = n+\sum_{i=0}^{n-2}\bigl((n-1)-i\bigr)k_i$$

$$\gamma(\mathbf{k}) = n-1+\sum_{i=0}^{n-2}(n-i)\,k_i$$

The scripts are implemented in **MATLAB**, which was selected due to the availability of one of the most efficient numerical implementations of the two-parameter Mittag-Leffler function developed within this framework.

## Authors and Financial Support

**Authors:**
* Carlos-Antonio Cruz-López (Universidad Autónoma Metropolitana)
* Marc Jornet (Escuela Superior de Ingeniería y Tecnología, Universidad Internacional de La Rioja)
* Gilberto Espinosa-Paredes (Universidad Autónoma Metropolitana)

**Financial Support**
The authors Carlos Antonio Cruz López and Gilberto Espinosa Paredes gratefully acknowledge the financial support received from the Secretaría de Ciencia, Humanidades, Tecnología e Innovación (SECIHTI, formerly known as CONAHCYT), through the program *Estancias Posdoctorales por México, 2022*, under the project entitled: *Desarrollo de modelos fenomenológicos energéticos de orden fraccional, para la optimización y simulación en reactores nucleares de potencia*, as well as the financial support from the Basic Science and Frontier Project 2023-2024, with the reference CBF-2023-2024-2023, also belonging to SECIHTI-CONAHCYT.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. You are free to use, modify, and distribute this software for academic and commercial purposes, provided that appropriate credit is given to the original authors.
