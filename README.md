# A Novel Power-series Solution of the Bateman Equations using the Mittag-Leffler Function

## Overview of the Repository

The present repository contains the **MATLAB** codes associated with a novel power-series formulation for the Bateman equations based on the two-parameter Mittag-Leffler function $E_{\alpha,\beta}(x)$ and recent advances in fractional calculus. 

The Bateman equations describe the evolution of nuclide concentrations under radioactive decay and transmutation. Although their classical analytical solution is well known, its direct numerical evaluation can suffer from catastrophic cancellation in certain regimes, particularly when the product $\lambda t$ is small or when two decay constants are very close in value. This new representation departs significantly from classical expansions, exhibiting a more compact structure that avoids deeply nested summations and facilitates its practical implementation.

The repository also includes a set of computational strategies (such as recursive evaluation, time discretization, and a superposition-based stepping scheme) to ensure numerical stability over extended time intervals. These optimized strategies significantly lower the execution time, achieving reductions of over 90% with respect to the original, non-optimized version.


## Authors and Financial Support

**Authors:**
* Carlos-Antonio Cruz-López (Universidad Autónoma Metropolitana)
* Marc Jornet (Escuela Superior de Ingeniería y Tecnología, Universidad Internacional de La Rioja)
* Gilberto Espinosa-Paredes (Universidad Autónoma Metropolitana)

**Financial Support**
The authors Carlos Antonio Cruz López and Gilberto Espinosa Paredes gratefully acknowledge the financial support received from the Secretaría de Ciencia, Humanidades, Tecnología e Innovación (SECIHTI, formerly known as CONAHCYT), through the program *Estancias Posdoctorales por México, 2022*, under the project entitled: *Desarrollo de modelos fenomenológicos energéticos de orden fraccional, para la optimización y simulación en reactores nucleares de potencia*, as well as the financial support from the Basic Science and Frontier Project 2023-2024, with the reference CBF-2023-2024-2023, also belonging to SECIHTI-CONAHCYT.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. You are free to use, modify, and distribute this software for academic and commercial purposes, provided that appropriate credit is given to the original authors.

## Web of the Repository

<div style="background:#f1f7ff; padding:12px 15px; border-left:4px solid #4a90e2; border-radius:8px; margin-top:20px;">
  <strong>📘 Full Documentation</strong><br>
  The complete project website (opens in a new page) is available here:<br><br>
  👉 <a href="https://cruz-lopez-carlos-antonio.github.io/Ramp_analytical_solution/"
        target="_blank" rel="noopener noreferrer">
        https://cruz-lopez-carlos-antonio.github.io/Ramp_analytical_solution/
      </a>
</div>
