---
layout: default
title: Main Equations
math: true
---

# The Bateman Equations
The codes developed in the present repository solves the following system of differential equations:

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\frac{dx_i(t)}{dt}
=\lambda_{i-1} x_{i-1}(t)-\lambda_i x_i(t), \qquad 1 \le i \le n,
$$
</div>

which models the linear following decay chain:

$$
x_1(t) \xrightarrow{\lambda_1} 
x_2(t) \xrightarrow{\lambda_2} 
\cdots 
\xrightarrow{\lambda_{n-1}} 
x_n(t) \xrightarrow{\lambda_n}\cdots,$$

Subject to the following initial conditions:
$x_i(0)=0$, with $i=2,3,\dots,n$.

---

## Analytical Solution using the Mittag-Leffler Function

The analytical solution obtained for the aforementioned system is expressed as follows:

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
%
$$
</div>

where the coefficients $a_i$ depends on the lambda parameters as follows:

$$
a_0=1, \ \ a_k =
\sum_{i_1=1}^{\,n-k+1}
\sum_{i_2=i_1+1}^{\,n-k+2}
\cdots \sum_{i_k=i_{k-1}+1}^{\,n}
\prod_{j=1}^{k}\lambda_{i_j}.
$$

and where the Mittag-Leffler function, $E_{\alpha,\beta}(z)$, is defined as

$$
E_{\alpha,\beta}(z)=
\sum_{k=0}^{\infty}
\frac{z^k}{\Gamma(\alpha k+\beta)},
\qquad
\Re(\alpha)>0,
\quad z,\beta\in\mathbb{C}.
$$

For the case of the derivatives, the following formulation is used in the work

$$
E_{\alpha,\beta}^{(m)}(z)
=m!\,E_{\alpha,\beta+\alpha m}^{\,m+1}(z).
$$

where $E_{\alpha,\beta}^{\gamma}(z)$ is the Prabhakar function. For the more details, see the related entro in this repository.
Finally, this expression can be written in a more compact vectorial form, which is very convenient to the computational implementation:
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
---

## Validation

For further details on the high-precision numerical tests and the evaluation of relative errors for these solutions, please see the **Validation** page.
