---
layout: post
title: "A Source-Controlled Nonlinear Kapitza Model: Exact Reduction, Pitchfork Bifurcation, and Physical Reconstruction"
date: 2026-07-23
categories: [PDE, Bifurcation, Heat Transfer, Numerical Analysis]
---

This note develops a revised one-dimensional nonlinear Kapitza transmission model in which the normalized intensity $\kappa$ of a positive volume heat source is the bifurcation parameter. The key change is physical as well as notational: $\kappa$ is not a bulk conductivity and is not an artificial sigmoid-slope parameter. It is the strength of the source acting in the right layer.

The two-layer Poisson problem can be reduced exactly to one scalar equilibrium equation,

$$
\boxed{x=\tanh(\kappa x+h)}.
$$

This reduction preserves enough structure to obtain the critical source strength, all symmetric branches, the imperfect-pitchfork folds, the local two-thirds scaling law, and the full piecewise temperature fields. It also produces a clean benchmark for matching Python and MATLAB implementations.

The complete supporting materials are available here:

- [Detailed English report (PDF)](/assets/kapitza-kappa/kapitza_kappa_report_en.pdf)
- [English presentation slides (PDF)](/assets/kapitza-kappa/kapitza_kappa_beamer_en.pdf)
- [Read the Python implementation with syntax highlighting](/assets/kapitza-kappa/code/viewer.html?lang=python)
- [Read the MATLAB implementation with syntax highlighting](/assets/kapitza-kappa/code/viewer.html?lang=matlab)
- [Download the Python source](/assets/kapitza-kappa/code/kapitza_kappa_model.py)
- [Download the MATLAB source](/assets/kapitza-kappa/code/kapitza_kappa_model.m)
- [Download the report LaTeX source](/assets/kapitza-kappa/code/kapitza_kappa_report_en.tex)
- [Download the Beamer LaTeX source](/assets/kapitza-kappa/code/kapitza_kappa_beamer_en.tex)

## 1. Why introduce a source-controlled model?

In the earlier reduced model, a parameter was introduced directly into an exponential conductance to control the slope of the sigmoid map. That construction is mathematically useful, but the control parameter is not inherited from the bulk heat equation.

The revised formulation starts from a physical source term. The bulk conductivities remain

$$
\lambda>0
\qquad\text{and}\qquad
\mu>0,
$$

while $\kappa>0$ controls a positive volume source supported in the right subdomain:

$$
f(z)=2\mu\kappa\,\chi_{\Omega_E}(z).
$$

After division by $\mu$, the normalized source amplitude is $2\kappa$. This distinction matters because several unrelated quantities can otherwise be denoted by similar symbols:

- $\lambda$ and $\mu$ are material conductivities;
- $\delta=\lambda/\mu$ is their ratio;
- $2\kappa$ is the normalized source amplitude;
- $\kappa$ is the bifurcation parameter;
- $h$ is a symmetry-breaking constitutive bias.

The PDE backbone follows the one-dimensional reduction in Section 7 of [Semilinear Elliptic Problems with NonLinear Contact Resistance in Composite Domains](https://hal.science/hal-05659943v1/document). The exponential contact law used below is a new constitutive choice built on that reduction; the cited work principally studies power-law conductances in its bifurcation examples.

## 2. The two-layer Poisson problem

Consider

$$
\Omega=[-1,1],
\qquad
\Omega_I=[-1,0],
\qquad
\Omega_E=[0,1].
$$

The exterior temperatures are homogeneous:

$$
u_I(-1)=0,
\qquad
u_E(1)=0.
$$

The left layer is source-free, while the right layer contains the positive source:

$$
-\lambda u_I''=0
\quad\text{in }\Omega_I,
\qquad
-\mu u_E''=2\mu\kappa
\quad\text{in }\Omega_E.
$$

Solving the bulk equations and imposing the endpoint conditions gives

$$
u_I(z)=a(z+1),
$$

$$
u_E(z)=b(z-1)+\kappa(1-z^2).
$$

The interface traces are therefore

$$
u_I(0)=a,
\qquad
u_E(0)=\kappa-b.
$$

The scalar coefficients $a$ and $b$ are not independent. Flux conservation yields

$$
\lambda a=\mu b.
$$

With

$$
\delta=\frac{\lambda}{\mu},
$$

we obtain

$$
b=\delta a.
$$

## 3. Exact elimination of the bulk fields

Use the jump convention

$$
s=[u]=u_E(0)-u_I(0).
$$

Substituting the two interface traces gives the exact relation

$$
\boxed{s=\kappa-(1+\delta)a.}
$$

If $\beta_{\mathrm{phys}}$ is the dimensional Kapitza conductance, define its normalization by

$$
\widehat{\beta}=\frac{\beta_{\mathrm{phys}}}{\mu}.
$$

The normalized common heat flux satisfies

$$
\widehat q=\delta a=\widehat{\beta}(s)s.
$$

Consequently, the admissible physical interval depends on the source strength:

$$
0<a<\gamma(\kappa),
\qquad
\gamma(\kappa)=\frac{\kappa}{1+\delta},
\qquad
0<s<\kappa.
$$

This moving interval is one of the main differences from a fixed-domain scalar model. Increasing $\kappa$ changes both the nonlinear equilibrium equation and the scale on which the physical interface trace $a$ lives.

## 4. Exponential nonlinear contact conductance

The constitutive law is chosen as

$$
\boxed{
\widehat{\beta}_{\kappa,h}(s)
=
\frac{\delta}{1+\delta}
\exp\!\left[2(\kappa-2s)+2h\right].
}
$$

For fixed $\kappa$ and $h$, it is positive, decreasing, and convex in the jump:

$$
\partial_s\widehat{\beta}_{\kappa,h}
=-4\widehat{\beta}_{\kappa,h}<0,
$$

$$
\partial_{ss}\widehat{\beta}_{\kappa,h}
=16\widehat{\beta}_{\kappa,h}>0.
$$

The bias $h$ rescales the conductance multiplicatively:

$$
\widehat{\beta}_{\kappa,h}(s)
=e^{2h}\widehat{\beta}_{\kappa,0}(s).
$$

Solving the interface balance for $a$ produces the physical fixed-point map

$$
\boxed{
a=F_{\kappa,h}(a)
=
\frac{\kappa\widehat{\beta}_{\kappa,h}(s)}
{\delta+(1+\delta)\widehat{\beta}_{\kappa,h}(s)},
\qquad
s=\kappa-(1+\delta)a.
}
$$

The following overview shows four equivalent descriptions of the same problem: fixed-point intersections, the conductance law, heat-flux intersections, and the spatial source.

![Four equivalent views of the source-controlled nonlinear Kapitza model](/assets/kapitza-kappa/model_summary.png)

## 5. Exact normalization

Normalize the interface trace relative to its moving physical interval:

$$
x=\frac{2a}{\gamma(\kappa)}-1
=\frac{2(1+\delta)a}{\kappa}-1.
$$

Using $s=\kappa-(1+\delta)a$, this is equivalently

$$
x=1-\frac{2s}{\kappa},
\qquad -1<x<1.
$$

The fixed-point map becomes

$$
F_{\kappa,h}(a)
=
\frac{\gamma(\kappa)}{2}
\left[1+\tanh(\kappa x+h)\right].
$$

Therefore the complete steady transmission problem is exactly equivalent to

$$
\boxed{x=\tanh(\kappa x+h)}.
$$

The material ratio $\delta$ disappears from the normalized root equation, but it remains essential when converting $x$ back to $a$, $b$, $s$, the heat flux, and the temperature fields.

## 6. Symmetric problem and the critical source strength

Set $h=0$ and define

$$
G(x)=\tanh(\kappa x)-x.
$$

The central root $x=0$ always exists, and

$$
G'(0)=\kappa-1.
$$

It follows that

$$
\boxed{
\begin{array}{ll}
0<\kappa\leq1: & \text{one steady state},\\
\kappa>1: & \text{three steady states}.
\end{array}}
$$

The critical source strength is

$$
\boxed{\kappa_c=1.}
$$

The central physical branch is

$$
x=0,
\qquad
a_c(\kappa)=\frac{\kappa}{2(1+\delta)},
\qquad
s_c(\kappa)=\frac{\kappa}{2}.
$$

The two nonzero branches admit the exact parameterization

$$
\kappa(x)=\frac{\operatorname{artanh}(x)}{x},
$$

$$
a(x)=\frac{\kappa(x)(1+x)}{2(1+\delta)},
\qquad 0<|x|<1.
$$

Near $x=0$,

$$
\frac{\operatorname{artanh}(x)}{x}
=1+\frac{x^2}{3}+O(x^4),
$$

so

$$
x\sim\pm\sqrt{3(\kappa-1)}.
$$

This is a supercritical pitchfork in normalized coordinates.

![Fixed-point intersections below, at, and above the critical source strength](/assets/kapitza-kappa/fixed_point_regimes.png)

## 7. Normalized and physical bifurcation geometry

The standard pitchfork appears in the $(\kappa,x)$ plane. In the physical $(\kappa,a)$ plane, however, the central branch is not horizontal because the admissible interval grows with $\kappa$:

$$
a_c(\kappa)=\frac{\kappa}{2(1+\delta)}.
$$

Thus the physical branches unfold around a moving center line. This is not a plotting artifact; it is a direct consequence of using the source strength as the control parameter.

![The standard normalized pitchfork and its representation on the moving physical interval](/assets/kapitza-kappa/symmetric_bifurcation.png)

## 8. Picard iteration stability

For the scalar iteration

$$
x_{n+1}=\tanh(\kappa x_n+h),
$$

the fixed-point multiplier is

$$
m_{\mathrm{FP}}
=
\kappa\operatorname{sech}^2(\kappa x+h).
$$

At a fixed point, $\tanh(\kappa x+h)=x$, hence

$$
\boxed{m_{\mathrm{FP}}=\kappa(1-x^2).}
$$

The notation $m_{\mathrm{FP}}$ avoids confusion with the conductivity $\lambda$. A root is locally attracting for Picard iteration when

$$
|m_{\mathrm{FP}}|<1.
$$

Therefore:

- the central branch is attracting for $\kappa<1$;
- it is neutral at $\kappa=1$;
- it is repelling for $\kappa>1$;
- the two outer symmetric branches are attracting.

The cobweb plots expose convergence below threshold, critical slowing at threshold, and branch selection above threshold.

![Cobweb paths across the subcritical, critical, and supercritical regimes](/assets/kapitza-kappa/cobweb_regimes.png)

This classification concerns the numerical fixed-point iteration. It is not a physical stability result for a time-dependent heat equation, which would require a separate dynamical model and linearization.

## 9. Bias and imperfect pitchfork unfolding

For a fixed nonzero bias $h$, continuation in the physical source parameter is parameterized by

$$
\boxed{
\kappa(x;h)
=
\frac{\operatorname{artanh}(x)-h}{x}.
}
$$

A positive bias favors one outer branch, while a negative bias favors its reflected counterpart. The exact reflection is

$$
(x,h)\mapsto(-x,-h).
$$

Thus a nonzero $h$ destroys the perfect pitchfork:

- one branch continues smoothly through the former critical neighborhood;
- the other two branches are created at a saddle-node fold;
- changing the sign of $h$ reflects the branch geometry.

![Imperfect bifurcation diagrams for small nonnegative biases](/assets/kapitza-kappa/small_positive_bias_bifurcations.png)

![Reflected imperfect bifurcation diagrams for small nonpositive biases](/assets/kapitza-kappa/small_negative_bias_bifurcations.png)

## 10. Exact saddle-node boundary

At a fold, the equilibrium and tangency conditions hold simultaneously. Eliminating $\kappa$ gives

$$
\boxed{
\kappa_{\mathrm{SN}}
=
\frac{1}{1-x_{\mathrm{SN}}^2},
}
$$

$$
\boxed{
h
=
\operatorname{artanh}(x_{\mathrm{SN}})
-\frac{x_{\mathrm{SN}}}{1-x_{\mathrm{SN}}^2}.
}
$$

Near the symmetric critical point,

$$
h=-\frac{2}{3}x_{\mathrm{SN}}^3+O(x_{\mathrm{SN}}^5),
$$

while

$$
\kappa_{\mathrm{SN}}-1
=x_{\mathrm{SN}}^2+O(x_{\mathrm{SN}}^4).
$$

Combining the expansions yields the cusp scaling law

$$
\boxed{
\kappa_{\mathrm{SN}}-1
\sim
\left(\frac{3|h|}{2}\right)^{2/3}.
}
$$

The exact trajectory and the local asymptotic expression agree closely near $h=0$.

![Exact saddle-node trajectory and the local two-thirds scaling law](/assets/kapitza-kappa/critical_point_trajectory.png)

The numerical fold data are also available as [CSV](/assets/kapitza-kappa/data/critical_trajectory.csv).

## 11. Fixed-source S-curve

Fix $\kappa>1$ and use $h$ as the continuation parameter. The equilibrium family is

$$
h(x)=\operatorname{artanh}(x)-\kappa x,
$$

$$
a(x;\kappa)
=
\frac{\kappa(1+x)}{2(1+\delta)}.
$$

The S-curve should be displayed in the physical $(h,a)$ plane. Its fold coordinates satisfy

$$
x_f=\sqrt{1-\frac{1}{\kappa}},
$$

$$
h_c=\kappa x_f-\operatorname{artanh}(x_f).
$$

For $\kappa=1.5$,

$$
h_c\approx0.2075464553.
$$

The same equilibria can be recovered as intersections of the bulk and interface heat-flux curves:

$$
\widehat q_{\mathrm{bulk}}(s;\kappa)
=
\frac{\delta}{1+\delta}(\kappa-s),
$$

$$
\widehat q_K(s;\kappa,h)
=
\widehat{\beta}_{\kappa,h}(s)s.
$$

![Fixed-source S-curve and equivalent heat-flux intersections](/assets/kapitza-kappa/s_curve.png)

At an interior biased slice, three roots coexist. Two are attracting for the Picard iteration and are separated by one repelling root.

![Biased fixed-point intersections and the two attracting cobweb basins](/assets/kapitza-kappa/biased_fixed_point_cobweb.png)

## 12. Reconstruction of the physical temperature fields

Once a normalized root $x$ is known, the physical coefficients are

$$
a=\frac{\kappa(1+x)}{2(1+\delta)},
\qquad
b=\delta a,
$$

and the jump is

$$
s=\frac{\kappa}{2}(1-x).
$$

The complete steady temperature field is

$$
u(z)=
\begin{cases}
a(z+1), & z\in[-1,0],\\
b(z-1)+\kappa(1-z^2), & z\in[0,1].
\end{cases}
$$

For $\delta=0.5$, $\kappa=1.5$, and $h=0$, the model has three steady states. Their normalized roots and interface traces are approximately

| Branch | $x$ | $a$ | $s$ | $m_{\mathrm{FP}}$ | Picard behavior |
|---|---:|---:|---:|---:|---|
| Lower | $-0.858560$ | $0.070720$ | $1.393920$ | $0.394313$ | Attracting |
| Central | $0$ | $0.500000$ | $0.750000$ | $1.500000$ | Repelling |
| Upper | $0.858560$ | $0.929280$ | $0.106080$ | $0.394313$ | Attracting |

The corresponding piecewise fields satisfy the same homogeneous exterior conditions but have different interface traces, jumps, and heat fluxes.

![Three reconstructed piecewise temperature fields at kappa equals 1.5](/assets/kapitza-kappa/temperature_profiles.png)

The complete benchmark table is available as [CSV](/assets/kapitza-kappa/data/steady_states.csv).

## 13. Numerical strategy

Both implementations use the same workflow:

1. Solve the normalized scalar equation in $(-1,1)$.
2. Recover $a$, $b$, $s$, $\widehat{\beta}$, and $\widehat q$.
3. Evaluate $m_{\mathrm{FP}}=\kappa(1-x^2)$.
4. Generate fixed-point, cobweb, bifurcation, fold, S-curve, and temperature-profile figures.
5. Export transparent CSV tables for independent checking.

In the symmetric case, odd symmetry is used directly. For $\kappa>1$, the positive nonzero root is isolated with a bracketed scalar solver and reflected to obtain the negative root.

For $h\neq0$, the residual

$$
R(x)=\tanh(\kappa x+h)-x
$$

is scanned on the physical interval and every sign-changing bracket is refined. The fold trajectory is computed independently from the analytic saddle-node equations, which provides a useful consistency check.

The Python workflow can be run with

```bash
python3 code/python/kapitza_kappa_model.py
```

and the MATLAB workflow with

```matlab
run('code/matlab/kapitza_kappa_model.m')
```

The default values are

$$
\delta=0.5,
\qquad
\kappa=1.5.
$$

## 14. What this revised formulation clarifies

The revised model separates three levels that should not be conflated.

First, the bulk PDE determines the source-dependent affine-quadratic temperature fields and the exact jump relation

$$
s=\kappa-(1+\delta)a.
$$

Second, the constitutive law determines how the jump and heat flux are coupled. The exponential law is selected so that the normalized equation has an exact sigmoid form.

Third, the scalar Picard map provides a numerical iteration and a local multiplier. Its attracting and repelling branches describe iteration behavior, not automatically the physical stability of an evolutionary heat-transfer problem.

The principal gain over the earlier parameterization is that the bifurcation parameter now has a direct physical interpretation. The threshold

$$
\kappa_c=1
$$

is a critical normalized source strength, not a statement that the two material conductivities are equal. Material contrast remains encoded by

$$
\delta=\frac{\lambda}{\mu}.
$$

## 15. Summary

The source-controlled transmission problem

$$
-\lambda u_I''=0,
\qquad
-\mu u_E''=2\mu\kappa,
$$

together with the exponential Kapitza law, reduces exactly to

$$
x=\tanh(\kappa x+h).
$$

For $h=0$, a supercritical pitchfork occurs at $\kappa_c=1$. A nonzero bias unfolds it into a smooth preferred branch and a saddle-node pair, with

$$
\kappa_{\mathrm{SN}}-1
\sim
\left(\frac{3|h|}{2}\right)^{2/3}.
$$

The normalized analysis can be converted back exactly to the moving physical interval, heat flux, interface jump, and full piecewise temperature fields. This makes the model a compact but physically interpretable benchmark for nonlinear contact resistance, bifurcation analysis, and reproducible scientific computation.
