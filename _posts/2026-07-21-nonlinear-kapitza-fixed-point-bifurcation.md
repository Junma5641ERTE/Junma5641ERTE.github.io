---
layout: post
title: "Nonlinear Kapitza Contact Resistance: Fixed-Point Reduction, Pitchfork Bifurcation, and S-Curve Unfolding"
date: 2026-07-21
categories: [PDE, Bifurcation, Heat Transfer, Research]
---

# Nonlinear Kapitza Contact Resistance: Fixed-Point Reduction, Pitchfork Bifurcation, and S-Curve Unfolding

These notes summarize a small but useful model problem that I developed while studying nonlinear Kapitza contact resistance in composite heat-conduction problems. The original folder contains French and Chinese understanding notes, a detailed English report, an English Beamer presentation, and reproducible Python/MATLAB codes. This post reorganizes the material into a single research note.

The main point is that a one-dimensional nonlinear transmission problem can be reduced exactly to a scalar fixed-point equation. With the exponential conductance used here, the reduced equation becomes

$$
x=\tanh\left(\frac{\alpha}{2}x+h\right).
$$

This equation is simple enough to analyze by hand, but rich enough to display the main mechanisms that appear in nonlinear interface laws: multiplicity, pitchfork bifurcation, imperfection under bias, saddle-node folds, S-curves, and the difference between numerical fixed-point stability and physical stability.

The complete English report and presentation are available here:

- [Detailed report: Nonlinear Kapitza Problem](/assets/kapitza/kapitza_solution_en.pdf)
- [Presentation slides](/assets/kapitza/kapitza_presentation_en.pdf)

## 1. The transmission problem

Consider a one-dimensional steady heat-conduction problem on a two-layer domain. The temperature is piecewise affine, and the two subdomains meet at an interface with a temperature jump. The unknown can be reduced to one scalar parameter \(a\), which determines the affine temperature profile.

The interface jump and heat flux are written as

$$
s=[u]=1-(1+\delta)a,
\qquad
q=\delta a=\beta(s)s,
\qquad \delta>0.
$$

Here:

- \(s=[u]\) is the temperature jump across the interface;
- \(q\) is the steady heat flux;
- \(\delta\) is the material-contrast parameter;
- \(\beta(s)\) is the nonlinear Kapitza conductance.

The physically admissible interval is

$$
0<a<\gamma,\qquad 0<s<1,\qquad \gamma=\frac{1}{1+\delta}.
$$

This interval is not an artificial assumption. It comes from the maximum-principle interpretation of the one-dimensional steady state: the interface traces stay between the prescribed boundary values, and positive conductance gives positive heat flow.

Solving the interface balance for \(a\) gives a scalar fixed-point map

$$
a=F(a)
=
\frac{\beta(1-(1+\delta)a)}
{\delta+(1+\delta)\beta(1-(1+\delta)a)}.
$$

Every fixed point of \(F\) corresponds to one steady transmission state.

## 2. The exponential Kapitza law

The conductance used in the notes follows the teacher's parameter convention:

$$
\boxed{
\beta_h(s)=
\frac{\delta}{1+\delta}
\exp\left[\alpha(1-2s)+2h\right].
}
$$

The symmetric law is obtained when \(h=0\). The parameter \(h\) is a bias parameter; it multiplies the whole conductance by \(e^{2h}\). The variable \(d\) used in the teacher's short Scilab code is the same quantity as \(\delta\) in these notes.

This nonlinear law is always positive. It is decreasing and convex as a function of the jump:

$$
\partial_s\beta_h=-2\alpha\beta_h<0,
\qquad
\partial_{ss}\beta_h=4\alpha^2\beta_h>0.
$$

Substitution into the general fixed-point map gives the sigmoid form

$$
F_{\alpha,h}(a)
=
\frac{\gamma}{2}
\left[
1+\tanh\left(\frac{\alpha}{2}(1-2s)+h\right)
\right].
$$

The important point is that the nonlinear interface law is not merely similar to a sigmoid map. It is exactly equivalent to one after the one-dimensional reduction.

For the teacher's subcritical example with \(\delta=0.5\) and \(\alpha=1.5\), the fixed-point map has a single intersection:

$$
a=\frac13,
\qquad
s=\frac12,
\qquad
\beta=\frac13,
\qquad
q=\frac16.
$$

![Teacher's subcritical fixed-point example reproduced in the project](/assets/kapitza/teacher_fixed_point.png)

## 3. Normalized state and exact fixed-point equation

Introduce the normalized variable

$$
x=\frac{2a}{\gamma}-1=2(1+\delta)a-1,
\qquad -1<x<1.
$$

Then

$$
s=\frac{1-x}{2},
\qquad 1-2s=x.
$$

The steady problem becomes

$$
\boxed{
x=\tanh\left(\frac{\alpha}{2}x+h\right).
}
$$

The material parameter \(\delta\) still matters when reconstructing the physical state \(a\), but it disappears from the normalized root equation. This is why the critical threshold in the teacher's \(\alpha\)-convention is universal:

$$
\alpha_c=2.
$$

## 4. Symmetric case: one root or three roots

Set \(h=0\), and define \(k=\alpha/2\). The residual is

$$
G(x)=\tanh(kx)-x.
$$

The function \(G\) is odd. Moreover,

$$
G'(x)=k\,\operatorname{sech}^2(kx)-1.
$$

If \(k\le 1\), or equivalently \(\alpha\le 2\), the origin is the unique root. If \(k>1\), then the origin remains a root, but two additional nonzero roots appear symmetrically. Thus

$$
\boxed{
\begin{array}{ll}
0<\alpha\le 2: & \text{one steady state},\\
\alpha>2: & \text{three steady states}.
\end{array}
}
$$

The central branch is

$$
x=0,
\qquad
a=\frac{1}{2(1+\delta)},
\qquad
s=\frac12.
$$

The two nonzero branches admit the exact parametrization

$$
\alpha(x)=2\frac{\operatorname{artanh}(x)}{x},
\qquad
a(x)=\frac{1+x}{2(1+\delta)},
\qquad 0<|x|<1.
$$

Since

$$
\frac{\operatorname{artanh}(x)}{x}
=1+\frac{x^2}{3}+O(x^4),
$$

the two outer branches emerge continuously from the central state at \(\alpha=2\). This is the classical geometry of a supercritical pitchfork.

![Fixed-point intersections below, at, and above the universal threshold](/assets/kapitza/fixed_point_regimes.png)

![Analytic pitchfork bifurcation for the symmetric nonlinear Kapitza model](/assets/kapitza/pitchfork_bifurcation.png)

## 5. Fixed-point iteration stability

The scalar iteration associated with the normalized equation is

$$
x_{n+1}
=
\tanh\left(\frac{\alpha}{2}x_n+h\right).
$$

At a fixed point, its multiplier is

$$
\lambda
=
\frac{\alpha}{2}
\operatorname{sech}^2\left(\frac{\alpha}{2}x+h\right)
=
\frac{\alpha}{2}(1-x^2).
$$

A fixed point is attracting for this numerical iteration if \(|\lambda|<1\), and repelling if \(|\lambda|>1\).

In the symmetric case:

- the central branch has multiplier \(\lambda=\alpha/2\);
- it is attracting for \(\alpha<2\);
- it is repelling for \(\alpha>2\);
- at \(\alpha=2\), the linearization is neutral;
- the outer branches are attracting.

At the critical point, the convergence is not geometric. Since

$$
\tanh x=x-\frac{x^3}{3}+O(x^5),
$$

nearby iterates approach the origin only algebraically.

This stability statement is deliberately limited. It classifies the scalar fixed-point algorithm, not automatically the physical stability of a time-dependent heat equation. A physical stability result would require an evolution model and a separate spectral or energy analysis.

![Cobweb paths showing convergence, critical slowing, and outer-branch selection](/assets/kapitza/cobweb_regimes.png)

## 6. Numerical example: \(\delta=0.5,\alpha=3\)

For the symmetric supercritical example with \(\delta=0.5\) and \(\alpha=3\), the physical interval is

$$
\gamma=\frac{1}{1+\delta}=\frac23.
$$

The computations recover three steady states:

| \(x\) | \(a\) | \(s\) | \(\beta(s)\) | \(q\) | \(\lambda\) | Fixed-point iteration |
|---:|---:|---:|---:|---:|---:|:---|
| -0.858560 | 0.047147 | 0.929280 | 0.025367 | 0.023573 | 0.394313 | attracting |
| 0.000000 | 0.333333 | 0.500000 | 0.333333 | 0.166667 | 1.500000 | repelling |
| 0.858560 | 0.619520 | 0.070720 | 4.380078 | 0.309760 | 0.394313 | attracting |

The identity

$$
q=\delta a=\beta(s)s
$$

holds at every steady state. The lower and upper outer states correspond to the two attracting branches; the central state is the symmetry-preserving but repelling state for the fixed-point iteration.

## 7. Bias and the imperfect pitchfork

The bias parameter enters through

$$
\beta_h(s)=e^{2h}\beta_0(s).
$$

Thus \(h\) rescales the entire interface conductance without changing its shape as a function of \(s\).

In the normalized equation, solve for \(h\):

$$
h(x)=\operatorname{artanh}(x)-kx,
\qquad k=\frac{\alpha}{2}.
$$

When \(k>1\), the fold condition is

$$
\frac{dh}{dx}
=
\frac{1}{1-x^2}-k=0.
$$

Therefore

$$
x_f=\sqrt{1-\frac1k},
\qquad
h_c=kx_f-\operatorname{artanh}(x_f)>0.
$$

The two folds are

$$
(h,x)=(h_c,-x_f),
\qquad
(h,x)=(-h_c,x_f).
$$

Consequently:

- \(|h|<h_c\): three distinct roots;
- \(|h|=h_c\): two distinct roots, one of them double;
- \(|h|>h_c\): one root.

For \(\alpha=3\), \(k=1.5\), and

$$
x_f=0.577350269190,
\qquad
h_c=0.207546455322.
$$

The sample biased case uses \(h=h_c/2=0.103773227661\), giving:

| \(x\) | \(a\) | \(s\) | \(\beta_h(s)\) | \(q\) | \(\lambda\) | Fixed-point iteration |
|---:|---:|---:|---:|---:|---:|:---|
| -0.798116 | 0.067295 | 0.899058 | 0.037425 | 0.033647 | 0.544516 | attracting |
| -0.214294 | 0.261902 | 0.607147 | 0.215683 | 0.130951 | 1.431117 | repelling |
| 0.894934 | 0.631645 | 0.052533 | 6.011867 | 0.315822 | 0.298641 | attracting |

These same states can be represented either as fixed-point intersections of \(F_{\alpha,h}\) with the identity or as intersections of the two heat-flux curves

$$
q_K(s)=\beta_h(s)s,
\qquad
q_{\mathrm{bulk}}(s)=\frac{\delta}{1+\delta}(1-s).
$$

![A fixed-h slice of the biased model and two cobweb basins](/assets/kapitza/s_curve_fixed_point_cobweb.png)

## 8. Continuing in \(\alpha\) with small fixed bias

To keep \(\alpha\) as the bifurcation parameter, fix a small nonzero \(h\) and solve the steady equation for \(\alpha\). For \(x\ne0\),

$$
\alpha(x;h)=
2\frac{\operatorname{artanh}(x)-h}{x}.
$$

A nonzero bias breaks the reflection symmetry \(x\mapsto -x\). The perfect pitchfork no longer exists. Instead:

- the branch favored by the sign of \(h\) continues smoothly;
- the other two branches are born in a saddle-node bifurcation;
- the saddle-node threshold moves to larger \(\alpha\) as \(|h|\) increases.

Combining the steady equation with the tangency condition gives

$$
\alpha_{\mathrm{SN}}
=
\frac{2}{1-x_{\mathrm{SN}}^2},
\qquad
h
=
\operatorname{artanh}(x_{\mathrm{SN}})
-
\frac{x_{\mathrm{SN}}}{1-x_{\mathrm{SN}}^2}.
$$

For

$$
h\in\{0,\ 0.002,\ 0.005,\ 0.010,\ 0.020,\ 0.040\},
$$

the corresponding multiplicity thresholds are

$$
\alpha_{\mathrm{SN}}
\in
\{2.000000,\ 2.041774,\ 2.077216,\ 2.123116,\ 2.196796,\ 2.315799\}.
$$

Thus even a very small bias destroys the exact pitchfork. The shift is not linear in \(h\). Near the symmetric critical state,

$$
\operatorname{artanh}(x)
=
x+\frac{x^3}{3}+\frac{x^5}{5}+O(x^7),
$$

while

$$
\frac{x}{1-x^2}
=
x+x^3+x^5+O(x^7).
$$

The linear terms cancel in the fold equation, leaving

$$
h=-\frac{2}{3}x_{\mathrm{SN}}^3+O(x_{\mathrm{SN}}^5).
$$

Hence

$$
x_{\mathrm{SN}}=O(|h|^{1/3}),
\qquad
\alpha_{\mathrm{SN}}-2=O(|h|^{2/3}).
$$

More precisely, for \(h\to0^+\),

$$
\alpha_{\mathrm{SN}}-2
\sim
2\left(\frac{3h}{2}\right)^{2/3}.
$$

This is the local cusp scaling law associated with the imperfect pitchfork.

![Continuation in alpha for small positive fixed biases](/assets/kapitza/small_bias_alpha_bifurcation.png)

## 9. Negative bias and reflection symmetry

The normalized equation is invariant under the simultaneous reflection

$$
(x,h)\mapsto(-x,-h).
$$

Because

$$
a=\frac{\gamma}{2}(1+x),
$$

changing the sign of the bias reflects the physical branches through the midpoint \(a=\gamma/2\):

$$
a(\alpha;-h)=\gamma-a(\alpha;h),
\qquad
\alpha_{\mathrm{SN}}(-h)=\alpha_{\mathrm{SN}}(h).
$$

Therefore negative bias is not a different mechanism. It favors the lower attracting branch and moves the saddle-node pair to the upper side, while the threshold for the same \(|h|\) is unchanged.

![Continuation in alpha for small negative fixed biases](/assets/kapitza/small_negative_bias_alpha_bifurcation.png)

The saddle-node trajectory makes the imperfect-bifurcation geometry clearer. In the \((h,\alpha)\)-plane, the exact multiplicity boundary is compared with the local \(|h|^{2/3}\) approximation; in the state plot, positive and negative biases move the saddle node to opposite sides of \(a=\gamma/2\).

![Trajectory of the saddle-node point replacing the pitchfork critical point](/assets/kapitza/critical_point_trajectory.png)

## 10. The S-curve with \(h\) as continuation parameter

The S-curve viewpoint fixes \(\alpha\) and uses \(h\) as the continuation parameter. The curve is

$$
h(x)=\operatorname{artanh}(x)-kx.
$$

For \(\alpha=3\), the two folds occur at \(h=\pm h_c\), with

$$
h_c=0.207546455322.
$$

The plotted interval is restricted to approximately

$$
[-2.2h_c,2.2h_c],
$$

so that the folds are not visually compressed.

The S-curve is the actual bifurcation diagram for \(h\)-continuation. The fixed-\(h\) cobweb figure is a slice through this family; its axes are \(a_n\) and \(a_{n+1}\) because once \(h\) is fixed, the iterated variable is the state \(a\), not the parameter \(h\).

![S-curve continuation and equivalent heat-flux intersections](/assets/kapitza/s_curve.png)

## 11. Fixed-state inversion for \(h\)

A complementary inverse calculation prescribes a target state \(a=a_*\) and recovers the compatible bias \(h=h_*\). With

$$
x_*=2(1+\delta)a_*-1,
$$

define

$$
R_{a_*}(h)
=
\tanh(kx_*+h)-x_*.
$$

The compatible bias is

$$
h_*=\operatorname{artanh}(x_*)-kx_*.
$$

The code also visualizes a relaxed fixed-point recovery procedure

$$
h_{n+1}
=
h_n+x_*-\tanh(kx_*+h_n).
$$

This is parameter inversion. It is not a physical time-dependent evolution, and its convergence should not be used as a physical stability criterion.

![Fixed-state recovery of the compatible bias parameter](/assets/kapitza/fixed_a_h_iteration.png)

## 12. Reproducible computation

Both the Python and MATLAB implementations follow the same four stages:

1. Evaluate \(\beta_h(s)\) and \(F_{\alpha,h}(a)\).
2. Solve the normalized residual

   $$
   R(x)=\tanh(\alpha x/2+h)-x.
   $$

3. Reconstruct \(a\), \(s\), \(\beta\), \(q\), and \(\lambda\).
4. Generate fixed-point intersections, cobwebs, bifurcation curves, S-curves, and CSV data.

For \(\alpha>2\), the residual has stationary points at

$$
x=\pm\sqrt{1-\frac{2}{\alpha}}.
$$

The code inserts these points into the search grid, so that each subinterval is monotone. Then it keeps roots at endpoints, applies a bracketed scalar solver on sign-changing intervals, sorts the candidates, and removes duplicates. This is important near folds, where a double root may otherwise be missed.

The Python command is

```bash
python3 code/python/kapitza_numerical.py \
  --delta 0.5 --alpha 3 --output-dir output/python
```

The MATLAB command is

```matlab
run('code/matlab/kapitza_numerical.m')
```

The two implementations were written to match the same mathematics and regenerate the same canonical figures. Python also writes CSV tables for the reconstructed states and iteration multipliers.

## 13. What this model clarifies

This model is small, but it clarifies several points that matter for nonlinear interface problems.

First, nonlinear Kapitza laws can create non-uniqueness even in a one-dimensional steady problem. The mechanism is visible after reduction to one scalar fixed-point equation.

Second, the threshold \(\alpha_c=2\) is not a numerical accident. In the teacher's parameter convention it follows from the slope of the normalized sigmoid at the origin.

Third, the symmetric pitchfork is structurally fragile. A nonzero bias \(h\) unfolds it into saddle-node behavior and an S-curve. The local displacement law \(\alpha_{\mathrm{SN}}-2\sim C|h|^{2/3}\) is a useful warning: the response to small bias is singular, not linear.

Fourth, fixed-point iteration stability should not be confused with physical stability. The multiplier \(\lambda\) tells us whether a numerical iteration converges to a fixed point. A time-dependent heat equation would require its own stability analysis.

Finally, the same mathematics can be made reproducible in both Python and MATLAB. The key is to keep the parameter convention explicit, reconstruct physical variables after solving the normalized equation, and distinguish carefully between continuation parameters and iterated state variables.

## 14. Short summary

The nonlinear conductance

$$
\beta_h(s)=
\frac{\delta}{1+\delta}
\exp[\alpha(1-2s)+2h]
$$

turns the one-dimensional Kapitza transmission problem into

$$
x=\tanh\left(\frac{\alpha}{2}x+h\right).
$$

For \(h=0\), the universal threshold is \(\alpha_c=2\). Below threshold there is one steady state; above threshold there are three, organized by a supercritical pitchfork. For \(h\ne0\), the pitchfork unfolds into saddle-node folds and an S-curve. The numerical codes reproduce the analytical branches, cobweb behavior, heat-flux intersections, fold trajectory, and fixed-state inversion.

These notes are a compact model for thinking about nonlinear contact resistance, bifurcation, and numerical branch reconstruction in more general composite-domain heat-transfer problems.
