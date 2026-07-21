#!/usr/bin/env python3
"""Numerically study the symmetric and imperfect nonlinear Kapitza models.

Every computational step is documented in English so that this file can also
be used as a teaching or discussion script.

The teacher's short Scilab example writes the material parameter as ``d``;
this project uses the paper notation ``delta`` for the same quantity.
"""

from __future__ import annotations  # Allow modern type hints without eager evaluation.

import argparse  # Parse model parameters supplied on the command line.
from pathlib import Path  # Handle output folders and filenames portably.

import matplotlib  # Configure the plotting backend before importing pyplot.
import numpy as np  # Provide vectorized arrays and elementary functions.
from scipy.optimize import brentq  # Solve bracketed scalar equations robustly.

matplotlib.use("Agg")  # Select a non-interactive backend suitable for batch runs.
import matplotlib.pyplot as plt  # Create and save all numerical figures.

plt.rcParams.update(  # Apply one publication-style typography policy to every Python figure.
    {
        "font.family": "serif",  # Use a serif face for titles, axes, ticks, and legends.
        "font.serif": ["Times New Roman"],  # Select the requested Times New Roman family.
        "mathtext.fontset": "stix",  # Use Times-compatible LaTeX-style mathematical glyphs.
        "axes.unicode_minus": False,  # Keep minus signs compatible with the selected serif font.
    }
)  # Finish the global plotting-style configuration.


def gamma(delta: float) -> float:
    """Return the upper physical bound gamma=1/(1+delta) for the state a."""
    return 1.0 / (1.0 + delta)  # Evaluate the definition of gamma.


def alpha_critical(delta: float) -> float:
    """Return the pitchfork threshold in the teacher's alpha convention."""
    _ = delta  # Keep delta in the public interface although the teacher's threshold is universal.
    return 2.0  # The central fixed-point slope alpha/2 equals one here.


def f_map(
    a: np.ndarray | float, alpha: float, delta: float, h: float = 0.0
) -> np.ndarray | float:
    """Evaluate the prescribed sigmoid fixed-point map F_{alpha,h}(a)."""
    gam = gamma(delta)  # Compute the admissible interval length for a.
    normalized_state = 2.0 * a / gam - 1.0  # Convert a from (0,gamma) to x in (-1,1).
    argument = 0.5 * alpha * normalized_state + h  # Use the teacher's normalized alpha.
    return 0.5 * gam * (1.0 + np.tanh(argument))  # Map the tanh range to (0,gamma).


def f_prime(
    a: np.ndarray | float, alpha: float, delta: float, h: float = 0.0
) -> np.ndarray | float:
    """Evaluate the derivative of the sigmoid map with respect to a."""
    gam = gamma(delta)  # Reuse the physical interval length.
    normalized_state = 2.0 * a / gam - 1.0  # Convert the physical state to x.
    z = 0.5 * alpha * normalized_state + h  # Evaluate the teacher-convention argument.
    return 0.5 * alpha / np.cosh(z) ** 2  # Differentiate the rescaled sigmoid map.


def beta(
    s: np.ndarray | float, alpha: float, delta: float, h: float = 0.0
) -> np.ndarray | float:
    """Evaluate the reconstructed nonlinear Kapitza conductance beta_h(s)."""
    return delta / (1.0 + delta) * np.exp(  # Multiply the exponential by its prefactor.
        alpha * (1.0 - 2.0 * s) + 2.0 * h
    )  # The factor exp(2h) changes the overall interface conductance amplitude.


def imperfect_roots(alpha: float, delta: float, h: float = 0.0) -> list[float]:
    """Return every root of x=tanh(kappa*x+h) in (-1,1)."""
    kappa = alpha / alpha_critical(delta)  # Convert alpha to the normalized slope.
    eps = 1.0e-12  # Stay slightly inside the open physical interval (-1,1).
    residual = lambda x: np.tanh(kappa * x + h) - x  # Define the root residual.

    split_points = [-1.0 + eps, 1.0 - eps]  # Start with the two interval endpoints.
    if kappa > 1.0:  # Add stationary points only in the multiple-solution regime.
        turning_x = np.sqrt(1.0 - 1.0 / kappa)  # Solve residual'(x)=0 analytically.
        split_points[1:1] = [-turning_x, turning_x]  # Split into monotone subintervals.

    roots: list[float] = []  # Collect ordinary roots and double roots at folds.
    tolerance = 2.0e-11  # Decide when a subinterval endpoint is already a root.
    for left, right in zip(split_points[:-1], split_points[1:]):  # Inspect each interval.
        f_left = residual(left)  # Evaluate the residual at the left endpoint.
        f_right = residual(right)  # Evaluate the residual at the right endpoint.
        if abs(f_left) < tolerance:  # Retain a root located exactly at a fold.
            roots.append(left)  # Store the left endpoint root.
        if f_left * f_right < 0.0:  # A sign change guarantees an interior root.
            roots.append(brentq(residual, left, right))  # Refine it with Brent's method.
        if abs(f_right) < tolerance:  # Retain a root at the right endpoint as well.
            roots.append(right)  # Store the right endpoint root.

    roots.sort()  # Put roots in increasing physical order.
    unique: list[float] = []  # Prepare a duplicate-free result list.
    for root in roots:  # Compare each candidate with the preceding accepted root.
        if not unique or abs(root - unique[-1]) > 1.0e-8:  # Reject numerical duplicates.
            unique.append(root)  # Keep a genuinely distinct root.
    return unique  # Return all distinct normalized steady states.


def normalized_roots(alpha: float, delta: float) -> list[float]:
    """Return the roots of the symmetric problem by setting h=0."""
    return imperfect_roots(alpha, delta, h=0.0)  # Reuse the general root solver.


def fold_data(alpha: float, delta: float) -> dict[str, float] | None:
    """Analytic folds for h(x)=atanh(x)-kappa*x."""
    kappa = alpha / alpha_critical(delta)  # Compute the normalized sigmoid slope.
    if kappa <= 1.0:  # The monotone subcritical relation has no folds.
        return None  # Signal that an S-curve cannot occur.
    x_fold = float(np.sqrt(1.0 - 1.0 / kappa))  # Evaluate the positive fold state.
    h_limit = float(kappa * x_fold - np.arctanh(x_fold))  # Evaluate |h| at a fold.
    return {"kappa": kappa, "x_fold": x_fold, "h_limit": h_limit}  # Package results.


def alpha_fold_for_bias(h: float) -> tuple[float, float]:
    """Return the saddle-node state and alpha value for one nonzero fixed bias."""
    if h == 0.0:  # Recover the symmetric pitchfork limit without a root search.
        return 0.0, 2.0  # The perfect pitchfork occurs at x=0 and alpha=2.
    sign = -np.sign(h)  # The saddle node lies on the side opposed to the bias.
    residual = lambda x: np.arctanh(x) - x / (1.0 - x * x) - h  # Combine both fold equations.
    if sign < 0.0:  # Bracket the fold on the negative-x side for a positive bias.
        x_fold = brentq(residual, -0.999999, -1.0e-12)  # Isolate the unique fold state.
    else:  # Bracket the reflected fold for a negative bias.
        x_fold = brentq(residual, 1.0e-12, 0.999999)  # Isolate the unique fold state.
    alpha_fold = 2.0 / (1.0 - x_fold * x_fold)  # Apply the tangency condition lambda=1.
    return float(x_fold), float(alpha_fold)  # Return the state and control-parameter threshold.


def physical_states(
    alpha: float, delta: float, h: float = 0.0
) -> list[dict[str, float]]:
    """Convert every normalized root into physical Kapitza state variables."""
    gam = gamma(delta)  # Compute the physical scale for a.
    states = []  # Store one dictionary for each steady state.
    kappa = alpha / alpha_critical(delta)  # Compute the normalized iteration slope.
    for x in imperfect_roots(alpha, delta, h):  # Reconstruct every detected branch.
        a = 0.5 * gam * (1.0 + x)  # Convert x back to the affine-solution parameter a.
        s = 0.5 * (1.0 - x)  # Convert x to the interface jump s=[u].
        b = float(beta(s, alpha, delta, h))  # Evaluate the conductance at this jump.
        states.append(
            {
                "x": x,  # Store the normalized state.
                "a": a,  # Store the affine-solution parameter.
                "s": s,  # Store the temperature jump.
                "beta": b,  # Store the nonlinear conductance.
                "q": delta * a,  # Store the heat flux q=delta*a=beta(s)*s.
                "iteration_slope": float(kappa * (1.0 - x * x)),  # Store dF/dx.
            }
        )  # Finish this state's dictionary and append it.
    return states  # Return all reconstructed physical states.


def write_states(path: Path, alpha: float, delta: float) -> None:
    """Write the symmetric steady states and iteration slopes to a CSV file."""
    lines = [  # Start the output with a descriptive comma-separated header.
        "branch,x,a,s,beta,q,iteration_slope,iteration_attractive",
    ]
    for index, state in enumerate(physical_states(alpha, delta), start=1):  # Number roots.
        attractive = abs(state["iteration_slope"]) < 1.0  # Apply the fixed-point criterion.
        lines.append(
            f"{index},{state['x']:.12g},{state['a']:.12g},"
            f"{state['s']:.12g},{state['beta']:.12g},{state['q']:.12g},"
            f"{state['iteration_slope']:.12g},{attractive}"
        )  # Serialize one physical state with twelve significant digits.
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")  # Save the CSV text.


def plot_summary(path: Path, delta: float, alpha: float) -> None:
    """Create the four-panel summary used in the report and presentation."""
    gam = gamma(delta)  # Compute the physical a interval.
    alpha_c = alpha_critical(delta)  # Compute the pitchfork threshold.
    alpha_values = [1.5, alpha_c, alpha]  # Include the teacher's example, threshold, and chosen case.
    colors = ["#3478a8", "#777777", "#c23b23"]  # Assign consistent curve colors.

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.2), constrained_layout=True)  # Build the panel grid.

    a_grid = np.linspace(0.0, gam, 1000)  # Sample the complete physical interval for a.
    ax = axes[0, 0]  # Select the fixed-point-intersection panel.
    ax.plot(a_grid, a_grid, "k--", linewidth=1.2, label=r"$y=a$")  # Draw the identity line.
    for alpha, color in zip(alpha_values, colors):  # Draw one sigmoid for each regime.
        ax.plot(
            a_grid,
            f_map(a_grid, alpha, delta),
            color=color,
            label=rf"$\alpha={alpha:.2f}$",
        )  # Plot the selected fixed-point map.
    ax.set(xlabel=r"$a$", ylabel=r"$f_\alpha(a)$", title="Fixed-point intersections")  # Label the panel.
    ax.legend(fontsize=8)  # Identify all fixed-point curves.
    ax.grid(alpha=0.25)  # Add a light reading grid.

    ax = axes[0, 1]  # Select the pitchfork-bifurcation panel.
    alpha_grid = np.linspace(0.25 * alpha_c, max(3.0 * alpha_c, 1.2 * alpha), 500)  # Set the alpha range.
    stable_central = alpha_grid <= alpha_c  # Identify the attracting part of the central branch.
    unstable_central = alpha_grid >= alpha_c  # Identify the repelling part after bifurcation.
    ax.plot(alpha_grid[stable_central], np.full(np.count_nonzero(stable_central), 0.5 * gam),
            color="#c23b23", linestyle="-", label="stable central")  # Draw it solid.
    ax.plot(alpha_grid[unstable_central], np.full(np.count_nonzero(unstable_central), 0.5 * gam),
            color="#c23b23", linestyle="--", label="unstable central")  # Draw it dashed.
    x_grid = np.linspace(1.0e-4, 0.985, 1000)  # Parameterize the positive outer branch by x.
    alpha_outer = alpha_c * np.arctanh(x_grid) / x_grid  # Evaluate the analytic branch formula.
    mask = alpha_outer <= alpha_grid[-1]  # Keep only points inside the displayed alpha range.
    ax.plot(alpha_outer[mask], 0.5 * gam * (1.0 + x_grid[mask]), color="#3478a8",
            label="stable outer")  # Draw the attracting upper branch as solid.
    ax.plot(alpha_outer[mask], 0.5 * gam * (1.0 - x_grid[mask]), color="#3478a8")  # Draw the lower branch.
    ax.axvline(alpha_c, color="k", linestyle="--", linewidth=1.0)  # Mark the critical parameter.
    ax.set(
        xlabel=r"$\alpha$",
        ylabel=r"$a$",
        title=rf"Bifurcation diagram ($\alpha_c={alpha_c:.2f}$)",
    )  # Label the bifurcation diagram.
    ax.legend(fontsize=7)  # Explain the solid-stable and dashed-unstable convention.
    ax.grid(alpha=0.25)  # Add a light reading grid.

    s_grid = np.linspace(0.0, 1.0, 1000)  # Sample all admissible temperature jumps.
    ax = axes[1, 0]  # Select the nonlinear-law panel.
    ax.plot(s_grid, beta(s_grid, alpha, delta), color="#3478a8")  # Plot beta(s).
    ax.set(
        xlabel=r"$s=[u]$",
        ylabel=r"$\beta(s)$",
        title=rf"Constructed Kapitza law ($\alpha={alpha:.2f}$)",
    )  # Label the conductance-law panel.
    ax.grid(alpha=0.25)  # Add a light reading grid.

    ax = axes[1, 1]  # Select the heat-flux-intersection panel.
    q_k = beta(s_grid, alpha, delta) * s_grid  # Compute the nonlinear interface flux.
    q_bulk = delta / (1.0 + delta) * (1.0 - s_grid)  # Compute the bulk-required flux.
    ax.plot(s_grid, q_k, color="#c23b23", label=r"$q_K=\beta(s)s$")  # Draw the interface flux.
    ax.plot(s_grid, q_bulk, color="#3478a8", label=r"$q_{\rm bulk}$")  # Draw the bulk flux.
    for state in physical_states(alpha, delta):  # Mark every steady-state intersection.
        ax.plot(state["s"], state["q"], "ko", markersize=4)  # Add a black root marker.
    ax.set(xlabel=r"$s=[u]$", ylabel=r"$q$", title="Flux intersections")  # Label the flux panel.
    ax.legend(fontsize=8)  # Identify the two flux curves.
    ax.grid(alpha=0.25)  # Add a light reading grid.

    fig.suptitle(rf"Nonlinear Kapitza problem, $\delta={delta:g}$", fontsize=14)  # Add the figure title.
    fig.savefig(path, dpi=220)  # Save a high-resolution PNG file.
    plt.close(fig)  # Release the Matplotlib figure from memory.


def cobweb_segments(
    a0: float, alpha: float, delta: float, iterations: int = 25, h: float = 0.0
) -> tuple[list[float], list[float]]:
    """Build the alternating vertical and horizontal segments of a cobweb plot."""
    xs = [a0]  # Start the polyline at the selected initial state.
    ys = [0.0]  # Start from the horizontal axis.
    a = a0  # Initialize the current fixed-point iterate.
    for _ in range(iterations):  # Generate the requested number of iterations.
        next_a = float(f_map(a, alpha, delta, h))  # Apply the possibly biased map once.
        xs.extend([a, next_a])  # Add vertical and horizontal x coordinates.
        ys.extend([next_a, next_a])  # Add the matching y coordinates.
        a = next_a  # Advance to the next iterate.
    return xs, ys  # Return coordinates ready for one plotting call.


def plot_cobweb(path: Path, delta: float, alpha: float) -> None:
    """Plot two fixed-point iterations started on opposite sides of the center."""
    gam = gamma(delta)  # Compute the admissible a interval.
    a_grid = np.linspace(0.0, gam, 1000)  # Sample the map smoothly.

    fig, ax = plt.subplots(figsize=(6.4, 6.0), constrained_layout=True)  # Create a square plotting area.
    ax.plot(a_grid, a_grid, "k--", linewidth=1.1, label=r"$y=a$")  # Draw the identity line.
    ax.plot(a_grid, f_map(a_grid, alpha, delta), color="#3478a8", label=r"$y=f_\alpha(a)$")  # Draw the map.
    for a0, color in [(0.18 * gam, "#2a8f55"), (0.82 * gam, "#c23b23")]:  # Choose two basins.
        xs, ys = cobweb_segments(a0, alpha, delta)  # Generate this iteration path.
        ax.plot(xs, ys, color=color, linewidth=0.9, alpha=0.9, label=rf"$a_0={a0:.3f}$")  # Plot it.
    for state in physical_states(alpha, delta):  # Mark all fixed points.
        ax.plot(state["a"], state["a"], "ko", markersize=4)  # Place each marker on y=a.
    ax.set(
        xlim=(0.0, gam),
        ylim=(0.0, gam),
        xlabel=r"$a_n$",
        ylabel=r"$a_{n+1}$",
        title=rf"Cobweb plot: $\delta={delta:g}$, $\alpha={alpha:.2f}$",
    )  # Set physical limits, labels, and the title.
    ax.legend(fontsize=8)  # Identify both iteration paths and reference curves.
    ax.grid(alpha=0.2)  # Add a light reading grid.
    fig.savefig(path, dpi=220)  # Save the cobweb plot.
    plt.close(fig)  # Release the figure from memory.


def plot_fixed_point_regimes(path: Path, delta: float) -> None:
    """Compare fixed-point intersections below, at, and above alpha_c=2."""
    gam = gamma(delta)  # Compute the complete physical interval for a.
    alpha_values = [1.5, 2.0, 3.0]  # Select subcritical, critical, and supercritical cases.
    a_grid = np.linspace(0.0, gam, 2000)  # Sample every map on a common dense grid.
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.4), constrained_layout=True)  # Build three panels.

    for ax, alpha_value in zip(axes, alpha_values):  # Visit each bifurcation regime once.
        states = physical_states(alpha_value, delta)  # Compute every physical fixed point.
        ax.plot(a_grid, f_map(a_grid, alpha_value, delta), color="#1857c9",
                linewidth=2.0, label=r"$F(a)$")  # Draw the nonlinear map.
        ax.plot(a_grid, a_grid, color="#d62728", linestyle="--",
                linewidth=1.4, label=r"$a$")  # Draw the identity line.
        ax.plot([state["a"] for state in states], [state["a"] for state in states],
                "ko", markersize=5, label="Fixed points")  # Mark every intersection.
        ax.set(xlim=(0.0, gam), ylim=(0.0, gam), xlabel=r"$a$", ylabel=r"$F(a)$",
               title=rf"$\alpha={alpha_value:.1f}$: {len(states)} fixed point(s)")  # Label the panel.
        ax.grid(alpha=0.25)  # Add a light grid for visual comparison.

    axes[0].legend(fontsize=8, loc="upper left")  # Explain the common line and marker styles once.
    fig.suptitle("Fixed-point geometry across the pitchfork threshold", fontsize=14)  # Add the message.
    fig.savefig(path, dpi=220)  # Save the fixed-point comparison figure.
    plt.close(fig)  # Release the Matplotlib resources.


def plot_cobweb_regimes(path: Path, delta: float) -> None:
    """Compare fixed-point iteration dynamics across the three alpha regimes."""
    gam = gamma(delta)  # Compute the physical interval length for a.
    alpha_values = [1.5, 2.0, 3.0]  # Select one case in each bifurcation regime.
    default_initial_values = [0.12 * gam, 0.88 * gam]  # Keep wide starts for the one-root regimes.
    path_colors = ["#1b8f4d", "#c43b2b"]  # Assign stable colors to the two initial states.
    a_grid = np.linspace(0.0, gam, 2000)  # Sample all maps on the same grid.
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.4), constrained_layout=True)  # Build three panels.

    for ax, alpha_value in zip(axes, alpha_values):  # Visit each selected alpha value.
        states = physical_states(alpha_value, delta)  # Compute all fixed points in this panel.
        if alpha_value > alpha_critical(delta):  # Detect the three-root supercritical panel.
            middle_fixed_point = states[1]["a"]  # Read the repelling central basin boundary.
            initial_offset = 0.08 * gam  # Stay close enough to reveal departure from the center.
            panel_initial_values = [  # Place one initial condition on each side of the separator.
                middle_fixed_point - initial_offset,  # Select the lower attracting basin.
                middle_fixed_point + initial_offset,  # Select the upper attracting basin.
            ]
        else:  # Preserve the original comparison for the one-root and critical panels.
            panel_initial_values = default_initial_values  # Use starts spanning most of the interval.
        ax.plot(a_grid, f_map(a_grid, alpha_value, delta), color="#1857c9",
                linewidth=1.8, label=r"$F(a)$")  # Draw the nonlinear map.
        ax.plot(a_grid, a_grid, "k--", linewidth=1.2, label=r"$a$")  # Draw the identity line.
        for a0, color, label in zip(panel_initial_values, path_colors,
                                    ["lower start", "upper start"]):  # Draw both basins.
            xs, ys = cobweb_segments(a0, alpha_value, delta, iterations=35)  # Build one path.
            ax.plot(xs, ys, color=color, linewidth=0.8, label=label)  # Draw the cobweb segments.
        ax.plot([state["a"] for state in states], [state["a"] for state in states],
                "ko", markersize=4, label="Fixed points")  # Mark all fixed points.
        ax.set(xlim=(0.0, gam), ylim=(0.0, gam), xlabel=r"$a_n$", ylabel=r"$a_{n+1}$",
               title=rf"$\alpha={alpha_value:.1f}$")  # Apply common bounds and labels.
        ax.grid(alpha=0.25)  # Add a light reading grid.

    axes[0].legend(fontsize=7, loc="upper left")  # Explain the common styles once.
    fig.suptitle("Cobweb iterations reveal stability and branch selection", fontsize=14)  # Add the message.
    fig.savefig(path, dpi=220)  # Save the cobweb-regime comparison.
    plt.close(fig)  # Release the Matplotlib resources.


def plot_small_bias_alpha_bifurcation(path: Path, delta: float) -> None:
    """Compare six small fixed biases while alpha remains the control parameter."""
    gam = gamma(delta)  # Compute the physical scale that maps x to a.
    alpha_min = 0.25  # Start below the symmetric critical value.
    alpha_max = 5.0  # Continue far enough to display all asymptotic branches.
    h_values = [0.0, 0.002, 0.005, 0.010, 0.020, 0.040]  # Resolve the progressive unfolding.
    stable_color = "#1857c9"  # Use one color for attracting fixed-point branches.
    unstable_color = "#c23b23"  # Use a contrasting color for repelling branches.
    fig, axes = plt.subplots(2, 3, figsize=(12.0, 7.6), sharex=True, sharey=True,
                             constrained_layout=True)  # Build six matched continuation panels.

    for index, (ax, h_small) in enumerate(zip(axes.flat, h_values)):  # Visit each selected bias.
        if h_small == 0.0:  # Draw the exact reflection-symmetric pitchfork as the baseline.
            alpha_c = alpha_critical(delta)  # Read the universal symmetric threshold.
            central_stable = np.linspace(alpha_min, alpha_c, 400)  # Sample stable center.
            central_unstable = np.linspace(alpha_c, alpha_max, 600)  # Sample unstable center.
            x_outer = np.linspace(1.0e-4, 0.995, 2200)  # Parameterize both outer branches.
            alpha_outer = 2.0 * np.arctanh(x_outer) / x_outer  # Recover alpha on the branches.
            outer_mask = alpha_outer <= alpha_max  # Retain the displayed part of each branch.
            ax.plot(central_stable, np.full_like(central_stable, 0.5 * gam),
                    color=stable_color, linewidth=1.8, label="attracting")  # Stable center.
            ax.plot(central_unstable, np.full_like(central_unstable, 0.5 * gam),
                    color=unstable_color, linestyle="--", linewidth=1.8,
                    label="repelling")  # Unstable center.
            ax.plot(alpha_outer[outer_mask], 0.5 * gam * (1.0 + x_outer[outer_mask]),
                    color=stable_color, linewidth=1.8)  # Stable upper branch.
            ax.plot(alpha_outer[outer_mask], 0.5 * gam * (1.0 - x_outer[outer_mask]),
                    color=stable_color, linewidth=1.8)  # Stable lower branch.
            ax.plot(alpha_c, 0.5 * gam, "ko", markersize=4.5,
                    label="critical point")  # Mark the perfect pitchfork.
            ax.set_title(r"$h=0$: perfect pitchfork")  # Identify the symmetric reference.
        else:  # Continue the imperfect branches at one positive fixed bias.
            x_negative = np.linspace(-0.995, -1.0e-4, 5000)  # Sample both negative branches.
            x_positive = np.linspace(np.tanh(h_small) + 1.0e-6, 0.995, 3000)  # Favored branch.
            alpha_negative = 2.0 * (np.arctanh(x_negative) - h_small) / x_negative  # Continue roots.
            alpha_positive = 2.0 * (np.arctanh(x_positive) - h_small) / x_positive  # Continue roots.
            lambda_negative = 0.5 * alpha_negative * (1.0 - x_negative * x_negative)  # Slopes.
            negative_visible = (alpha_negative >= alpha_min) & (alpha_negative <= alpha_max)  # Clip.
            negative_stable = negative_visible & (lambda_negative < 1.0)  # Stable lower branch.
            negative_unstable = negative_visible & (lambda_negative >= 1.0)  # Repelling branch.
            positive_visible = (alpha_positive >= alpha_min) & (alpha_positive <= alpha_max)  # Clip.
            ax.plot(alpha_positive[positive_visible],
                    0.5 * gam * (1.0 + x_positive[positive_visible]),
                    color=stable_color, linewidth=1.8, label="attracting")  # Favored branch.
            ax.plot(alpha_negative[negative_stable],
                    0.5 * gam * (1.0 + x_negative[negative_stable]),
                    color=stable_color, linewidth=1.8)  # Stable branch born at the fold.
            ax.plot(alpha_negative[negative_unstable],
                    0.5 * gam * (1.0 + x_negative[negative_unstable]),
                    color=unstable_color, linestyle="--", linewidth=1.8,
                    label="repelling")  # Unstable partner branch.
            x_fold, alpha_fold = alpha_fold_for_bias(h_small)  # Locate this saddle node.
            a_fold = 0.5 * gam * (1.0 + x_fold)  # Convert the fold state to physical a.
            ax.plot(alpha_fold, a_fold, "ko", markersize=4.5,
                    label="saddle node")  # Mark the shifted multiplicity threshold.
            ax.set_title(rf"$h={h_small:.3f}$, $\alpha_{{\rm SN}}={alpha_fold:.3f}$")  # Show shift.
        ax.set_xlim(alpha_min, alpha_max)  # Use one alpha range for all six panels.
        ax.set_ylim(0.0, gam)  # Use one physical state range for all six panels.
        ax.set_xlabel(r"$\alpha$")  # Keep alpha visibly identified as the control parameter.
        ax.set_ylabel(r"$a$")  # Label the continued physical state.
        ax.grid(alpha=0.25)  # Add a light common reading grid.
        if index == 0:  # Explain line styles once without repeating six legends.
            ax.legend(fontsize=7, loc="center right")  # Place the shared legend in the baseline panel.

    fig.suptitle(rf"Progressive unfolding at six small fixed biases "
                 rf"($\delta={delta:g}$; $\alpha$ is the control parameter)", fontsize=14)  # Main result.
    fig.savefig(path, dpi=220)  # Save the comparison for the report and presentation.
    plt.close(fig)  # Release the Matplotlib resources.


def plot_small_negative_bias_alpha_bifurcation(path: Path, delta: float) -> None:
    """Compare six small nonpositive biases with alpha as the control parameter."""
    gam = gamma(delta)  # Compute the physical scale that maps x to a.
    alpha_min = 0.25  # Start below the symmetric critical value.
    alpha_max = 5.0  # Continue far enough to display all asymptotic branches.
    h_values = [0.0, -0.002, -0.005, -0.010, -0.020, -0.040]  # Reflect the positive sequence.
    stable_color = "#1857c9"  # Use one color for attracting fixed-point branches.
    unstable_color = "#c23b23"  # Use a contrasting color for repelling branches.
    fig, axes = plt.subplots(2, 3, figsize=(12.0, 7.6), sharex=True, sharey=True,
                             constrained_layout=True)  # Build six matched continuation panels.

    for index, (ax, h_small) in enumerate(zip(axes.flat, h_values)):  # Visit each selected bias.
        if h_small == 0.0:  # Draw the exact reflection-symmetric pitchfork as the baseline.
            alpha_c = alpha_critical(delta)  # Read the universal symmetric threshold.
            central_stable = np.linspace(alpha_min, alpha_c, 400)  # Sample stable center.
            central_unstable = np.linspace(alpha_c, alpha_max, 600)  # Sample unstable center.
            x_outer = np.linspace(1.0e-4, 0.995, 2200)  # Parameterize both outer branches.
            alpha_outer = 2.0 * np.arctanh(x_outer) / x_outer  # Recover alpha on the branches.
            outer_mask = alpha_outer <= alpha_max  # Retain the displayed part of each branch.
            ax.plot(central_stable, np.full_like(central_stable, 0.5 * gam),
                    color=stable_color, linewidth=1.8, label="attracting")  # Stable center.
            ax.plot(central_unstable, np.full_like(central_unstable, 0.5 * gam),
                    color=unstable_color, linestyle="--", linewidth=1.8,
                    label="repelling")  # Unstable center.
            ax.plot(alpha_outer[outer_mask], 0.5 * gam * (1.0 + x_outer[outer_mask]),
                    color=stable_color, linewidth=1.8)  # Stable upper branch.
            ax.plot(alpha_outer[outer_mask], 0.5 * gam * (1.0 - x_outer[outer_mask]),
                    color=stable_color, linewidth=1.8)  # Stable lower branch.
            ax.plot(alpha_c, 0.5 * gam, "ko", markersize=4.5,
                    label="critical point")  # Mark the perfect pitchfork.
            ax.set_title(r"$h=0$: perfect pitchfork")  # Identify the symmetric reference.
        else:  # Continue the imperfect branches at one negative fixed bias.
            x_negative = np.linspace(-0.995, np.tanh(h_small) - 1.0e-6, 3000)  # Favored branch.
            x_positive = np.linspace(1.0e-4, 0.995, 5000)  # Sample the upper branch pair.
            alpha_negative = 2.0 * (np.arctanh(x_negative) - h_small) / x_negative  # Continue roots.
            alpha_positive = 2.0 * (np.arctanh(x_positive) - h_small) / x_positive  # Continue roots.
            lambda_positive = 0.5 * alpha_positive * (1.0 - x_positive * x_positive)  # Slopes.
            negative_visible = (alpha_negative >= alpha_min) & (alpha_negative <= alpha_max)  # Clip.
            positive_visible = (alpha_positive >= alpha_min) & (alpha_positive <= alpha_max)  # Clip.
            positive_stable = positive_visible & (lambda_positive < 1.0)  # Stable upper branch.
            positive_unstable = positive_visible & (lambda_positive >= 1.0)  # Repelling branch.
            ax.plot(alpha_negative[negative_visible],
                    0.5 * gam * (1.0 + x_negative[negative_visible]),
                    color=stable_color, linewidth=1.8, label="attracting")  # Favored lower branch.
            ax.plot(alpha_positive[positive_stable],
                    0.5 * gam * (1.0 + x_positive[positive_stable]),
                    color=stable_color, linewidth=1.8)  # Stable branch born at the fold.
            ax.plot(alpha_positive[positive_unstable],
                    0.5 * gam * (1.0 + x_positive[positive_unstable]),
                    color=unstable_color, linestyle="--", linewidth=1.8,
                    label="repelling")  # Unstable partner branch.
            x_fold, alpha_fold = alpha_fold_for_bias(h_small)  # Locate the reflected saddle node.
            a_fold = 0.5 * gam * (1.0 + x_fold)  # Convert the fold state to physical a.
            ax.plot(alpha_fold, a_fold, "ko", markersize=4.5,
                    label="saddle node")  # Mark the shifted multiplicity threshold.
            ax.set_title(rf"$h={h_small:.3f}$, $\alpha_{{\rm SN}}={alpha_fold:.3f}$")  # Show shift.
        ax.set_xlim(alpha_min, alpha_max)  # Use one alpha range for all six panels.
        ax.set_ylim(0.0, gam)  # Use one physical state range for all six panels.
        ax.set_xlabel(r"$\alpha$")  # Keep alpha visibly identified as the control parameter.
        ax.set_ylabel(r"$a$")  # Label the continued physical state.
        ax.grid(alpha=0.25)  # Add a light common reading grid.
        if index == 0:  # Explain line styles once without repeating six legends.
            ax.legend(fontsize=7, loc="center right")  # Place the shared legend in the baseline panel.

    fig.suptitle(rf"Mirror unfolding at six small nonpositive biases "
                 rf"($\delta={delta:g}$; $\alpha$ is the control parameter)", fontsize=14)  # Main result.
    fig.savefig(path, dpi=220)  # Save the reflected comparison for direct inspection.
    plt.close(fig)  # Release the Matplotlib resources.


def plot_critical_point_trajectory(path: Path, delta: float) -> None:
    """Trace the saddle-node point that emanates from the pitchfork limit."""
    gam = gamma(delta)  # Compute the physical scale that maps x to a.
    x_fold = np.linspace(-0.45, 0.45, 2401)  # Parameterize both sides of the local cusp.
    h_exact = np.arctanh(x_fold) - x_fold / (1.0 - x_fold * x_fold)  # Fold bias.
    alpha_exact = 2.0 / (1.0 - x_fold * x_fold)  # Fold alpha from the tangency condition.
    a_exact = 0.5 * gam * (1.0 + x_fold)  # Convert each normalized fold state to physical a.
    order = np.argsort(h_exact)  # Sort the parametric curve from negative to positive bias.
    h_sorted = h_exact[order]  # Reorder the horizontal continuation coordinate.
    alpha_sorted = alpha_exact[order]  # Reorder the exact saddle-node threshold.
    a_sorted = a_exact[order]  # Reorder the corresponding physical fold state.
    alpha_asymptotic = 2.0 + 2.0 * (1.5 * np.abs(h_sorted)) ** (2.0 / 3.0)  # Local cusp law.
    h_samples_positive = np.array([0.002, 0.005, 0.010, 0.020, 0.040])  # Match Figure 5.
    h_samples = np.concatenate((-h_samples_positive[::-1], [0.0], h_samples_positive))  # Mirror them.
    x_samples = []  # Store the exact fold state for every highlighted bias.
    alpha_samples = []  # Store the exact threshold for every highlighted bias.
    for h_value in h_samples:  # Evaluate all positive, negative, and zero sample points.
        x_value, alpha_value = alpha_fold_for_bias(float(h_value))  # Solve the fold equations.
        x_samples.append(x_value)  # Record the normalized saddle-node state.
        alpha_samples.append(alpha_value)  # Record the corresponding alpha threshold.
    x_samples = np.asarray(x_samples)  # Convert the state list to a plotting array.
    alpha_samples = np.asarray(alpha_samples)  # Convert threshold values to a plotting array.
    a_samples = 0.5 * gam * (1.0 + x_samples)  # Convert sampled fold states to physical a.

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.8), constrained_layout=True)  # Build two projections.
    axes[0].plot(h_sorted, alpha_sorted, color="#1857c9", linewidth=2.2,
                 label="exact fold trajectory")  # Draw the exact cusp boundary.
    axes[0].plot(h_sorted, alpha_asymptotic, color="#c23b23", linestyle="--",
                 linewidth=1.8, label=r"$2+2(3|h|/2)^{2/3}$")  # Draw the local approximation.
    axes[0].plot(h_samples, alpha_samples, "ko", markersize=4.5,
                 label="sampled biases")  # Mark the values used in the six-panel figure.
    axes[0].plot(0.0, 2.0, "o", color="#d627b5", markersize=6,
                 label="pitchfork limit")  # Mark the cusp tip at the symmetric problem.
    axes[0].set(xlabel=r"fixed bias $h$", ylabel=r"$\alpha_{\rm SN}$",
                title="Multiplicity threshold in parameter space")  # Label the first projection.
    axes[0].legend(fontsize=8, loc="upper center")  # Identify exact, asymptotic, and sample data.
    axes[0].grid(alpha=0.25)  # Add a light reading grid.

    axes[1].plot(h_sorted, a_sorted, color="#1857c9", linewidth=2.2,
                 label="exact fold state")  # Draw the state coordinate of the moving fold.
    axes[1].plot(h_samples, a_samples, "ko", markersize=4.5,
                 label="sampled biases")  # Mark the same selected biases in state space.
    axes[1].plot(0.0, 0.5 * gam, "o", color="#d627b5", markersize=6,
                 label="pitchfork limit")  # Mark the symmetric critical state.
    axes[1].axhline(0.5 * gam, color="#777777", linestyle=":", linewidth=1.0)  # Show symmetry center.
    axes[1].text(-0.075, 0.725 * gam, r"$h<0$: upper saddle node")  # Label the upper trajectory.
    axes[1].text(0.018, 0.285 * gam, r"$h>0$: lower saddle node")  # Label the lower trajectory.
    axes[1].set(xlabel=r"fixed bias $h$", ylabel=r"$a_{\rm SN}$",
                title="Critical state moves away from the center")  # Label the state projection.
    axes[1].legend(fontsize=8, loc="center right")  # Identify the exact path and highlighted points.
    axes[1].grid(alpha=0.25)  # Add a light reading grid.

    fig.suptitle(rf"Saddle-node trajectory from the pitchfork limit ($\delta={delta:g}$)",
                 fontsize=14)  # State the geometric interpretation concisely.
    fig.savefig(path, dpi=220)  # Save the trajectory figure for documents and comparison.
    plt.close(fig)  # Release the Matplotlib resources.


def plot_s_curve_fixed_point_cobweb(path: Path, delta: float, alpha: float) -> None:
    """Show the three biased fixed points and their two attracting cobweb basins."""
    folds = fold_data(alpha, delta)  # Compute the fold data for the selected supercritical alpha.
    if folds is None:  # Reject parameters that cannot support a three-root biased state.
        raise ValueError("The S-curve cobweb requires alpha > 2")  # Explain the requirement.
    h_sample = 0.5 * folds["h_limit"]  # Select a bias strictly inside the three-root region.
    gam = gamma(delta)  # Compute the physical interval length for a.
    a_grid = np.linspace(0.0, gam, 2000)  # Sample the biased map smoothly.
    states = physical_states(alpha, delta, h_sample)  # Compute all three biased fixed points.
    middle_fixed_point = states[1]["a"]  # Read the repelling fixed point that separates the basins.
    initial_offset = 0.08 * gam  # Stay close enough to the separator to reveal the transient paths.
    initial_values = [  # Place one initial condition on each side of the basin boundary.
        middle_fixed_point - initial_offset,  # Start in the lower attracting basin.
        middle_fixed_point + initial_offset,  # Start in the upper attracting basin.
    ]
    path_colors = ["#1b8f4d", "#c43b2b"]  # Assign colors to the lower and upper paths.
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.8), constrained_layout=True)  # Build two panels.

    axes[0].plot(a_grid, f_map(a_grid, alpha, delta, h_sample), color="#1857c9",
                 linewidth=2.0, label=r"$F_{\alpha,h}(a)$")  # Draw the biased map.
    axes[0].plot(a_grid, a_grid, color="#d62728", linestyle="--",
                 linewidth=1.4, label=r"$a$")  # Draw the identity line.
    axes[0].plot([state["a"] for state in states], [state["a"] for state in states],
                 "ko", markersize=6, label="Fixed points")  # Mark all three intersections.
    axes[0].set(xlim=(0.0, gam), ylim=(0.0, gam), xlabel=r"$a$",
                ylabel=r"$F_{\alpha,h}(a)$", title="Three fixed-point intersections")  # Label it.
    axes[0].legend(fontsize=8, loc="upper left")  # Identify the geometry elements.
    axes[0].grid(alpha=0.25)  # Add a light grid.

    axes[1].plot(a_grid, f_map(a_grid, alpha, delta, h_sample), color="#1857c9",
                 linewidth=1.8, label=r"$F_{\alpha,h}(a)$")  # Draw the same biased map.
    axes[1].plot(a_grid, a_grid, "k--", linewidth=1.2, label=r"$a$")  # Draw the identity line.
    for a0, color, label in zip(initial_values, path_colors,
                                ["lower start", "upper start"]):  # Draw both basins.
        xs, ys = cobweb_segments(a0, alpha, delta, iterations=35, h=h_sample)  # Build one path.
        axes[1].plot(xs, ys, color=color, linewidth=0.8, label=label)  # Draw it.
    axes[1].plot([state["a"] for state in states], [state["a"] for state in states],
                 "ko", markersize=5, label="Fixed points")  # Mark all fixed points.
    axes[1].set(xlim=(0.0, gam), ylim=(0.0, gam), xlabel=r"$a_n$",
                ylabel=r"$a_{n+1}$", title="Two attracting basins")  # Label the dynamics.
    axes[1].legend(fontsize=7, loc="upper left")  # Identify paths and reference curves.
    axes[1].grid(alpha=0.25)  # Add a light grid.

    fig.suptitle(rf"Fixed-$h$ slice of the biased model: $\delta={delta:g}$, "
                 rf"$\alpha={alpha:g}$, $h={h_sample:.4f}$", fontsize=14)  # Identify h as fixed.
    fig.savefig(path, dpi=220)  # Save the fixed-point and cobweb figure.
    plt.close(fig)  # Release the Matplotlib resources.


def plot_fixed_a_h_iteration(path: Path, delta: float, alpha: float) -> None:
    """Recover the bias h that makes one prescribed state a_star steady."""
    folds = fold_data(alpha, delta)  # Compute the focused h range and a three-root sample.
    if folds is None:  # Reject a parameter value without the biased three-root regime.
        raise ValueError("The fixed-a h iteration requires alpha > 2")  # Explain the requirement.
    h_star = 0.5 * folds["h_limit"]  # Reuse the bias selected in the S-curve examples.
    target_state = physical_states(alpha, delta, h_star)[-1]  # Select the upper attracting state.
    x_star = target_state["x"]  # Read its normalized state, which is held fixed here.
    a_star = target_state["a"]  # Record the corresponding physical affine parameter.
    kappa = alpha / alpha_critical(delta)  # Convert to the normalized sigmoid slope.
    h_grid = np.linspace(-2.2 * folds["h_limit"], 2.2 * folds["h_limit"], 2000)  # Focus on folds.
    residual = np.tanh(kappa * x_star + h_grid) - x_star  # Test whether a_star is steady.
    recovery_map = h_grid + x_star - np.tanh(kappa * x_star + h_grid)  # Relax the residual.
    initial_values = [-1.6 * folds["h_limit"], 1.8 * folds["h_limit"]]  # Start on both sides.
    path_colors = ["#1b8f4d", "#c43b2b"]  # Distinguish the lower and upper h starts.
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.8), constrained_layout=True)  # Build two panels.

    axes[0].plot(h_grid, residual, color="#1857c9", linewidth=2.0,
                 label=r"$R_{a_*}(h)$")  # Draw the scalar parameter-recovery residual.
    axes[0].axhline(0.0, color="#d62728", linestyle="--", linewidth=1.4,
                    label=r"$R_{a_*}=0$")  # Draw the root condition.
    axes[0].plot(h_star, 0.0, "ko", markersize=6, label=r"$h_*$")  # Mark the recovered bias.
    axes[0].set(xlabel=r"$h$", ylabel=r"$R_{a_*}(h)$",
                title="A fixed state determines one bias")  # Label the residual panel.
    axes[0].legend(fontsize=8, loc="upper left")  # Identify the residual and its root.
    axes[0].grid(alpha=0.25)  # Add a light grid.

    axes[1].plot(h_grid, recovery_map, color="#1857c9", linewidth=1.8,
                 label=r"$\mathcal{H}_{a_*}(h)$")  # Draw the relaxed h iteration map.
    axes[1].plot(h_grid, h_grid, "k--", linewidth=1.2, label=r"$h$")  # Draw the identity line.
    for h0, color, label in zip(initial_values, path_colors,
                                ["lower start", "upper start"]):  # Draw both h iterations.
        path_x = [h0]  # Start the cobweb at the selected bias guess.
        path_y = [h0]  # Place the initial point on the identity line.
        h_current = h0  # Initialize the relaxed parameter-recovery iteration.
        for _ in range(35):  # Generate enough steps to display convergence clearly.
            h_next = h_current + x_star - np.tanh(kappa * x_star + h_current)  # Update h.
            path_x.extend([h_current, h_next])  # Add vertical and horizontal coordinates.
            path_y.extend([h_next, h_next])  # Complete one cobweb step on the map and identity.
            h_current = h_next  # Advance to the next recovered-bias iterate.
        axes[1].plot(path_x, path_y, color=color, linewidth=0.8, label=label)  # Draw one path.
    axes[1].plot(h_star, h_star, "ko", markersize=5, label=r"$h_*$")  # Mark the fixed bias.
    axes[1].set(xlabel=r"$h_n$", ylabel=r"$h_{n+1}$",
                title="Relaxed iteration recovers the bias")  # Label the cobweb panel.
    axes[1].legend(fontsize=7, loc="upper left")  # Identify the map, starts, and recovered bias.
    axes[1].grid(alpha=0.25)  # Add a light grid.

    fig.suptitle(rf"Fixed-$a$ parameter recovery: $\delta={delta:g}$, $\alpha={alpha:g}$, "
                 rf"$a_*={a_star:.6f}$, $h_*={h_star:.6f}$", fontsize=14)  # State the target.
    fig.savefig(path, dpi=220)  # Save the fixed-a parameter-recovery figure.
    plt.close(fig)  # Release the Matplotlib resources.


def write_s_curve_data(path: Path, alpha: float, delta: float) -> None:
    """Write fold points and representative imperfect-model roots to CSV."""
    folds = fold_data(alpha, delta)  # Compute analytic fold information if it exists.
    lines = ["kind,h,x,a,s,beta,q,iteration_slope,iteration_attractive"]  # Add the header.
    if folds is None:  # Handle a subcritical case without an S-curve.
        sample_h_values = [0.0]  # Record only the symmetric state set.
    else:  # Process the two-fold S-curve.
        x_fold = folds["x_fold"]  # Read the positive fold state.
        h_limit = folds["h_limit"]  # Read the positive critical bias.
        gam = gamma(delta)  # Compute the physical scale for a.
        for h_value, x_value in [(h_limit, -x_fold), (-h_limit, x_fold)]:  # Visit both folds.
            a = 0.5 * gam * (1.0 + x_value)  # Reconstruct a at the fold.
            s = 0.5 * (1.0 - x_value)  # Reconstruct the jump at the fold.
            lines.append(
                f"fold,{h_value:.12g},{x_value:.12g},{a:.12g},{s:.12g},"
                f"{beta(s, alpha, delta, h_value):.12g},{delta * a:.12g},"
                "1,False"
            )  # Store the fold, whose iteration slope is exactly one.
        sample_h_values = [-1.2 * h_limit, -0.5 * h_limit, 0.0,
                           0.5 * h_limit, 1.2 * h_limit]  # Sample one- and three-root regimes.

    for h_value in sample_h_values:  # Loop over every selected interface-law amplitude.
        for state in physical_states(alpha, delta, h_value):  # Loop over all roots at this bias.
            attractive = abs(state["iteration_slope"]) < 1.0  # Classify fixed-point attraction.
            lines.append(
                f"state,{h_value:.12g},{state['x']:.12g},{state['a']:.12g},"
                f"{state['s']:.12g},{state['beta']:.12g},{state['q']:.12g},"
                f"{state['iteration_slope']:.12g},{attractive}"
            )  # Serialize one imperfect-model state.
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")  # Save all records.


def plot_s_curve(path: Path, delta: float, alpha: float) -> None:
    """Create the four-panel summary of the imperfect S-shaped bifurcation."""
    folds = fold_data(alpha, delta)  # Compute the two analytic folds.
    if folds is None:  # Reject parameters that cannot produce an S-curve.
        raise ValueError("The S-curve requires alpha > alpha_c")  # Explain the condition.

    kappa = folds["kappa"]  # Read the normalized map slope.
    x_fold = folds["x_fold"]  # Read the positive fold state.
    h_limit = folds["h_limit"]  # Read the positive critical bias.
    gam = gamma(delta)  # Compute the physical a interval.
    x_grid = np.linspace(-0.995, 0.995, 2400)  # Parameterize the S-curve away from singular endpoints.
    h_grid = np.arctanh(x_grid) - kappa * x_grid  # Evaluate the exact continuation curve h(x).
    lower_stable = x_grid < -x_fold  # Identify the lower attracting outer branch.
    middle_unstable = np.abs(x_grid) <= x_fold  # Identify the repelling middle branch.
    upper_stable = x_grid > x_fold  # Identify the upper attracting outer branch.
    h_sample = 0.5 * h_limit  # Choose a bias strictly inside the three-root region.
    sample_states = physical_states(alpha, delta, h_sample)  # Compute its three states.

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.2), constrained_layout=True)  # Build the panel grid.

    ax = axes[0, 0]  # Select the continuation-curve panel.
    ax.plot(h_grid[lower_stable], x_grid[lower_stable], color="#3478a8",
            linewidth=2.0, label=r"$|\partial_x\tanh(\kappa x+h)|<1$")  # Plot the lower stable branch.
    ax.plot(h_grid[upper_stable], x_grid[upper_stable], color="#3478a8",
            linewidth=2.0)  # Plot the upper stable branch.
    ax.plot(h_grid[middle_unstable], x_grid[middle_unstable],
            color="#c23b23", linestyle="--", linewidth=2.0,
            label="repelling branch")  # Plot the unstable middle branch.
    ax.plot([h_limit, -h_limit], [-x_fold, x_fold], "ko", label="folds")  # Mark both folds.
    ax.axvline(h_sample, color="#777777", linestyle="--", linewidth=1.0)  # Mark the sample bias.
    for state in sample_states:  # Mark all sample roots on the continuation curve.
        ax.plot(h_sample, state["x"], "o", color="#6b3fa0", markersize=5)  # Add one root marker.
    ax.set(
        xlabel=r"$h$",
        ylabel=r"$x$",
        title=rf"S-curve, $\kappa={kappa:.2f}$",
        xlim=(-2.2 * h_limit, 2.2 * h_limit),
    )  # Label and bound the continuation panel.
    ax.legend(fontsize=7)  # Identify branch types and folds.
    ax.grid(alpha=0.25)  # Add a light reading grid.

    a_grid = np.linspace(0.0, gam, 1000)  # Sample the physical interval for a.
    ax = axes[0, 1]  # Select the biased fixed-point-map panel.
    ax.plot(a_grid, a_grid, "k--", linewidth=1.1, label=r"$y=a$")  # Draw the identity line.
    for h_value, color in [
        (-1.2 * h_limit, "#3478a8"),
        (h_sample, "#6b3fa0"),
        (1.2 * h_limit, "#c23b23"),
    ]:  # Compare one-root and three-root biases.
        ax.plot(
            a_grid,
            f_map(a_grid, alpha, delta, h_value),
            color=color,
            label=rf"$h={h_value:.3f}$",
        )  # Draw one biased fixed-point map.
    ax.set(xlabel=r"$a$", ylabel=r"$f_{\alpha,h}(a)$",
           title="One-root and three-root regimes")  # Label the map panel.
    ax.legend(fontsize=7)  # Identify the three bias values.
    ax.grid(alpha=0.25)  # Add a light reading grid.

    s_grid = np.linspace(0.0, 1.0, 1000)  # Sample the full jump interval.
    ax = axes[1, 0]  # Select the conductance-amplitude panel.
    for h_value, color in [
        (-h_limit, "#3478a8"),
        (0.0, "#777777"),
        (h_limit, "#c23b23"),
    ]:  # Compare the law at the two folds and at h=0.
        ax.plot(
            s_grid,
            beta(s_grid, alpha, delta, h_value),
            color=color,
            label=rf"$h={h_value:.3f}$",
        )  # Draw one scaled Kapitza conductance law.
    ax.set(
        xlabel=r"$s=[u]$",
        ylabel=r"$\beta_h(s)$",
        title=r"$\beta_h=e^{2h}\beta_0$",
    )  # Label the nonlinear-law panel.
    ax.legend(fontsize=7)  # Identify all conductance amplitudes.
    ax.grid(alpha=0.25)  # Add a light reading grid.

    ax = axes[1, 1]  # Select the imperfect heat-flux panel.
    q_k = beta(s_grid, alpha, delta, h_sample) * s_grid  # Compute the biased interface flux.
    q_bulk = delta / (1.0 + delta) * (1.0 - s_grid)  # Compute the unchanged bulk flux.
    ax.plot(s_grid, q_k, color="#c23b23", label=r"$q_K=\beta_h(s)s$")  # Draw the interface flux.
    ax.plot(s_grid, q_bulk, color="#3478a8", label=r"$q_{\rm bulk}$")  # Draw the bulk flux.
    for state in sample_states:  # Mark all three heat-flux intersections.
        ax.plot(state["s"], state["q"], "ko", markersize=4)  # Add one root marker.
    ax.set(
        xlabel=r"$s=[u]$",
        ylabel=r"$q$",
        title=rf"Three flux intersections, $h={h_sample:.3f}$",
    )  # Label the biased flux panel.
    ax.legend(fontsize=7)  # Identify the two heat-flux curves.
    ax.grid(alpha=0.25)  # Add a light reading grid.

    fig.suptitle(
        rf"Imperfect Kapitza bifurcation: $\delta={delta:g}$, $\alpha={alpha:g}$",
        fontsize=14,
    )  # Add a title describing the selected parameters.
    fig.savefig(path, dpi=220)  # Save the complete S-curve summary.
    plt.close(fig)  # Release the figure from memory.


def main() -> None:
    """Parse inputs, run all computations, and save every requested output."""
    parser = argparse.ArgumentParser()  # Create the command-line parser.
    parser.add_argument(
        "--delta", "--d", dest="delta", type=float, default=0.5
    )  # Use paper notation while accepting the teacher's shorthand --d as an alias.
    parser.add_argument("--alpha", type=float, default=None)  # Read alpha in the teacher's convention.
    parser.add_argument(
        "--output-dir", type=Path, default=Path("output/python")
    )  # Keep Python-generated figures and data in their dedicated output folder.
    args = parser.parse_args()  # Convert command-line text into typed values.

    if args.delta <= 0.0:  # Enforce the physical sign assumption on delta.
        raise ValueError("delta must be positive")  # Stop with a clear error message.
    alpha = (
        args.alpha
        if args.alpha is not None
        else 1.5 * alpha_critical(args.delta)
    )  # Use a supercritical default when the user omits alpha.
    if alpha <= 0.0:  # Enforce a positive sigmoid steepness.
        raise ValueError("alpha must be positive")  # Stop with a clear error message.

    args.output_dir.mkdir(parents=True, exist_ok=True)  # Create missing output folders safely.
    write_states(args.output_dir / "kapitza_states.csv", alpha, args.delta)  # Save symmetric roots.
    plot_summary(args.output_dir / "kapitza_summary.png", args.delta, alpha)  # Save the summary plot.
    plot_cobweb(args.output_dir / "kapitza_cobweb.png", args.delta, alpha)  # Save the iteration plot.
    plot_fixed_point_regimes(
        args.output_dir / "fixed_point_regimes.png", args.delta
    )  # Save the one-critical-three fixed-point comparison.
    plot_cobweb_regimes(
        args.output_dir / "cobweb_regimes.png", args.delta
    )  # Save the matching cobweb comparison.
    plot_small_bias_alpha_bifurcation(
        args.output_dir / "small_bias_alpha_bifurcation.png", args.delta
    )  # Show how a small fixed bias unfolds the alpha-controlled pitchfork.
    plot_small_negative_bias_alpha_bifurcation(
        args.output_dir / "small_negative_bias_alpha_bifurcation.png", args.delta
    )  # Show the reflected unfolding for the matching negative bias values.
    plot_critical_point_trajectory(
        args.output_dir / "critical_point_trajectory.png", args.delta
    )  # Trace the exact saddle-node path and its local cusp approximation.
    if alpha > alpha_critical(args.delta):  # Generate S-curve outputs only above threshold.
        write_s_curve_data(
            args.output_dir / "kapitza_s_curve_states.csv", alpha, args.delta
        )  # Save fold and imperfect-branch data.
        plot_s_curve(args.output_dir / "kapitza_s_curve.png", args.delta, alpha)  # Save the S-curve plot.
        plot_s_curve_fixed_point_cobweb(
            args.output_dir / "s_curve_fixed_point_cobweb.png", args.delta, alpha
        )  # Save biased fixed-point intersections and cobweb basins.
        plot_fixed_a_h_iteration(
            args.output_dir / "fixed_a_h_iteration.png", args.delta, alpha
        )  # Save the fixed-state parameter-recovery residual and h cobweb.

    print(f"delta={args.delta:g}")  # Report the selected delta value.
    print(f"alpha_c={alpha_critical(args.delta):.12g}")  # Report the analytic threshold.
    print(f"alpha={alpha:.12g}")  # Report the selected alpha value.
    folds = fold_data(alpha, args.delta)  # Compute folds for terminal reporting.
    if folds is not None:  # Print fold coordinates only when they exist.
        print(
            f"S-curve folds: (h,x)=({folds['h_limit']:.12g},"
            f"{-folds['x_fold']:.12g}) and "
            f"({-folds['h_limit']:.12g},{folds['x_fold']:.12g})"
        )  # Display both analytic fold points.
    for index, state in enumerate(physical_states(alpha, args.delta), start=1):  # Report each root.
        print(
            f"state {index}: a={state['a']:.12g}, s={state['s']:.12g}, "
            f"beta={state['beta']:.12g}, q={state['q']:.12g}, "
            f"f'={state['iteration_slope']:.12g}"
        )  # Display all reconstructed physical quantities.


if __name__ == "__main__":  # Run the workflow only when this file is executed directly.
    main()  # Enter the command-line program.
