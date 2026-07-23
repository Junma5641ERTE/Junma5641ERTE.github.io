#!/usr/bin/env python3
"""Study the source-controlled one-dimensional nonlinear Kapitza model.

The script is intentionally self-contained and extensively commented in
English.  It reconstructs the one-dimensional temperature fields, computes
all fixed points, classifies their Picard-iteration stability, and produces
the figures and CSV tables used by the accompanying report and Beamer slides.
"""

from __future__ import annotations  # Enable modern annotations on older Python releases.

import argparse  # Read optional model parameters from the command line.
import csv  # Write transparent, software-independent numerical tables.
from pathlib import Path  # Construct output paths without platform-specific separators.

import matplotlib  # Configure Matplotlib before importing its pyplot interface.
import numpy as np  # Evaluate vectorized model formulas and parameterizations.
from scipy.optimize import brentq  # Isolate scalar roots with a robust bracketed method.

matplotlib.use("Agg")  # Select a non-interactive backend for reproducible batch execution.
import matplotlib.pyplot as plt  # Create publication-quality static figures.

plt.rcParams.update(  # Apply one typography policy to every generated figure.
    {
        "font.family": "serif",  # Use a serif family for prose, ticks, and legends.
        "font.serif": ["Times New Roman", "Times", "STIXGeneral"],  # Prefer Times.
        "mathtext.fontset": "stix",  # Use Times-compatible mathematical glyphs.
        "font.size": 11,  # Keep labels readable in reports and presentation slides.
        "axes.titlesize": 12,  # Use compact subplot titles.
        "axes.labelsize": 11,  # Match axis-label size to the body text.
        "legend.fontsize": 9,  # Keep legends informative without obscuring curves.
        "axes.unicode_minus": False,  # Render minus signs consistently with the serif font.
        "savefig.bbox": "tight",  # Remove unnecessary margins around saved figures.
        "pdf.fonttype": 42,  # Embed TrueType-compatible fonts in vector PDF figures.
        "ps.fonttype": 42,  # Preserve editable text if figures are converted later.
    }
)  # Finish the global plotting configuration.

BLUE = "#1f5fbf"  # Encode attracting branches and primary maps in blue.
RED = "#c43d2b"  # Encode repelling branches and upper cobweb starts in red.
GREEN = "#278c5a"  # Encode lower cobweb starts in green.
BLACK = "#111111"  # Encode identities, roots, and critical points in black.
GRAY = "#6b6b6b"  # Encode reference curves and secondary annotations in gray.


def gamma(kappa: float, delta: float) -> float:
    """Return the upper physical bound gamma(kappa)=kappa/(1+delta)."""
    return kappa / (1.0 + delta)  # Evaluate the moving physical interval exactly.


def beta_hat(
    s: np.ndarray | float, kappa: float, delta: float, h: float = 0.0
) -> np.ndarray | float:
    """Return the interface conductance normalized by the right conductivity mu."""
    exponent = 2.0 * (kappa - 2.0 * s) + 2.0 * h  # Build the biased exponential argument.
    return delta / (1.0 + delta) * np.exp(exponent)  # Apply the normalized prefactor.


def kapitza_map(
    a: np.ndarray | float, kappa: float, delta: float, h: float = 0.0
) -> np.ndarray | float:
    """Evaluate the physical fixed-point map a -> F_{kappa,h}(a)."""
    s = kappa - (1.0 + delta) * a  # Reconstruct the interface jump from the bulk solution.
    conductance = beta_hat(s, kappa, delta, h)  # Evaluate the nonlinear interface law.
    denominator = delta + (1.0 + delta) * conductance  # Assemble the interface balance.
    return kappa * conductance / denominator  # Solve the balance for the next state a.


def normalized_map(x: np.ndarray | float, kappa: float, h: float = 0.0):
    """Evaluate the equivalent normalized map x -> tanh(kappa*x+h)."""
    return np.tanh(kappa * x + h)  # Apply the exact reduced sigmoid map.


def normalized_roots(kappa: float, h: float = 0.0) -> list[float]:
    """Find every root of x=tanh(kappa*x+h) in the physical interval (-1,1)."""
    if abs(h) < 1.0e-14:  # Exploit exact odd symmetry for the perfect problem.
        if kappa <= 1.0:  # The central state is the unique root through the threshold.
            return [0.0]  # Avoid misclassifying the flat critical residual numerically.
        positive = brentq(lambda x: float(normalized_map(x, kappa, 0.0) - x), 1.0e-10, 1.0 - 1.0e-10)  # Isolate the outer root.
        return [-float(positive), 0.0, float(positive)]  # Restore the exact symmetric triplet.
    residual = lambda x: float(normalized_map(x, kappa, h) - x)  # Define the root residual.
    grid = np.linspace(-1.0 + 1.0e-10, 1.0 - 1.0e-10, 20001)  # Scan the open interval finely.
    values = normalized_map(grid, kappa, h) - grid  # Evaluate all residual samples at once.
    roots: list[float] = []  # Collect sign-changing and nearly tangent roots.
    for index in range(grid.size - 1):  # Inspect every neighboring pair in the scan.
        left = float(grid[index])  # Extract the left endpoint as a scalar.
        right = float(grid[index + 1])  # Extract the right endpoint as a scalar.
        f_left = float(values[index])  # Read the residual at the left endpoint.
        f_right = float(values[index + 1])  # Read the residual at the right endpoint.
        if f_left * f_right < 0.0:  # Detect an ordinary simple root by a sign change.
            roots.append(float(brentq(residual, left, right)))  # Refine it to high precision.
    roots.sort()  # Put all candidate roots in increasing state order.
    unique: list[float] = []  # Prepare a duplicate-free root list.
    for root in roots:  # Compare each candidate with the previous accepted root.
        if not unique or abs(root - unique[-1]) > 2.0e-6:  # Reject repeated scan detections.
            unique.append(root)  # Keep a genuinely distinct root.
    return unique  # Return all physical normalized steady states.


def physical_states(kappa: float, delta: float, h: float = 0.0) -> list[dict[str, float]]:
    """Reconstruct physical variables and stability data for every fixed point."""
    states: list[dict[str, float]] = []  # Store one dictionary per steady solution.
    scale = gamma(kappa, delta)  # Compute the current physical state interval.
    for x in normalized_roots(kappa, h):  # Convert every normalized solution separately.
        a = 0.5 * scale * (1.0 + x)  # Recover the left interface temperature u_I(0).
        s = 0.5 * kappa * (1.0 - x)  # Recover the temperature jump [u].
        b = delta * a  # Recover the right affine coefficient from flux continuity.
        conductance = float(beta_hat(s, kappa, delta, h))  # Evaluate normalized beta.
        q_hat = delta * a  # Evaluate the heat flux normalized by mu.
        multiplier = kappa * (1.0 - x * x)  # Differentiate the normalized Picard map.
        states.append(  # Package all quantities needed by tables and plots.
            {
                "x": float(x),  # Store the normalized fixed point.
                "a": float(a),  # Store the left interface trace.
                "b": float(b),  # Store the right affine coefficient.
                "s": float(s),  # Store the temperature jump.
                "beta_hat": conductance,  # Store the normalized conductance.
                "q_hat": float(q_hat),  # Store the normalized common heat flux.
                "multiplier": float(multiplier),  # Store the Picard multiplier.
                "attracting": abs(multiplier) < 1.0,  # Apply the local iteration criterion.
            }
        )  # Finish one physical-state record.
    return states  # Return every reconstructed steady state.


def fold_for_bias(h: float) -> tuple[float, float]:
    """Return the saddle-node state x_SN and source parameter kappa_SN."""
    if abs(h) < 1.0e-15:  # Handle the perfect pitchfork limit analytically.
        return 0.0, 1.0  # Return the symmetric critical point.
    residual = lambda x: np.arctanh(x) - x / (1.0 - x * x) - h  # Combine fold equations.
    if h > 0.0:  # A positive bias produces the fold on the negative-state side.
        x_fold = brentq(residual, -0.999999, -1.0e-12)  # Isolate the unique negative root.
    else:  # A negative bias produces the reflected positive-state fold.
        x_fold = brentq(residual, 1.0e-12, 0.999999)  # Isolate the unique positive root.
    kappa_fold = 1.0 / (1.0 - x_fold * x_fold)  # Apply the tangency condition exactly.
    return float(x_fold), float(kappa_fold)  # Return the fold state and control parameter.


def fold_data_at_kappa(kappa: float) -> dict[str, float] | None:
    """Return the two fixed-kappa S-curve folds when kappa exceeds one."""
    if kappa <= 1.0:  # A subcritical or critical slice has no pair of folds.
        return None  # Signal that an S-curve is absent.
    x_fold = float(np.sqrt(1.0 - 1.0 / kappa))  # Compute the positive fold coordinate.
    h_limit = float(kappa * x_fold - np.arctanh(x_fold))  # Compute the positive bias magnitude.
    return {"x_fold": x_fold, "h_limit": h_limit}  # Return the analytic fold data.


def cobweb_path(
    a0: float, kappa: float, delta: float, h: float = 0.0, iterations: int = 18
) -> tuple[np.ndarray, np.ndarray]:
    """Build alternating vertical and horizontal segments for a physical cobweb path."""
    x_coordinates = [a0]  # Start the path on the horizontal axis at the chosen initial state.
    y_coordinates = [0.0]  # Set the initial vertical coordinate to zero.
    current = float(a0)  # Initialize the current physical iterate.
    for _ in range(iterations):  # Apply the map the requested number of times.
        following = float(kapitza_map(current, kappa, delta, h))  # Compute a_{n+1}.
        x_coordinates.extend([current, following])  # Add vertical then horizontal x values.
        y_coordinates.extend([following, following])  # Add the matching y values.
        current = following  # Advance the iteration state.
    return np.asarray(x_coordinates), np.asarray(y_coordinates)  # Return plot-ready arrays.


def save_figure(fig: plt.Figure, output_directory: Path, stem: str) -> None:
    """Save one figure as both a high-resolution PNG and a vector PDF."""
    fig.savefig(output_directory / f"{stem}.png", dpi=240)  # Save the separated raster figure.
    fig.savefig(output_directory / f"{stem}.pdf")  # Save the matching vector PDF figure.
    plt.close(fig)  # Release the figure and its memory after both exports.


def style_axis(axis: plt.Axes) -> None:
    """Apply the shared light-grid plotting style to one axis."""
    axis.grid(True, color="#777777", alpha=0.18, linewidth=0.7)  # Add a subtle reading grid.
    axis.tick_params(direction="in", top=True, right=True)  # Use compact inward ticks.


def plot_fixed_point_regimes(output_directory: Path, delta: float) -> None:
    """Compare fixed-point intersections below, at, and above kappa_c=1."""
    kappas = [0.75, 1.0, 1.5]  # Select representative subcritical, critical, and supercritical cases.
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.1), constrained_layout=True)  # Build aligned panels.
    for axis, kappa in zip(axes, kappas):  # Draw one source-strength regime per panel.
        scale = gamma(kappa, delta)  # Compute this panel's moving physical interval.
        a_grid = np.linspace(0.0, scale, 1400)  # Sample the complete admissible interval.
        roots = physical_states(kappa, delta)  # Compute every physical fixed point.
        axis.plot(a_grid, kapitza_map(a_grid, kappa, delta), color=BLUE, lw=2.0, label=r"$F_\kappa(a)$")  # Draw the map.
        axis.plot(a_grid, a_grid, color=RED, lw=1.5, ls="--", label=r"$a$")  # Draw the identity.
        axis.plot([state["a"] for state in roots], [state["a"] for state in roots], "o", color=BLACK, ms=5, label="Fixed points")  # Mark roots.
        axis.set(xlabel=r"$a$", ylabel=r"$F_\kappa(a)$", title=rf"$\kappa={kappa:.2f}$: {len(roots)} fixed point(s)")  # Label the panel.
        axis.set_xlim(0.0, scale)  # Display exactly the physical horizontal interval.
        axis.set_ylim(0.0, scale)  # Use the same scale vertically for geometric clarity.
        style_axis(axis)  # Apply the common visual style.
    axes[0].legend(loc="upper left")  # Explain the common encoding only once.
    fig.suptitle(r"Fixed-point geometry as the source strength crosses $\kappa_c=1$")  # State the comparison.
    save_figure(fig, output_directory, "fixed_point_regimes")  # Export both figure formats.


def plot_cobweb_regimes(output_directory: Path, delta: float) -> None:
    """Show convergence, critical slowing, and outer-branch selection."""
    kappas = [0.75, 1.0, 1.5]  # Match the fixed-point-regime figure.
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.1), constrained_layout=True)  # Build aligned panels.
    for axis, kappa in zip(axes, kappas):  # Draw one cobweb experiment per source strength.
        scale = gamma(kappa, delta)  # Compute the moving physical interval.
        a_grid = np.linspace(0.0, scale, 1400)  # Sample the fixed-point map smoothly.
        axis.plot(a_grid, kapitza_map(a_grid, kappa, delta), color=BLUE, lw=2.0, label=r"$F_\kappa(a)$")  # Draw the map.
        axis.plot(a_grid, a_grid, color=BLACK, lw=1.3, ls="--", label=r"$a$")  # Draw the identity.
        if kappa < 1.0:  # Use widely separated starts when the center is attracting.
            starts = [0.14 * scale, 0.86 * scale]  # Reveal convergence from both sides.
        elif kappa == 1.0:  # Avoid starts too far away at the slowly convergent threshold.
            starts = [0.22 * scale, 0.78 * scale]  # Keep the transient visible.
        else:  # Place starts near, but not on, the repelling central state.
            starts = [0.43 * scale, 0.57 * scale]  # Reveal selection of both outer branches.
        for start, color, label in zip(starts, [GREEN, RED], ["lower start", "upper start"]):  # Draw both paths.
            path_x, path_y = cobweb_path(start, kappa, delta, iterations=22)  # Generate the path.
            axis.plot(path_x, path_y, color=color, lw=0.9, label=label)  # Display the transient.
        roots = physical_states(kappa, delta)  # Compute the fixed points for reference markers.
        axis.plot([state["a"] for state in roots], [state["a"] for state in roots], "o", color=BLACK, ms=4.5, label="Fixed points")  # Mark all roots.
        axis.set(xlabel=r"$a_n$", ylabel=r"$a_{n+1}$", title=rf"$\kappa={kappa:.2f}$")  # Label the iteration plane.
        axis.set_xlim(0.0, scale)  # Use the physical horizontal interval.
        axis.set_ylim(0.0, scale)  # Match the vertical interval.
        style_axis(axis)  # Apply the common visual style.
    axes[0].legend(loc="upper left", fontsize=8)  # Explain the shared encoding once.
    fig.suptitle("Cobweb iterations reveal stability and branch selection")  # State the main message.
    save_figure(fig, output_directory, "cobweb_regimes")  # Export both formats.


def plot_symmetric_bifurcation(output_directory: Path, delta: float) -> None:
    """Plot the perfect pitchfork in both normalized and physical coordinates."""
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.5), constrained_layout=True)  # Compare coordinates.
    kappa_grid = np.linspace(0.2, 4.5, 1000)  # Sample the central branch over the display range.
    stable_center = kappa_grid <= 1.0  # Select the attracting central segment.
    unstable_center = kappa_grid >= 1.0  # Select the repelling central segment.
    for axis, physical in zip(axes, [False, True]):  # Draw normalized then physical geometry.
        center = np.zeros_like(kappa_grid) if not physical else gamma(kappa_grid, delta) / 2.0  # Build center.
        axis.plot(kappa_grid[stable_center], center[stable_center], color=BLUE, lw=2.1, label="Attracting")  # Solid stable center.
        axis.plot(kappa_grid[unstable_center], center[unstable_center], color=RED, lw=2.1, ls="--", label="Repelling")  # Dashed unstable center.
        x_parameter = np.linspace(1.0e-5, 0.995, 5000)  # Parameterize both outer branches.
        kappa_outer = np.arctanh(x_parameter) / x_parameter  # Evaluate the exact branch relation.
        visible = kappa_outer <= kappa_grid[-1]  # Clip the branches to the displayed range.
        upper = x_parameter if not physical else kappa_outer * (1.0 + x_parameter) / (2.0 * (1.0 + delta))  # Convert upper branch.
        lower = -x_parameter if not physical else kappa_outer * (1.0 - x_parameter) / (2.0 * (1.0 + delta))  # Convert lower branch.
        axis.plot(kappa_outer[visible], upper[visible], color=BLUE, lw=2.1)  # Draw the attracting upper branch.
        axis.plot(kappa_outer[visible], lower[visible], color=BLUE, lw=2.1)  # Draw the attracting lower branch.
        axis.plot(1.0, 0.0 if not physical else 1.0 / (2.0 * (1.0 + delta)), "o", color=BLACK, ms=5, label="Critical point")  # Mark threshold.
        axis.axvline(1.0, color=GRAY, lw=1.0, ls=":")  # Add a vertical critical reference.
        axis.set(xlabel=r"Source strength $\kappa$", ylabel=r"Normalized state $x$" if not physical else r"Interface trace $a$")  # Label coordinates.
        axis.set_title("Standard pitchfork" if not physical else "Physical moving interval")  # Explain the coordinate effect.
        style_axis(axis)  # Apply the common visual style.
    axes[0].legend(loc="best")  # Define the line-style convention once.
    fig.suptitle(r"Supercritical pitchfork generated by $x=\tanh(\kappa x)$")  # State the exact equation.
    save_figure(fig, output_directory, "symmetric_bifurcation")  # Export both formats.


def plot_bias_family(output_directory: Path, delta: float, signs: int) -> None:
    """Plot six positive or negative fixed-bias unfoldings with kappa as control."""
    magnitudes = [0.0, 0.002, 0.005, 0.010, 0.020, 0.040]  # Use six progressively larger small biases.
    biases = [signs * value for value in magnitudes]  # Reflect the family when negative bias is requested.
    fig, axes = plt.subplots(2, 3, figsize=(12.2, 7.2), constrained_layout=True)  # Build six aligned panels.
    for axis, h in zip(axes.flat, biases):  # Continue every selected fixed-bias slice.
        if h == 0.0:  # Draw the exact perfect pitchfork as the reference panel.
            kappa_grid = np.linspace(0.25, 5.0, 1000)  # Sample the central branch.
            center_a = gamma(kappa_grid, delta) / 2.0  # Convert the center to physical a.
            axis.plot(kappa_grid[kappa_grid <= 1.0], center_a[kappa_grid <= 1.0], color=BLUE, lw=1.8, label="Attracting")  # Stable center.
            axis.plot(kappa_grid[kappa_grid >= 1.0], center_a[kappa_grid >= 1.0], color=RED, lw=1.8, ls="--", label="Repelling")  # Unstable center.
            x = np.linspace(1.0e-5, 0.9999999, 8000)  # Reach large kappa values on the outer branches.
            kappa_curve = np.arctanh(x) / x  # Compute the exact source parameter.
            visible = kappa_curve <= 5.0  # Clip to the panel range.
            axis.plot(kappa_curve[visible], kappa_curve[visible] * (1.0 + x[visible]) / (2.0 * (1.0 + delta)), color=BLUE, lw=1.8)  # Upper branch.
            axis.plot(kappa_curve[visible], kappa_curve[visible] * (1.0 - x[visible]) / (2.0 * (1.0 + delta)), color=BLUE, lw=1.8)  # Lower branch.
            axis.plot(1.0, 1.0 / (2.0 * (1.0 + delta)), "o", color=BLACK, ms=4.5, label="Critical point")  # Mark pitchfork.
            axis.set_title(r"$h=0$: perfect pitchfork")  # Identify the symmetric limit.
        else:  # Draw the imperfect continuation parameterized by the state x.
            x = np.concatenate([np.linspace(-0.9999999, -1.0e-4, 16000), np.linspace(1.0e-4, 0.9999999, 16000)])  # Reach kappa=5 while avoiding x=0.
            kappa_curve = (np.arctanh(x) - h) / x  # Recover kappa from the fixed-point equation.
            a_curve = kappa_curve * (1.0 + x) / (2.0 * (1.0 + delta))  # Convert to physical a.
            multiplier = kappa_curve * (1.0 - x * x)  # Evaluate Picard stability along the curve.
            visible = (kappa_curve >= 0.25) & (kappa_curve <= 5.0) & (a_curve >= 0.0)  # Keep physical display points.
            stable = visible & (np.abs(multiplier) < 1.0)  # Select attracting pieces.
            unstable = visible & ~stable  # Select repelling pieces.
            plot_masked_segments(axis, kappa_curve, a_curve, stable, BLUE, "-", "Attracting")  # Draw stable pieces.
            plot_masked_segments(axis, kappa_curve, a_curve, unstable, RED, "--", "Repelling")  # Draw unstable pieces.
            x_fold, kappa_fold = fold_for_bias(h)  # Compute the unique finite saddle node.
            a_fold = kappa_fold * (1.0 + x_fold) / (2.0 * (1.0 + delta))  # Convert its state.
            axis.plot(kappa_fold, a_fold, "o", color=BLACK, ms=4.5, label="Saddle node")  # Mark the fold.
            axis.set_title(rf"$h={h:+.3f}$, $\kappa_{{\rm SN}}={kappa_fold:.3f}$")  # Report its shift.
        axis.set(xlim=(0.25, 5.0), ylim=(0.0, 3.35), xlabel=r"Source strength $\kappa$", ylabel=r"$a$")  # Use shared axes.
        style_axis(axis)  # Apply the common visual style.
    axes.flat[0].legend(loc="center right", fontsize=8)  # Define solid and dashed encodings once.
    direction = "positive" if signs > 0 else "negative"  # Build a readable family description.
    fig.suptitle(rf"Progressive unfolding at six small {direction} biases ($\delta={delta:g}$; $\kappa$ is the control parameter)")  # State the experiment.
    save_figure(fig, output_directory, f"small_{direction}_bias_bifurcations")  # Export both formats.


def plot_masked_segments(axis, x, y, mask, color, linestyle, label) -> None:
    """Plot contiguous masked curve pieces without connecting across gaps."""
    indices = np.flatnonzero(mask)  # Locate all samples retained by the mask.
    if indices.size == 0:  # Skip an empty branch selection safely.
        return  # Leave the axis unchanged when no points are visible.
    breaks = np.where(np.diff(indices) > 1)[0] + 1  # Find gaps between disjoint curve pieces.
    groups = np.split(indices, breaks)  # Split retained indices into contiguous groups.
    first = True  # Add the legend label only to the first visible piece.
    for group in groups:  # Draw each continuous branch separately.
        if group.size > 1:  # Require at least two samples for a line segment.
            axis.plot(x[group], y[group], color=color, ls=linestyle, lw=1.8, label=label if first else None)  # Draw one piece.
            first = False  # Suppress duplicate legend entries on later pieces.


def plot_critical_trajectory(output_directory: Path, delta: float) -> None:
    """Trace both saddle-node trajectories and verify the local two-thirds law."""
    h_positive = np.geomspace(1.0e-5, 0.25, 260)  # Resolve the singular small-bias regime logarithmically.
    h_values = np.concatenate([-h_positive[::-1], [0.0], h_positive])  # Include both reflected signs.
    x_values = []  # Store the saddle-node state for every bias.
    kappa_values = []  # Store the corresponding source threshold.
    a_values = []  # Store the physical fold trace.
    for h in h_values:  # Evaluate the analytic fold equations point by point.
        x_fold, kappa_fold = fold_for_bias(float(h))  # Compute one saddle node.
        x_values.append(x_fold)  # Retain its normalized state.
        kappa_values.append(kappa_fold)  # Retain its source parameter.
        a_values.append(kappa_fold * (1.0 + x_fold) / (2.0 * (1.0 + delta)))  # Convert to a.
    x_values = np.asarray(x_values)  # Convert lists to vectorized arrays.
    kappa_values = np.asarray(kappa_values)  # Convert threshold values to an array.
    a_values = np.asarray(a_values)  # Convert physical folds to an array.
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.0), constrained_layout=True)  # Show trajectory and scaling.
    axes[0].plot(h_values, kappa_values, color=BLUE, lw=2.0)  # Draw the cusp-like threshold trajectory.
    axes[0].plot(0.0, 1.0, "o", color=BLACK, ms=5)  # Mark the perfect pitchfork limit.
    axes[0].set(xlabel=r"Bias $h$", ylabel=r"$\kappa_{\rm SN}$", title="Critical source strength")  # Label the first view.
    axes[1].plot(h_values, a_values, color=BLUE, lw=2.0)  # Draw the physical fold-state trajectory.
    axes[1].plot(0.0, 1.0 / (2.0 * (1.0 + delta)), "o", color=BLACK, ms=5)  # Mark its limit.
    axes[1].set(xlabel=r"Bias $h$", ylabel=r"$a_{\rm SN}$", title="Moving critical state")  # Label the second view.
    positive_mask = h_values > 0.0  # Select one symmetric half for a log-log comparison.
    absolute_h = h_values[positive_mask]  # Extract positive magnitudes.
    shift = kappa_values[positive_mask] - 1.0  # Compute the exact threshold displacement.
    asymptotic = (1.5 * absolute_h) ** (2.0 / 3.0)  # Evaluate the local two-thirds law.
    axes[2].loglog(absolute_h, shift, color=BLUE, lw=2.0, label="Exact fold")  # Draw exact displacement.
    axes[2].loglog(absolute_h, asymptotic, color=RED, lw=1.6, ls="--", label=r"$(3|h|/2)^{2/3}$")  # Draw asymptotic law.
    axes[2].set(xlabel=r"$|h|$", ylabel=r"$\kappa_{\rm SN}-1$", title="Local scaling law")  # Label log axes.
    axes[2].legend(loc="best")  # Identify exact and asymptotic curves.
    for axis in axes:  # Apply the shared style to every panel.
        style_axis(axis)  # Add light grids and inward ticks.
    fig.suptitle("Saddle-node trajectory created by imperfect pitchfork unfolding")  # State the purpose.
    save_figure(fig, output_directory, "critical_point_trajectory")  # Export both formats.


def plot_s_curve(output_directory: Path, delta: float, kappa: float) -> None:
    """Plot the fixed-kappa S-curve and the corresponding heat-flux intersections."""
    fold = fold_data_at_kappa(kappa)  # Obtain the analytic folds for this source slice.
    if fold is None:  # Protect the function from an invalid subcritical request.
        raise ValueError("The S-curve requires kappa>1.")  # Explain the mathematical requirement.
    x = np.linspace(-0.995, 0.995, 6000)  # Parameterize the complete physical branch.
    h_curve = np.arctanh(x) - kappa * x  # Recover the bias required at each state.
    a_curve = kappa * (1.0 + x) / (2.0 * (1.0 + delta))  # Convert to physical interface trace.
    multiplier = kappa * (1.0 - x * x)  # Evaluate Picard stability along the S-curve.
    stable = np.abs(multiplier) < 1.0  # Select the two attracting outer portions.
    unstable = ~stable  # Select the repelling middle portion.
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.4), constrained_layout=True)  # Pair continuation and flux views.
    plot_masked_segments(axes[0], h_curve, a_curve, stable, BLUE, "-", "Attracting branches")  # Draw stable branches.
    plot_masked_segments(axes[0], h_curve, a_curve, unstable, RED, "--", "Repelling middle branch")  # Draw unstable branch.
    x_fold = fold["x_fold"]  # Read the positive fold state.
    h_limit = fold["h_limit"]  # Read the positive fold-bias magnitude.
    fold_xs = np.array([-x_fold, x_fold])  # Assemble both reflected fold states.
    fold_hs = np.array([h_limit, -h_limit])  # Match each state to its bias sign.
    fold_as = kappa * (1.0 + fold_xs) / (2.0 * (1.0 + delta))  # Convert both folds to a.
    axes[0].plot(fold_hs, fold_as, "o", color=BLACK, ms=5, label="Saddle nodes")  # Mark both folds.
    axes[0].set(xlabel=r"Bias $h$", ylabel=r"Interface trace $a$", title=rf"S-curve at fixed $\kappa={kappa:.2f}$")  # Label continuation.
    axes[0].legend(loc="best")  # Explain stability and folds.
    selected_h = 0.5 * h_limit  # Choose a three-root slice halfway to the fold.
    s_grid = np.linspace(0.0, kappa, 2500)  # Sample the complete jump interval.
    q_interface = beta_hat(s_grid, kappa, delta, selected_h) * s_grid  # Compute Kapitza flux.
    q_bulk = delta / (1.0 + delta) * (kappa - s_grid)  # Compute bulk-required flux.
    states = physical_states(kappa, delta, selected_h)  # Compute all intersections on the slice.
    axes[1].plot(s_grid, q_interface, color=RED, lw=2.0, label=r"$\widehat q_K=\widehat\beta(s)s$")  # Draw interface flux.
    axes[1].plot(s_grid, q_bulk, color=BLUE, lw=2.0, label=r"$\widehat q_{\rm bulk}$")  # Draw bulk flux.
    axes[1].plot([state["s"] for state in states], [state["q_hat"] for state in states], "o", color=BLACK, ms=5, label="Steady states")  # Mark roots.
    axes[1].set(xlabel=r"Jump $s=[u]$", ylabel=r"Normalized heat flux $\widehat q$", title=rf"Three intersections at $h={selected_h:.4f}$")  # Label flux view.
    axes[1].legend(loc="best")  # Explain both flux laws.
    for axis in axes:  # Apply the shared style to both panels.
        style_axis(axis)  # Add light grids and inward ticks.
    fig.suptitle(rf"Imperfect bifurcation and heat-flux geometry ($\delta={delta:g}$)")  # State the connection.
    save_figure(fig, output_directory, "s_curve")  # Export both formats.


def plot_biased_cobweb(output_directory: Path, delta: float, kappa: float) -> None:
    """Show three fixed points and two clearly visible attracting basins."""
    fold = fold_data_at_kappa(kappa)  # Obtain the admissible three-root bias range.
    if fold is None:  # Reject an invalid source slice.
        raise ValueError("Biased cobwebs require kappa>1.")  # Explain the requirement.
    h = 0.5 * fold["h_limit"]  # Select an interior slice with three simple roots.
    states = physical_states(kappa, delta, h)  # Compute all roots on this slice.
    if len(states) != 3:  # Guard against accidental parameter or solver regressions.
        raise RuntimeError("Expected exactly three fixed points for the biased cobweb figure.")  # Stop clearly.
    scale = gamma(kappa, delta)  # Compute the physical plotting interval.
    a_grid = np.linspace(0.0, scale, 1800)  # Sample the map smoothly.
    middle = states[1]["a"]  # Identify the repelling basin boundary.
    lower_start = middle - 0.08 * scale  # Start visibly below the repelling root.
    upper_start = middle + 0.08 * scale  # Start visibly above the repelling root.
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.5), constrained_layout=True)  # Pair intersections and cobwebs.
    for axis in axes:  # Draw common map geometry in both panels.
        axis.plot(a_grid, kapitza_map(a_grid, kappa, delta, h), color=BLUE, lw=2.0, label=r"$F_{\kappa,h}(a)$")  # Draw the map.
        axis.plot(a_grid, a_grid, color=BLACK, lw=1.3, ls="--", label=r"$a$")  # Draw the identity.
        axis.plot([state["a"] for state in states], [state["a"] for state in states], "o", color=BLACK, ms=5, label="Fixed points")  # Mark roots.
        axis.set_xlim(0.0, scale)  # Use the physical horizontal interval.
        axis.set_ylim(0.0, scale)  # Match the vertical scale.
        style_axis(axis)  # Apply the common visual style.
    axes[0].set(xlabel=r"$a$", ylabel=r"$F_{\kappa,h}(a)$", title="Three fixed-point intersections")  # Label the geometry panel.
    axes[0].legend(loc="upper left")  # Explain the fixed-point construction.
    for start, color, label in [(lower_start, GREEN, "lower start"), (upper_start, RED, "upper start")]:  # Draw both basins.
        path_x, path_y = cobweb_path(start, kappa, delta, h, iterations=16)  # Generate a visible transient.
        axes[1].plot(path_x, path_y, color=color, lw=0.95, label=label)  # Draw the path.
    axes[1].set(xlabel=r"$a_n$", ylabel=r"$a_{n+1}$", title="Two attracting basins")  # Label the iteration panel.
    axes[1].legend(loc="upper left")  # Identify both starting sides.
    fig.suptitle(rf"Fixed-$h$ slice: $\delta={delta:g}$, $\kappa={kappa:.2f}$, $h={h:.4f}$")  # State all fixed parameters.
    save_figure(fig, output_directory, "biased_fixed_point_cobweb")  # Export both formats.


def plot_temperature_profiles(output_directory: Path, delta: float, kappa: float) -> None:
    """Reconstruct the paper's piecewise temperature fields on all symmetric branches."""
    states = physical_states(kappa, delta, h=0.0)  # Compute all symmetric states at the chosen source.
    z_left = np.linspace(-1.0, 0.0, 500)  # Sample the left material interval.
    z_right = np.linspace(0.0, 1.0, 500)  # Sample the right material interval.
    fig, axis = plt.subplots(figsize=(7.4, 4.6), constrained_layout=True)  # Create one profile comparison.
    colors = [GREEN, RED, BLUE] if len(states) == 3 else [BLUE]  # Distinguish the three branches.
    for index, (state, color) in enumerate(zip(states, colors), start=1):  # Reconstruct each solution.
        a = state["a"]  # Read the left affine coefficient.
        b = state["b"]  # Read the right affine coefficient.
        u_left = a * (z_left + 1.0)  # Evaluate u_I(z)=a(z+1).
        u_right = b * (z_right - 1.0) + kappa * (1.0 - z_right * z_right)  # Evaluate u_E(z).
        stability = "attracting" if state["attracting"] else "repelling"  # Describe Picard stability.
        linestyle = "-" if state["attracting"] else "--"  # Encode stable solutions as solid.
        axis.plot(z_left, u_left, color=color, ls=linestyle, lw=2.0, label=rf"Branch {index}: {stability}")  # Draw left profile.
        axis.plot(z_right, u_right, color=color, ls=linestyle, lw=2.0)  # Continue the same branch on the right.
        axis.plot([0.0, 0.0], [state["a"], kappa - state["b"]], color=color, ls=":", lw=1.1)  # Display the interface jump.
    axis.axvline(0.0, color=BLACK, lw=0.9)  # Mark the material interface.
    axis.axhline(0.0, color=GRAY, lw=0.8)  # Mark the homogeneous boundary level.
    axis.set(xlabel=r"Position $z$", ylabel=r"Temperature $u(z)$", title=rf"Piecewise solutions at $\kappa={kappa:.2f}$ and $h=0$")  # Label the profile plot.
    axis.text(-0.52, axis.get_ylim()[1] * 0.94, r"$\Omega_I$", ha="center")  # Identify the left material.
    axis.text(0.52, axis.get_ylim()[1] * 0.94, r"$\Omega_E$", ha="center")  # Identify the right material.
    axis.legend(loc="best")  # Identify every solution branch.
    style_axis(axis)  # Apply the common visual style.
    save_figure(fig, output_directory, "temperature_profiles")  # Export both formats.


def plot_model_summary(output_directory: Path, delta: float, kappa: float) -> None:
    """Create a compact four-panel summary of the PDE and bifurcation model."""
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.0), constrained_layout=True)  # Build the summary grid.
    scale = gamma(kappa, delta)  # Compute the current physical interval.
    a_grid = np.linspace(0.0, scale, 1600)  # Sample the physical fixed-point map.
    axes[0, 0].plot(a_grid, kapitza_map(a_grid, kappa, delta), color=BLUE, lw=2.0, label=r"$F_\kappa(a)$")  # Draw map.
    axes[0, 0].plot(a_grid, a_grid, color=RED, lw=1.4, ls="--", label=r"$a$")  # Draw identity.
    states = physical_states(kappa, delta)  # Compute all symmetric states.
    axes[0, 0].plot([state["a"] for state in states], [state["a"] for state in states], "o", color=BLACK, ms=5)  # Mark roots.
    axes[0, 0].set(xlabel=r"$a$", ylabel=r"$F_\kappa(a)$", title="Fixed-point intersections")  # Label panel.
    axes[0, 0].legend(loc="best")  # Explain curves.
    s_grid = np.linspace(0.0, kappa, 1600)  # Sample every physical jump.
    axes[0, 1].semilogy(s_grid, beta_hat(s_grid, kappa, delta), color=BLUE, lw=2.0)  # Draw exponential conductance.
    axes[0, 1].set(xlabel=r"Jump $s$", ylabel=r"$\widehat\beta_{\kappa,0}(s)$", title="Normalized Kapitza conductance")  # Label panel.
    axes[1, 0].plot(s_grid, beta_hat(s_grid, kappa, delta) * s_grid, color=RED, lw=2.0, label=r"$\widehat q_K$")  # Draw interface flux.
    axes[1, 0].plot(s_grid, delta * (kappa - s_grid) / (1.0 + delta), color=BLUE, lw=2.0, label=r"$\widehat q_{\rm bulk}$")  # Draw bulk flux.
    axes[1, 0].plot([state["s"] for state in states], [state["q_hat"] for state in states], "o", color=BLACK, ms=5)  # Mark steady states.
    axes[1, 0].set(xlabel=r"Jump $s$", ylabel=r"Normalized flux $\widehat q$", title="Flux intersections")  # Label panel.
    axes[1, 0].legend(loc="best")  # Explain flux curves.
    z = np.linspace(0.0, 1.0, 500)  # Sample the source-supported right layer.
    axes[1, 1].fill_between(z, 0.0, 2.0 * kappa, color=BLUE, alpha=0.18)  # Display the constant source support.
    axes[1, 1].plot(z, np.full_like(z, 2.0 * kappa), color=BLUE, lw=2.0)  # Draw its amplitude.
    axes[1, 1].set(xlim=(-1.0, 1.0), ylim=(0.0, 2.25 * kappa), xlabel=r"Position $z$", ylabel=r"$f(z)/\mu$", title=r"Source $2\kappa\chi_{\Omega_E}$")  # Label source panel.
    axes[1, 1].axvline(0.0, color=BLACK, lw=0.9)  # Mark the interface.
    for axis in axes.flat:  # Apply the shared style to all panels.
        style_axis(axis)  # Add light grids and inward ticks.
    fig.suptitle(rf"Source-controlled nonlinear Kapitza model ($\delta={delta:g}$, $\kappa={kappa:.2f}$)")  # State parameters.
    save_figure(fig, output_directory, "model_summary")  # Export both formats.


def write_state_table(path: Path, delta: float, kappas: list[float], h: float = 0.0) -> None:
    """Write steady states for selected source parameters to a CSV file."""
    with path.open("w", newline="", encoding="utf-8") as stream:  # Open a portable UTF-8 CSV file.
        writer = csv.writer(stream)  # Create a standards-compliant CSV writer.
        writer.writerow(["kappa", "h", "branch", "x", "a", "b", "s", "beta_hat", "q_hat", "multiplier", "attracting"])  # Write headers.
        for kappa in kappas:  # Process every requested source parameter.
            for branch, state in enumerate(physical_states(kappa, delta, h), start=1):  # Number its roots.
                writer.writerow([kappa, h, branch, state["x"], state["a"], state["b"], state["s"], state["beta_hat"], state["q_hat"], state["multiplier"], state["attracting"]])  # Write one state.


def write_fold_table(path: Path, biases: list[float], delta: float) -> None:
    """Write the imperfect-bifurcation saddle-node trajectory to CSV."""
    with path.open("w", newline="", encoding="utf-8") as stream:  # Open the destination table.
        writer = csv.writer(stream)  # Create the CSV serializer.
        writer.writerow(["h", "x_SN", "kappa_SN", "a_SN", "local_asymptotic_shift"])  # Describe every column.
        for h in biases:  # Evaluate each requested bias independently.
            x_fold, kappa_fold = fold_for_bias(h)  # Compute the exact fold.
            a_fold = kappa_fold * (1.0 + x_fold) / (2.0 * (1.0 + delta))  # Convert to physical a.
            approximation = (1.5 * abs(h)) ** (2.0 / 3.0) if h != 0.0 else 0.0  # Evaluate the local law.
            writer.writerow([h, x_fold, kappa_fold, a_fold, approximation])  # Save one trajectory point.


def main() -> None:
    """Run the complete source-controlled Kapitza workflow."""
    parser = argparse.ArgumentParser(description=__doc__)  # Build the command-line interface.
    parser.add_argument("--delta", type=float, default=0.5, help="Conductivity ratio lambda/mu.")  # Expose delta.
    parser.add_argument("--kappa", type=float, default=1.5, help="Source parameter for detailed figures.")  # Expose kappa.
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[2] / "output" / "python", help="Output directory.")  # Expose output path.
    arguments = parser.parse_args()  # Parse all supplied options.
    if arguments.delta <= 0.0:  # Enforce the positive-conductivity-ratio assumption.
        parser.error("delta must be positive")  # Report an invalid material ratio clearly.
    if arguments.kappa <= 1.0:  # Require the detailed slice to possess three symmetric roots.
        parser.error("kappa must exceed one for the S-curve and three-state figures")  # Report the requirement.
    arguments.output.mkdir(parents=True, exist_ok=True)  # Create the output directory if necessary.
    plot_fixed_point_regimes(arguments.output, arguments.delta)  # Generate intersection regimes.
    plot_cobweb_regimes(arguments.output, arguments.delta)  # Generate matching cobweb regimes.
    plot_symmetric_bifurcation(arguments.output, arguments.delta)  # Generate the perfect pitchfork.
    plot_bias_family(arguments.output, arguments.delta, signs=1)  # Generate positive-bias unfoldings.
    plot_bias_family(arguments.output, arguments.delta, signs=-1)  # Generate negative-bias unfoldings.
    plot_critical_trajectory(arguments.output, arguments.delta)  # Generate the fold trajectory and scaling.
    plot_s_curve(arguments.output, arguments.delta, arguments.kappa)  # Generate the fixed-kappa S-curve.
    plot_biased_cobweb(arguments.output, arguments.delta, arguments.kappa)  # Generate biased cobweb basins.
    plot_temperature_profiles(arguments.output, arguments.delta, arguments.kappa)  # Reconstruct PDE solutions.
    plot_model_summary(arguments.output, arguments.delta, arguments.kappa)  # Generate the compact summary.
    write_state_table(arguments.output / "steady_states.csv", arguments.delta, [0.75, 1.0, arguments.kappa])  # Save benchmark roots.
    table_biases = [-0.04, -0.02, -0.01, -0.005, -0.002, 0.0, 0.002, 0.005, 0.01, 0.02, 0.04]  # Select fold samples.
    write_fold_table(arguments.output / "critical_trajectory.csv", table_biases, arguments.delta)  # Save exact folds.
    print(f"delta=lambda/mu={arguments.delta:g}")  # Report the selected material ratio.
    print(f"symmetric threshold kappa_c=1")  # Report the analytic critical source strength.
    print(f"detailed source parameter kappa={arguments.kappa:g}")  # Report the detailed slice.
    for state in physical_states(arguments.kappa, arguments.delta):  # Print every detailed steady state.
        print(f"x={state['x']:+.10f}, a={state['a']:.10f}, s={state['s']:.10f}, multiplier={state['multiplier']:.10f}, attracting={state['attracting']}")  # Summarize one state.
    print(f"outputs written to {arguments.output}")  # Confirm the output location.


if __name__ == "__main__":  # Run the workflow only when this file is executed directly.
    main()  # Invoke the complete reproducible calculation.
