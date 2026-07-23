%KAPITZA_KAPPA_MODEL Complete source-controlled nonlinear Kapitza study.
% This single MATLAB script reproduces the Python calculations, generates all
% report figures, writes numerical tables, and documents every computational
% step in English.  The bifurcation parameter kappa is the normalized positive
% heat-source strength; delta=lambda/mu is the conductivity ratio.

clearvars;                                      % Remove variables left by previous sessions.
close all;                                      % Close figures left by previous sessions.
clc;                                            % Clear the Command Window for readable output.

set(groot, 'defaultAxesFontName', 'Times New Roman'); % Use Times for ticks and axis text.
set(groot, 'defaultTextFontName', 'Times New Roman'); % Use Times for titles and annotations.
set(groot, 'defaultLegendFontName', 'Times New Roman'); % Use Times for every legend.
set(groot, 'defaultAxesTickLabelInterpreter', 'tex'); % Render Greek symbols with TeX syntax.
set(groot, 'defaultTextInterpreter', 'tex');     % Keep prose in Times and formulas italicized.
set(groot, 'defaultLegendInterpreter', 'tex');   % Apply the same policy inside legends.
set(groot, 'defaultAxesFontSize', 11);           % Set one readable base font size.

delta = 0.5;                                    % Set delta=lambda/mu as in the paper example.
kappaDetail = 1.5;                              % Select a supercritical source-strength slice.
outputDirectory = kapitza_output_directory();   % Locate and create output/matlab robustly.

fprintf('Running the source-controlled nonlinear Kapitza workflow.\n'); % Announce execution.
workflowStage = getenv('KAPITZA_STAGE');         % Read an optional batch-stage selector.
if isempty(workflowStage)                        % Use one complete run during normal execution.
    workflowStage = 'all';                       % Select every output by default.
end                                              % Finish the default-stage selection.
if strcmpi(workflowStage, 'all') || strcmpi(workflowStage, 'early') % Select the first figure group.
    plot_fixed_point_regimes(outputDirectory, delta); % Compare intersections across kappa_c=1.
    plot_cobweb_regimes(outputDirectory, delta); % Compare Picard iteration across regimes.
    plot_symmetric_bifurcation(outputDirectory, delta); % Draw normalized and physical pitchforks.
    plot_bias_family(outputDirectory, delta, 1); % Draw six small positive-bias unfoldings.
    plot_bias_family(outputDirectory, delta, -1); % Draw six reflected negative-bias unfoldings.
    plot_critical_trajectory(outputDirectory, delta); % Trace folds and verify the two-thirds law.
end                                              % Finish the first batch group.
if strcmpi(workflowStage, 'all') || strcmpi(workflowStage, 'late') % Select the second figure group.
    plot_s_curve(outputDirectory, delta, kappaDetail); % Draw fixed-kappa S-curve and flux roots.
    plot_biased_cobweb(outputDirectory, delta, kappaDetail); % Show two attracting basins.
    plot_temperature_profiles(outputDirectory, delta, kappaDetail); % Reconstruct piecewise PDE states.
    plot_model_summary(outputDirectory, delta, kappaDetail); % Create a compact four-panel overview.
    write_state_table(outputDirectory, delta, [0.75, 1.0, kappaDetail]); % Save benchmark states.
    write_fold_table(outputDirectory, delta);    % Save representative saddle-node data.
end                                              % Finish the second batch group.

states = physical_states(kappaDetail, delta, 0); % Recompute detailed states for console output.
fprintf('delta=lambda/mu=%.12g\n', delta);        % Report the conductivity ratio.
fprintf('symmetric threshold kappa_c=1\n');       % Report the analytic source threshold.
fprintf('detailed source parameter kappa=%.12g\n', kappaDetail); % Report the selected slice.
for index = 1:numel(states)                      % Print every steady state on the selected slice.
    state = states(index);                       % Read one state structure.
    fprintf('x=%+.10f, a=%.10f, s=%.10f, multiplier=%.10f, attracting=%d\n', ...
        state.x, state.a, state.s, state.multiplier, state.attracting); % Summarize it.
end                                             % Finish the state-reporting loop.
fprintf('Outputs written to %s\n', outputDirectory); % Confirm the destination directory.


function value = gamma_value(kappa, delta)      % Evaluate the moving physical upper bound.
%GAMMA_VALUE Return gamma(kappa)=kappa/(1+delta).
value = kappa ./ (1 + delta);                   % Apply the definition element by element.
end                                             % Finish the gamma helper.


function value = beta_hat(s, kappa, delta, h)   % Evaluate normalized interface conductance.
%BETA_HAT Return beta_phys/mu for the exponential Kapitza law.
if nargin < 4                                  % Check whether the optional bias was omitted.
    h = 0;                                     % Use the symmetric interface law by default.
end                                             % Finish the optional-argument check.
exponent = 2 .* (kappa - 2 .* s) + 2 .* h;     % Assemble the biased exponential argument.
value = delta ./ (1 + delta) .* exp(exponent); % Apply the normalized material prefactor.
end                                             % Finish the conductance helper.


function value = kapitza_map(a, kappa, delta, h) % Evaluate the physical fixed-point map.
%KAPITZA_MAP Compute F_{kappa,h}(a) from the interface balance.
if nargin < 4                                  % Check whether the optional bias was omitted.
    h = 0;                                     % Use the perfect symmetric law by default.
end                                             % Finish the optional-argument check.
s = kappa - (1 + delta) .* a;                  % Reconstruct the interface jump from a.
conductance = beta_hat(s, kappa, delta, h);    % Evaluate the nonlinear interface law.
denominator = delta + (1 + delta) .* conductance; % Assemble the scalar balance denominator.
value = kappa .* conductance ./ denominator;   % Solve the interface balance for the next a.
end                                             % Finish the fixed-point map.


function roots = normalized_roots(kappa, h)    % Find all roots of x=tanh(kappa*x+h).
%NORMALIZED_ROOTS Return every physical normalized fixed point in (-1,1).
if nargin < 2                                  % Check whether the optional bias was omitted.
    h = 0;                                     % Solve the symmetric problem by default.
end                                             % Finish the optional-argument check.
if abs(h) < 1e-14                              % Exploit exact odd symmetry when h is zero.
    if kappa <= 1                              % Detect the unique-root regime and threshold.
        roots = 0;                             % Return the exact central root only.
    else                                       % Handle the three-root supercritical regime.
        residual = @(x) tanh(kappa .* x) - x;  % Define the positive-root residual.
        positive = fzero(residual, [1e-10, 1 - 1e-10]); % Isolate the positive outer root.
        roots = [-positive, 0, positive];       % Restore the exact symmetric triplet.
    end                                         % Finish the symmetric-regime split.
    return                                      % Skip the generic scanning algorithm.
end                                             % Finish the exact-symmetry branch.
grid = linspace(-1 + 1e-10, 1 - 1e-10, 20001); % Sample the complete open interval.
values = tanh(kappa .* grid + h) - grid;       % Evaluate the residual vectorially.
candidates = zeros(1, 3);                      % Preallocate the maximum number of roots.
count = 0;                                     % Initialize the detected-root counter.
residual = @(x) tanh(kappa .* x + h) - x;      % Create a scalar residual for fzero.
for index = 1:(numel(grid) - 1)                % Inspect every neighboring scan pair.
    if values(index) * values(index + 1) < 0   % Detect a simple root by a sign change.
        count = count + 1;                     % Advance the root counter.
        candidates(count) = fzero(residual, [grid(index), grid(index + 1)]); % Refine it.
    end                                         % Finish the sign-change test.
end                                             % Finish the interval scan.
roots = unique(round(candidates(1:count), 12)); % Remove duplicate numerical detections.
end                                             % Finish the root solver.


function states = physical_states(kappa, delta, h) % Reconstruct all physical state variables.
%PHYSICAL_STATES Convert normalized roots into a, b, s, flux, and stability.
if nargin < 3                                  % Check whether the optional bias was omitted.
    h = 0;                                     % Use the symmetric model by default.
end                                             % Finish the optional-argument check.
xRoots = normalized_roots(kappa, h);           % Compute every normalized fixed point.
template = struct('x', 0, 'a', 0, 'b', 0, 's', 0, ... % Define a consistent record layout.
    'betaHat', 0, 'qHat', 0, 'multiplier', 0, 'attracting', false); % Complete the layout.
states = repmat(template, 1, numel(xRoots));    % Preallocate one record per root.
scale = gamma_value(kappa, delta);              % Compute the moving physical interval.
for index = 1:numel(xRoots)                    % Reconstruct every branch independently.
    x = xRoots(index);                         % Read one normalized state.
    a = 0.5 * scale * (1 + x);                 % Recover the left interface trace.
    s = 0.5 * kappa * (1 - x);                 % Recover the temperature jump.
    b = delta * a;                             % Recover the right affine coefficient.
    conductance = beta_hat(s, kappa, delta, h); % Evaluate normalized interface conductance.
    multiplier = kappa * (1 - x ^ 2);          % Differentiate the normalized Picard map.
    states(index).x = x;                       % Store the normalized fixed point.
    states(index).a = a;                       % Store the left interface trace.
    states(index).b = b;                       % Store the right affine coefficient.
    states(index).s = s;                       % Store the jump.
    states(index).betaHat = conductance;       % Store normalized conductance.
    states(index).qHat = delta * a;            % Store normalized common heat flux.
    states(index).multiplier = multiplier;     % Store the Picard multiplier.
    states(index).attracting = abs(multiplier) < 1; % Apply local iteration stability.
end                                             % Finish the reconstruction loop.
end                                             % Finish the physical-state helper.


function [xFold, kappaFold] = fold_for_bias(h) % Evaluate one imperfect saddle node.
%FOLD_FOR_BIAS Solve the exact fold equations for a fixed nonzero bias.
if abs(h) < 1e-15                              % Recover the perfect-pitchfork limit exactly.
    xFold = 0;                                 % Set the critical normalized state.
    kappaFold = 1;                             % Set the critical source strength.
    return                                      % Skip the nonlinear root search.
end                                             % Finish the zero-bias case.
residual = @(x) atanh(x) - x ./ (1 - x .^ 2) - h; % Combine stationarity and fixed point.
if h > 0                                      % Locate a positive-bias fold on negative x.
    xFold = fzero(residual, [-0.999999, -1e-12]); % Isolate the unique fold state.
else                                           % Locate the reflected negative-bias fold.
    xFold = fzero(residual, [1e-12, 0.999999]); % Isolate the unique positive state.
end                                             % Finish the sign-dependent bracketing.
kappaFold = 1 / (1 - xFold ^ 2);              % Apply the exact tangency condition.
end                                             % Finish the fold helper.


function [pathX, pathY] = cobweb_path(a0, kappa, delta, h, iterations)
%COBWEB_PATH Build alternating vertical and horizontal Picard segments.
pathX = a0;                                    % Start on the horizontal axis at a0.
pathY = 0;                                     % Set the initial vertical coordinate to zero.
current = a0;                                  % Initialize the current physical iterate.
for index = 1:iterations                       % Apply the fixed-point map repeatedly.
    following = kapitza_map(current, kappa, delta, h); % Compute the next iterate.
    pathX = [pathX, current, following]; %#ok<AGROW> % Append vertical and horizontal x values.
    pathY = [pathY, following, following]; %#ok<AGROW> % Append matching y values.
    current = following;                       % Advance the iteration state.
end                                             % Finish the requested iterations.
end                                             % Finish the cobweb helper.


function plot_fixed_point_regimes(outputDirectory, delta)
%PLOT_FIXED_POINT_REGIMES Compare intersections below, at, and above kappa_c.
kappas = [0.75, 1.0, 1.5];                    % Select three representative source strengths.
fig = figure('Color', 'w', 'Position', [80, 80, 1300, 430]); % Create a wide panel figure.
layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact'); % Align panels.
title(layout, 'Fixed-point geometry as the source strength crosses \kappa_c=1'); % State comparison.
for panel = 1:numel(kappas)                    % Draw one regime per panel.
    kappa = kappas(panel);                     % Read this panel's source strength.
    scale = gamma_value(kappa, delta);         % Compute its physical interval.
    aGrid = linspace(0, scale, 1400);          % Sample the complete interval.
    states = physical_states(kappa, delta, 0); % Compute every fixed point.
    ax = nexttile(layout);                     % Select the next aligned axis.
    plot(ax, aGrid, kapitza_map(aGrid, kappa, delta, 0), 'Color', [0.12, 0.37, 0.75], 'LineWidth', 2, 'DisplayName', 'F_{\kappa}(a)'); % Draw map.
    hold(ax, 'on');                            % Retain the axis for reference geometry.
    plot(ax, aGrid, aGrid, '--', 'Color', [0.77, 0.24, 0.17], 'LineWidth', 1.5, 'DisplayName', 'a'); % Draw identity.
    rootA = [states.a];                        % Collect all physical fixed points.
    plot(ax, rootA, rootA, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'DisplayName', 'Fixed points'); % Mark roots.
    xlabel(ax, 'a');                           % Label the physical state axis.
    ylabel(ax, 'F_{\kappa}(a)');              % Label the map axis.
    title(ax, sprintf('\\kappa=%.2f: %d fixed point(s)', kappa, numel(states))); % Report count.
    xlim(ax, [0, scale]);                      % Use the exact physical horizontal range.
    ylim(ax, [0, scale]);                      % Match the vertical range.
    style_axis(ax);                            % Apply the shared light-grid style.
    if panel == 1                              % Explain shared encodings only once.
        legend(ax, 'Location', 'northwest');   % Place the common legend unobtrusively.
    end                                         % Finish the legend condition.
end                                             % Finish the regime loop.
save_figure(fig, outputDirectory, 'fixed_point_regimes'); % Export PNG and PDF.
end                                             % Finish the intersection figure.


function plot_cobweb_regimes(outputDirectory, delta)
%PLOT_COBWEB_REGIMES Show convergence, critical slowing, and branch selection.
kappas = [0.75, 1.0, 1.5];                    % Match the fixed-point-regime cases.
fig = figure('Color', 'w', 'Position', [80, 80, 1300, 430]); % Create a wide panel figure.
layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact'); % Align panels.
title(layout, 'Cobweb iterations reveal stability and branch selection'); % State the message.
for panel = 1:numel(kappas)                    % Draw one source regime per panel.
    kappa = kappas(panel);                     % Read the current source strength.
    scale = gamma_value(kappa, delta);         % Compute its physical state interval.
    aGrid = linspace(0, scale, 1400);          % Sample the fixed-point map smoothly.
    ax = nexttile(layout);                     % Select the next aligned axis.
    plot(ax, aGrid, kapitza_map(aGrid, kappa, delta, 0), 'Color', [0.12, 0.37, 0.75], 'LineWidth', 2, 'DisplayName', 'F_{\kappa}(a)'); % Draw map.
    hold(ax, 'on');                            % Retain the axis for all paths.
    plot(ax, aGrid, aGrid, 'k--', 'LineWidth', 1.3, 'DisplayName', 'a'); % Draw identity.
    if kappa < 1                               % Choose broad starts in the attracting-center case.
        starts = [0.14, 0.86] .* scale;       % Reveal convergence from both sides.
    elseif kappa == 1                          % Choose moderate starts at criticality.
        starts = [0.22, 0.78] .* scale;       % Keep slow transients visible.
    else                                       % Choose starts around the repelling center.
        starts = [0.43, 0.57] .* scale;       % Reveal selection of both outer branches.
    end                                         % Finish the start selection.
    colors = {[0.15, 0.55, 0.35], [0.77, 0.24, 0.17]}; % Assign lower and upper path colors.
    labels = {'lower start', 'upper start'};   % Name both initial-condition sides.
    for pathIndex = 1:2                        % Draw both cobweb paths.
        [pathX, pathY] = cobweb_path(starts(pathIndex), kappa, delta, 0, 22); % Generate path.
        plot(ax, pathX, pathY, 'Color', colors{pathIndex}, 'LineWidth', 0.9, 'DisplayName', labels{pathIndex}); % Draw path.
    end                                         % Finish the cobweb-path loop.
    states = physical_states(kappa, delta, 0); % Compute fixed points for markers.
    rootA = [states.a];                        % Collect their physical coordinates.
    plot(ax, rootA, rootA, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4.5, 'DisplayName', 'Fixed points'); % Mark roots.
    xlabel(ax, 'a_n');                         % Label the current iterate.
    ylabel(ax, 'a_{n+1}');                     % Label the next iterate.
    title(ax, sprintf('\\kappa=%.2f', kappa)); % Identify the source regime.
    xlim(ax, [0, scale]);                      % Use the physical horizontal interval.
    ylim(ax, [0, scale]);                      % Match the vertical interval.
    style_axis(ax);                            % Apply the shared axis style.
    if panel == 1                              % Explain shared encodings only once.
        legend(ax, 'Location', 'northwest', 'FontSize', 8); % Place the common legend.
    end                                         % Finish the legend condition.
end                                             % Finish the regime loop.
save_figure(fig, outputDirectory, 'cobweb_regimes'); % Export PNG and PDF.
end                                             % Finish the cobweb-regime figure.


function plot_symmetric_bifurcation(outputDirectory, delta)
%PLOT_SYMMETRIC_BIFURCATION Compare standard and physical pitchfork geometry.
fig = figure('Color', 'w', 'Position', [100, 100, 1080, 470]); % Create two-panel figure.
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); % Align panels.
title(layout, 'Supercritical pitchfork generated by x=tanh(\kappa x)'); % State exact equation.
kappaGrid = linspace(0.2, 4.5, 1000);         % Sample the complete central branch.
xParameter = linspace(1e-5, 0.9998, 6000);    % Parameterize the positive outer state.
kappaOuter = atanh(xParameter) ./ xParameter; % Evaluate the exact branch relation.
visible = kappaOuter <= 4.5;                  % Clip the branches to the display range.
for panel = 1:2                                % Draw normalized then physical coordinates.
    ax = nexttile(layout);                     % Select one comparison axis.
    hold(ax, 'on');                            % Retain all branch segments.
    if panel == 1                              % Construct the standard normalized pitchfork.
        center = zeros(size(kappaGrid));       % Keep the central state at x=0.
        upper = xParameter;                    % Use x directly on the upper branch.
        lower = -xParameter;                   % Reflect x for the lower branch.
        criticalState = 0;                     % Set the critical normalized state.
        yLabel = 'Normalized state x';         % Label the normalized vertical axis.
        panelTitle = 'Standard pitchfork';     % Explain this coordinate choice.
    else                                       % Construct the physical moving-interval view.
        center = gamma_value(kappaGrid, delta) ./ 2; % Convert the center to a(kappa).
        upper = kappaOuter .* (1 + xParameter) ./ (2 * (1 + delta)); % Convert upper branch.
        lower = kappaOuter .* (1 - xParameter) ./ (2 * (1 + delta)); % Convert lower branch.
        criticalState = 1 / (2 * (1 + delta)); % Convert the critical point.
        yLabel = 'Interface trace a';          % Label the physical vertical axis.
        panelTitle = 'Physical moving interval'; % Explain the tilted center branch.
    end                                         % Finish the coordinate selection.
    plot(ax, kappaGrid(kappaGrid <= 1), center(kappaGrid <= 1), 'Color', [0.12, 0.37, 0.75], 'LineWidth', 2.1, 'DisplayName', 'Attracting'); % Stable center.
    plot(ax, kappaGrid(kappaGrid >= 1), center(kappaGrid >= 1), '--', 'Color', [0.77, 0.24, 0.17], 'LineWidth', 2.1, 'DisplayName', 'Repelling'); % Unstable center.
    plot(ax, kappaOuter(visible), upper(visible), 'Color', [0.12, 0.37, 0.75], 'LineWidth', 2.1, 'HandleVisibility', 'off'); % Upper branch.
    plot(ax, kappaOuter(visible), lower(visible), 'Color', [0.12, 0.37, 0.75], 'LineWidth', 2.1, 'HandleVisibility', 'off'); % Lower branch.
    plot(ax, 1, criticalState, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'DisplayName', 'Critical point'); % Mark threshold.
    xline(ax, 1, ':', 'Color', [0.42, 0.42, 0.42], 'HandleVisibility', 'off'); % Add reference.
    xlabel(ax, 'Source strength \kappa');      % Label the control parameter.
    ylabel(ax, yLabel);                        % Label the selected state coordinate.
    title(ax, panelTitle);                     % Identify the coordinate system.
    style_axis(ax);                            % Apply the shared style.
    if panel == 1                              % Explain stability encodings once.
        legend(ax, 'Location', 'best');        % Place the common legend.
    end                                         % Finish the legend condition.
end                                             % Finish the coordinate-comparison loop.
save_figure(fig, outputDirectory, 'symmetric_bifurcation'); % Export PNG and PDF.
end                                             % Finish the pitchfork figure.


function plot_bias_family(outputDirectory, delta, signValue)
%PLOT_BIAS_FAMILY Draw six positive or negative fixed-h continuations in kappa.
magnitudes = [0, 0.002, 0.005, 0.010, 0.020, 0.040]; % Select six small magnitudes.
biases = signValue .* magnitudes;               % Reflect the family for negative bias.
fig = figure('Color', 'w', 'Position', [50, 50, 1250, 760]); % Create six-panel figure.
layout = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact'); % Align panels.
if signValue > 0                               % Build a descriptive positive-family title.
    direction = 'positive';                    % Store the sign word for title and filename.
else                                           % Build the reflected negative-family title.
    direction = 'negative';                    % Store the sign word for title and filename.
end                                             % Finish the direction selection.
title(layout, sprintf('Progressive unfolding at six small %s biases (\\delta=%.1f; \\kappa is the control parameter)', direction, delta)); % State experiment.
for panel = 1:numel(biases)                    % Draw one fixed-bias continuation per panel.
    h = biases(panel);                         % Read this panel's bias.
    ax = nexttile(layout);                     % Select the next aligned axis.
    hold(ax, 'on');                            % Retain all branch pieces.
    if h == 0                                  % Draw the perfect pitchfork reference panel.
        kappaGrid = linspace(0.25, 5, 1200);  % Sample the central branch.
        centerA = gamma_value(kappaGrid, delta) ./ 2; % Convert center to physical a.
        plot(ax, kappaGrid(kappaGrid <= 1), centerA(kappaGrid <= 1), 'Color', [0.12, 0.37, 0.75], 'LineWidth', 1.8, 'DisplayName', 'Attracting'); % Stable center.
        plot(ax, kappaGrid(kappaGrid >= 1), centerA(kappaGrid >= 1), '--', 'Color', [0.77, 0.24, 0.17], 'LineWidth', 1.8, 'DisplayName', 'Repelling'); % Unstable center.
        x = linspace(1e-5, 0.9999999, 9000);  % Reach kappa=5 on the outer branches.
        kappaCurve = atanh(x) ./ x;           % Evaluate the exact branch relation.
        visible = kappaCurve <= 5;            % Clip to the common source range.
        plot(ax, kappaCurve(visible), kappaCurve(visible) .* (1 + x(visible)) ./ (2 * (1 + delta)), 'Color', [0.12, 0.37, 0.75], 'LineWidth', 1.8, 'HandleVisibility', 'off'); % Upper branch.
        plot(ax, kappaCurve(visible), kappaCurve(visible) .* (1 - x(visible)) ./ (2 * (1 + delta)), 'Color', [0.12, 0.37, 0.75], 'LineWidth', 1.8, 'HandleVisibility', 'off'); % Lower branch.
        plot(ax, 1, 1 / (2 * (1 + delta)), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4.5, 'DisplayName', 'Critical point'); % Mark pitchfork.
        title(ax, 'h=0: perfect pitchfork');   % Identify the symmetric limit.
    else                                       % Draw an imperfect continuation in parametric form.
        x = [linspace(-0.9999999, -1e-4, 16000), linspace(1e-4, 0.9999999, 16000)]; % Avoid x=0.
        kappaCurve = (atanh(x) - h) ./ x;      % Recover kappa from the fixed-point equation.
        aCurve = kappaCurve .* (1 + x) ./ (2 * (1 + delta)); % Convert to physical a.
        multiplier = kappaCurve .* (1 - x .^ 2); % Evaluate Picard stability.
        visible = kappaCurve >= 0.25 & kappaCurve <= 5 & aCurve >= 0; % Keep physical points.
        stable = visible & abs(multiplier) < 1; % Select attracting branch pieces.
        unstable = visible & ~stable;          % Select repelling branch pieces.
        stableA = aCurve;                      % Copy the state curve for masked plotting.
        stableA(~stable) = NaN;                % Break the line outside stable pieces.
        unstableA = aCurve;                    % Copy the state curve for unstable masking.
        unstableA(~unstable) = NaN;            % Break the line outside unstable pieces.
        plot(ax, kappaCurve, stableA, 'Color', [0.12, 0.37, 0.75], 'LineWidth', 1.8, 'DisplayName', 'Attracting'); % Draw stable pieces.
        plot(ax, kappaCurve, unstableA, '--', 'Color', [0.77, 0.24, 0.17], 'LineWidth', 1.8, 'DisplayName', 'Repelling'); % Draw unstable pieces.
        [xFold, kappaFold] = fold_for_bias(h); % Compute the unique finite saddle node.
        aFold = kappaFold * (1 + xFold) / (2 * (1 + delta)); % Convert its state.
        plot(ax, kappaFold, aFold, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4.5, 'DisplayName', 'Saddle node'); % Mark fold.
        title(ax, sprintf('h=%+.3f, \\kappa_{SN}=%.3f', h, kappaFold)); % Report its shift.
    end                                         % Finish the perfect/imperfect split.
    xlim(ax, [0.25, 5]);                       % Use a common source interval.
    ylim(ax, [0, 3.35]);                       % Use a common physical state interval.
    xlabel(ax, 'Source strength \kappa');      % Label the control parameter.
    ylabel(ax, 'a');                           % Label the physical state.
    style_axis(ax);                            % Apply the common visual style.
    if panel == 1                              % Explain stability encoding only once.
        legend(ax, 'Location', 'east', 'FontSize', 8); % Place the shared legend.
    end                                         % Finish the legend condition.
end                                             % Finish the bias-family loop.
save_figure(fig, outputDirectory, sprintf('small_%s_bias_bifurcations', direction)); % Export both formats.
end                                             % Finish the bias-family figure.


function plot_critical_trajectory(outputDirectory, delta)
%PLOT_CRITICAL_TRAJECTORY Trace folds and compare with the two-thirds law.
hPositive = logspace(-5, log10(0.25), 260);    % Resolve small biases logarithmically.
hValues = [-fliplr(hPositive), 0, hPositive];  % Include both reflected signs and the origin.
xValues = zeros(size(hValues));                % Preallocate fold states.
kappaValues = zeros(size(hValues));            % Preallocate critical source strengths.
aValues = zeros(size(hValues));                % Preallocate physical fold traces.
for index = 1:numel(hValues)                   % Evaluate every exact fold point.
    [xValues(index), kappaValues(index)] = fold_for_bias(hValues(index)); % Solve fold equations.
    aValues(index) = kappaValues(index) * (1 + xValues(index)) / (2 * (1 + delta)); % Convert state.
end                                             % Finish the trajectory evaluation.
fig = figure('Color', 'w', 'Position', [70, 70, 1320, 430]); % Create three-panel figure.
layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact'); % Align views.
title(layout, 'Saddle-node trajectory created by imperfect pitchfork unfolding'); % State purpose.
ax = nexttile(layout);                         % Select the critical-source panel.
plot(ax, hValues, kappaValues, 'Color', [0.12, 0.37, 0.75], 'LineWidth', 2); % Draw cusp trajectory.
hold(ax, 'on');                                % Retain the critical marker.
plot(ax, 0, 1, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5); % Mark perfect limit.
xlabel(ax, 'Bias h');                          % Label the unfolding parameter.
ylabel(ax, '\kappa_{SN}');                    % Label the critical source strength.
title(ax, 'Critical source strength');         % Explain this view.
style_axis(ax);                                % Apply common style.
ax = nexttile(layout);                         % Select the physical-state panel.
plot(ax, hValues, aValues, 'Color', [0.12, 0.37, 0.75], 'LineWidth', 2); % Draw moving state.
hold(ax, 'on');                                % Retain the critical marker.
plot(ax, 0, 1 / (2 * (1 + delta)), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5); % Mark limit.
xlabel(ax, 'Bias h');                          % Label the unfolding parameter.
ylabel(ax, 'a_{SN}');                          % Label the physical fold trace.
title(ax, 'Moving critical state');            % Explain this view.
style_axis(ax);                                % Apply common style.
ax = nexttile(layout);                         % Select the log-log scaling panel.
positiveMask = hValues > 0;                    % Select one symmetric half.
absoluteH = hValues(positiveMask);             % Extract positive bias magnitudes.
shift = kappaValues(positiveMask) - 1;         % Compute exact threshold displacement.
asymptotic = (1.5 .* absoluteH) .^ (2 / 3);    % Evaluate the local two-thirds law.
loglog(ax, absoluteH, shift, 'Color', [0.12, 0.37, 0.75], 'LineWidth', 2, 'DisplayName', 'Exact fold'); % Draw exact shift.
hold(ax, 'on');                                % Retain the asymptotic comparison.
loglog(ax, absoluteH, asymptotic, '--', 'Color', [0.77, 0.24, 0.17], 'LineWidth', 1.6, 'DisplayName', '(3|h|/2)^{2/3}'); % Draw local law.
xlabel(ax, '|h|');                             % Label bias magnitude.
ylabel(ax, '\kappa_{SN}-1');                  % Label threshold displacement.
title(ax, 'Local scaling law');                % Explain the comparison.
legend(ax, 'Location', 'northwest');           % Identify exact and asymptotic curves.
style_axis(ax);                                % Apply common style.
save_figure(fig, outputDirectory, 'critical_point_trajectory'); % Export both formats.
end                                             % Finish the trajectory figure.


function plot_s_curve(outputDirectory, delta, kappa)
%PLOT_S_CURVE Draw the fixed-kappa S-curve and heat-flux intersections.
xFold = sqrt(1 - 1 / kappa);                  % Compute the positive fold state analytically.
hLimit = kappa * xFold - atanh(xFold);         % Compute the positive bias magnitude.
x = linspace(-0.995, 0.995, 6000);             % Parameterize the complete S-curve.
hCurve = atanh(x) - kappa .* x;                % Recover the bias at each fixed point.
aCurve = kappa .* (1 + x) ./ (2 * (1 + delta)); % Convert to physical interface trace.
multiplier = kappa .* (1 - x .^ 2);            % Evaluate Picard stability.
stable = abs(multiplier) < 1;                  % Select both attracting outer branches.
stableA = aCurve;                              % Copy the curve for stable masking.
stableA(~stable) = NaN;                        % Break it across the unstable middle.
unstableA = aCurve;                            % Copy the curve for unstable masking.
unstableA(stable) = NaN;                       % Retain only the repelling middle.
fig = figure('Color', 'w', 'Position', [90, 90, 1120, 470]); % Create paired panels.
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); % Align panels.
title(layout, sprintf('Imperfect bifurcation and heat-flux geometry (\\delta=%.1f)', delta)); % State relation.
ax = nexttile(layout);                         % Select the S-curve panel.
plot(ax, hCurve, stableA, 'Color', [0.12, 0.37, 0.75], 'LineWidth', 2, 'DisplayName', 'Attracting branches'); % Draw stable pieces.
hold(ax, 'on');                                % Retain middle branch and folds.
plot(ax, hCurve, unstableA, '--', 'Color', [0.77, 0.24, 0.17], 'LineWidth', 2, 'DisplayName', 'Repelling middle branch'); % Draw middle.
foldXs = [-xFold, xFold];                      % Assemble both reflected fold states.
foldHs = [hLimit, -hLimit];                    % Match each state to its bias.
foldAs = kappa .* (1 + foldXs) ./ (2 * (1 + delta)); % Convert both folds to a.
plot(ax, foldHs, foldAs, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'DisplayName', 'Saddle nodes'); % Mark folds.
xlabel(ax, 'Bias h');                          % Label continuation parameter.
ylabel(ax, 'Interface trace a');               % Label physical state.
title(ax, sprintf('S-curve at fixed \\kappa=%.2f', kappa)); % Identify slice.
legend(ax, 'Location', 'northwest');           % Explain stability and folds.
style_axis(ax);                                % Apply common style.
selectedH = 0.5 * hLimit;                      % Choose an interior three-root bias.
sGrid = linspace(0, kappa, 2500);              % Sample the full jump interval.
qInterface = beta_hat(sGrid, kappa, delta, selectedH) .* sGrid; % Compute Kapitza flux.
qBulk = delta ./ (1 + delta) .* (kappa - sGrid); % Compute bulk-required flux.
states = physical_states(kappa, delta, selectedH); % Compute all steady intersections.
ax = nexttile(layout);                         % Select the heat-flux panel.
plot(ax, sGrid, qInterface, 'Color', [0.77, 0.24, 0.17], 'LineWidth', 2, 'DisplayName', 'q_K=\beta(s)s'); % Draw interface flux.
hold(ax, 'on');                                % Retain bulk flux and roots.
plot(ax, sGrid, qBulk, 'Color', [0.12, 0.37, 0.75], 'LineWidth', 2, 'DisplayName', 'q_{bulk}'); % Draw bulk flux.
plot(ax, [states.s], [states.qHat], 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'DisplayName', 'Steady states'); % Mark roots.
xlabel(ax, 'Jump s=[u]');                      % Label the interface jump.
ylabel(ax, 'Normalized heat flux q');          % Label normalized flux.
title(ax, sprintf('Three intersections at h=%.4f', selectedH)); % State the selected bias.
legend(ax, 'Location', 'northeast');           % Explain flux laws.
style_axis(ax);                                % Apply common style.
save_figure(fig, outputDirectory, 's_curve');  % Export both formats.
end                                             % Finish the S-curve figure.


function plot_biased_cobweb(outputDirectory, delta, kappa)
%PLOT_BIASED_COBWEB Show biased intersections and two visible attracting basins.
xFold = sqrt(1 - 1 / kappa);                  % Compute the positive fold state.
hLimit = kappa * xFold - atanh(xFold);         % Compute the three-root bias limit.
h = 0.5 * hLimit;                             % Select a slice halfway to the fold.
states = physical_states(kappa, delta, h);     % Compute all three fixed points.
scale = gamma_value(kappa, delta);             % Compute the physical plotting interval.
aGrid = linspace(0, scale, 1800);              % Sample the map smoothly.
middle = states(2).a;                          % Identify the repelling basin boundary.
starts = [middle - 0.08 * scale, middle + 0.08 * scale]; % Place visible starts on both sides.
fig = figure('Color', 'w', 'Position', [90, 90, 1120, 470]); % Create paired panels.
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); % Align panels.
title(layout, sprintf('Fixed-h slice: \\delta=%.1f, \\kappa=%.2f, h=%.4f', delta, kappa, h)); % State parameters.
for panel = 1:2                                % Draw common fixed-point geometry twice.
    ax = nexttile(layout);                     % Select one panel.
    plot(ax, aGrid, kapitza_map(aGrid, kappa, delta, h), 'Color', [0.12, 0.37, 0.75], 'LineWidth', 2, 'DisplayName', 'F_{\kappa,h}(a)'); % Draw map.
    hold(ax, 'on');                            % Retain identity, roots, and paths.
    plot(ax, aGrid, aGrid, 'k--', 'LineWidth', 1.3, 'DisplayName', 'a'); % Draw identity.
    plot(ax, [states.a], [states.a], 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'DisplayName', 'Fixed points'); % Mark roots.
    xlim(ax, [0, scale]);                      % Use the physical state interval.
    ylim(ax, [0, scale]);                      % Match vertical scale.
    style_axis(ax);                            % Apply common style.
    if panel == 1                              % Label the intersection-only panel.
        xlabel(ax, 'a');                       % Label physical state.
        ylabel(ax, 'F_{\kappa,h}(a)');        % Label map output.
        title(ax, 'Three fixed-point intersections'); % State geometry.
        legend(ax, 'Location', 'northwest');   % Explain curves and roots.
    else                                       % Add the two cobweb paths to the second panel.
        colors = {[0.15, 0.55, 0.35], [0.77, 0.24, 0.17]}; % Assign path colors.
        labels = {'lower start', 'upper start'}; % Name both sides of the repeller.
        for pathIndex = 1:2                    % Draw both attracting-basin transients.
            [pathX, pathY] = cobweb_path(starts(pathIndex), kappa, delta, h, 16); % Generate path.
            plot(ax, pathX, pathY, 'Color', colors{pathIndex}, 'LineWidth', 0.95, 'DisplayName', labels{pathIndex}); % Draw path.
        end                                     % Finish the path loop.
        xlabel(ax, 'a_n');                     % Label current iterate.
        ylabel(ax, 'a_{n+1}');                 % Label next iterate.
        title(ax, 'Two attracting basins');    % State the dynamical interpretation.
        legend(ax, 'Location', 'northwest');   % Identify both starts and reference curves.
    end                                         % Finish the panel-specific content.
end                                             % Finish the two-panel loop.
save_figure(fig, outputDirectory, 'biased_fixed_point_cobweb'); % Export both formats.
end                                             % Finish the biased cobweb figure.


function plot_temperature_profiles(outputDirectory, delta, kappa)
%PLOT_TEMPERATURE_PROFILES Reconstruct the paper's piecewise steady solutions.
states = physical_states(kappa, delta, 0);     % Compute all symmetric states.
zLeft = linspace(-1, 0, 500);                 % Sample the left material interval.
zRight = linspace(0, 1, 500);                 % Sample the right material interval.
colors = {[0.15, 0.55, 0.35], [0.77, 0.24, 0.17], [0.12, 0.37, 0.75]}; % Assign branch colors.
fig = figure('Color', 'w', 'Position', [150, 100, 800, 500]); % Create one profile axis.
ax = axes(fig);                                % Obtain the plotting axis explicitly.
hold(ax, 'on');                                % Retain every branch and jump segment.
for index = 1:numel(states)                    % Reconstruct each steady solution.
    state = states(index);                     % Read one physical state record.
    uLeft = state.a .* (zLeft + 1);            % Evaluate u_I(z)=a(z+1).
    uRight = state.b .* (zRight - 1) + kappa .* (1 - zRight .^ 2); % Evaluate u_E(z).
    if state.attracting                        % Encode Picard-attracting branches as solid.
        lineStyle = '-';                       % Select a solid line.
        stability = 'attracting';              % Build a readable legend label.
    else                                       % Encode the repelling branch as dashed.
        lineStyle = '--';                      % Select a dashed line.
        stability = 'repelling';               % Build a readable legend label.
    end                                         % Finish the stability-style selection.
    plot(ax, zLeft, uLeft, lineStyle, 'Color', colors{index}, 'LineWidth', 2, 'DisplayName', sprintf('Branch %d: %s', index, stability)); % Draw left profile.
    plot(ax, zRight, uRight, lineStyle, 'Color', colors{index}, 'LineWidth', 2, 'HandleVisibility', 'off'); % Draw right profile.
    plot(ax, [0, 0], [state.a, kappa - state.b], ':', 'Color', colors{index}, 'LineWidth', 1.1, 'HandleVisibility', 'off'); % Show jump.
end                                             % Finish the profile loop.
xline(ax, 0, 'k-', 'LineWidth', 0.9, 'HandleVisibility', 'off'); % Mark material interface.
yline(ax, 0, '-', 'Color', [0.42, 0.42, 0.42], 'LineWidth', 0.8, 'HandleVisibility', 'off'); % Mark boundary level.
xlabel(ax, 'Position z');                      % Label the one-dimensional coordinate.
ylabel(ax, 'Temperature u(z)');                % Label the temperature field.
title(ax, sprintf('Piecewise solutions at \\kappa=%.2f and h=0', kappa)); % Identify slice.
text(ax, -0.52, 1.52, '\Omega_I', 'HorizontalAlignment', 'center', 'FontSize', 13); % Label left layer.
text(ax, 0.52, 1.52, '\Omega_E', 'HorizontalAlignment', 'center', 'FontSize', 13); % Label right layer.
legend(ax, 'Location', 'northeast');           % Identify all three branches.
style_axis(ax);                                % Apply common style.
save_figure(fig, outputDirectory, 'temperature_profiles'); % Export both formats.
end                                             % Finish the PDE-profile figure.


function plot_model_summary(outputDirectory, delta, kappa)
%PLOT_MODEL_SUMMARY Create a compact four-panel source-model overview.
fig = figure('Color', 'w', 'Position', [90, 50, 1000, 800]); % Create a square summary figure.
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); % Align panels.
title(layout, sprintf('Source-controlled nonlinear Kapitza model (\\delta=%.1f, \\kappa=%.2f)', delta, kappa)); % State parameters.
scale = gamma_value(kappa, delta);             % Compute the physical state interval.
aGrid = linspace(0, scale, 1600);              % Sample the fixed-point map.
states = physical_states(kappa, delta, 0);     % Compute all symmetric roots.
ax = nexttile(layout);                         % Select the fixed-point panel.
plot(ax, aGrid, kapitza_map(aGrid, kappa, delta, 0), 'Color', [0.12, 0.37, 0.75], 'LineWidth', 2, 'DisplayName', 'F_\kappa(a)'); % Draw map.
hold(ax, 'on');                                % Retain identity and roots.
plot(ax, aGrid, aGrid, '--', 'Color', [0.77, 0.24, 0.17], 'LineWidth', 1.4, 'DisplayName', 'a'); % Draw identity.
plot(ax, [states.a], [states.a], 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'HandleVisibility', 'off'); % Mark roots.
xlabel(ax, 'a'); ylabel(ax, 'F_\kappa(a)'); title(ax, 'Fixed-point intersections'); % Label panel.
legend(ax, 'Location', 'best'); style_axis(ax); % Explain curves and style axis.
sGrid = linspace(0, kappa, 1600);              % Sample the complete jump interval.
ax = nexttile(layout);                         % Select the conductance panel.
semilogy(ax, sGrid, beta_hat(sGrid, kappa, delta, 0), 'Color', [0.12, 0.37, 0.75], 'LineWidth', 2); % Draw beta.
xlabel(ax, 'Jump s'); ylabel(ax, '\beta_{\kappa,0}(s)'); title(ax, 'Normalized Kapitza conductance'); style_axis(ax); % Label panel.
ax = nexttile(layout);                         % Select the flux panel.
plot(ax, sGrid, beta_hat(sGrid, kappa, delta, 0) .* sGrid, 'Color', [0.77, 0.24, 0.17], 'LineWidth', 2, 'DisplayName', 'q_K'); % Draw interface flux.
hold(ax, 'on');                                % Retain bulk flux and roots.
plot(ax, sGrid, delta ./ (1 + delta) .* (kappa - sGrid), 'Color', [0.12, 0.37, 0.75], 'LineWidth', 2, 'DisplayName', 'q_{bulk}'); % Draw bulk flux.
plot(ax, [states.s], [states.qHat], 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'HandleVisibility', 'off'); % Mark intersections.
xlabel(ax, 'Jump s'); ylabel(ax, 'Normalized flux q'); title(ax, 'Flux intersections'); legend(ax, 'Location', 'best'); style_axis(ax); % Label panel.
ax = nexttile(layout);                         % Select the source-support panel.
fill(ax, [0, 1, 1, 0], [0, 0, 2 * kappa, 2 * kappa], [0.12, 0.37, 0.75], 'FaceAlpha', 0.18, 'EdgeColor', 'none'); % Shade source support.
hold(ax, 'on');                                % Retain the source-amplitude line.
plot(ax, [0, 1], [2 * kappa, 2 * kappa], 'Color', [0.12, 0.37, 0.75], 'LineWidth', 2); % Draw source amplitude.
xline(ax, 0, 'k-', 'LineWidth', 0.9);          % Mark the material interface.
xlim(ax, [-1, 1]); ylim(ax, [0, 2.25 * kappa]); % Frame both layers and source height.
xlabel(ax, 'Position z'); ylabel(ax, 'f(z)/\mu'); title(ax, 'Source 2\kappa\chi_{\Omega_E}'); style_axis(ax); % Label panel.
save_figure(fig, outputDirectory, 'model_summary'); % Export both formats.
end                                             % Finish the summary figure.


function write_state_table(outputDirectory, delta, kappas)
%WRITE_STATE_TABLE Save selected steady states and stability data to CSV.
rows = [];                                     % Initialize a numeric row matrix.
for kappa = kappas                             % Process every requested source parameter.
    states = physical_states(kappa, delta, 0); % Compute all symmetric states.
    for branch = 1:numel(states)               % Serialize every branch independently.
        state = states(branch);                % Read one state record.
        rows = [rows; kappa, 0, branch, state.x, state.a, state.b, state.s, ... %#ok<AGROW>
            state.betaHat, state.qHat, state.multiplier, state.attracting]; % Append one row.
    end                                         % Finish the branch loop.
end                                             % Finish the source-parameter loop.
tableData = array2table(rows, 'VariableNames', {'kappa', 'h', 'branch', 'x', 'a', 'b', 's', 'beta_hat', 'q_hat', 'multiplier', 'attracting'}); % Add headers.
writetable(tableData, fullfile(outputDirectory, 'steady_states.csv')); % Write the portable table.
end                                             % Finish the state-table writer.


function write_fold_table(outputDirectory, delta)
%WRITE_FOLD_TABLE Save representative positive and negative saddle nodes.
biases = [-0.04, -0.02, -0.01, -0.005, -0.002, 0, 0.002, 0.005, 0.01, 0.02, 0.04]; % Select samples.
rows = zeros(numel(biases), 5);                % Preallocate one row per bias.
for index = 1:numel(biases)                   % Evaluate every exact fold.
    h = biases(index);                         % Read the current bias.
    [xFold, kappaFold] = fold_for_bias(h);     % Solve the fold equations.
    aFold = kappaFold * (1 + xFold) / (2 * (1 + delta)); % Convert to physical state.
    if h == 0                                  % Evaluate the local approximation at the origin.
        approximation = 0;                    % Use the exact limiting displacement.
    else                                       % Evaluate the nonzero local asymptotic law.
        approximation = (1.5 * abs(h)) ^ (2 / 3); % Compute the two-thirds prediction.
    end                                         % Finish the zero/nonzero split.
    rows(index, :) = [h, xFold, kappaFold, aFold, approximation]; % Store one fold row.
end                                             % Finish the fold loop.
tableData = array2table(rows, 'VariableNames', {'h', 'x_SN', 'kappa_SN', 'a_SN', 'local_asymptotic_shift'}); % Add headers.
writetable(tableData, fullfile(outputDirectory, 'critical_trajectory.csv')); % Write the table.
end                                             % Finish the fold-table writer.


function style_axis(ax)                        % Apply one visual style to a MATLAB axis.
%STYLE_AXIS Add a subtle grid, box, and inward ticks.
grid(ax, 'on');                                % Enable the reading grid.
ax.GridAlpha = 0.18;                           % Keep grid lines visually secondary.
ax.GridColor = [0.45, 0.45, 0.45];            % Use a neutral gray grid.
ax.TickDir = 'in';                             % Draw compact inward ticks.
ax.Box = 'on';                                 % Frame each panel clearly.
ax.Layer = 'top';                              % Keep ticks visible above plotted objects.
end                                             % Finish the axis-style helper.


function save_figure(fig, outputDirectory, stem)
%SAVE_FIGURE Export one MATLAB figure as separated PNG and vector PDF files.
drawnow;                                       % Force all text and graphics to finish rendering.
pngPath = fullfile(outputDirectory, [stem, '.png']); % Build the raster destination.
pdfPath = fullfile(outputDirectory, [stem, '.pdf']); % Build the vector destination.
exportgraphics(fig, pngPath, 'Resolution', 240); % Save a high-resolution PNG figure.
set(fig, 'PaperPositionMode', 'auto');          % Preserve the on-screen dimensions in PDF output.
print(fig, pdfPath, '-dpdf', '-painters', '-bestfit'); % Save a vector PDF through the stable print path.
close(fig);                                    % Release the completed figure from memory.
end                                             % Finish the figure-export helper.


function outputDirectory = kapitza_output_directory()
%KAPITZA_OUTPUT_DIRECTORY Resolve output/matlab relative to this script.
scriptPath = mfilename('fullpath');             % Obtain the absolute script path.
matlabDirectory = fileparts(scriptPath);        % Remove the script filename.
codeDirectory = fileparts(matlabDirectory);     % Move from matlab/ to code/.
projectDirectory = fileparts(codeDirectory);    % Move from code/ to the project root.
outputDirectory = fullfile(projectDirectory, 'output', 'matlab'); % Build destination.
if ~isfolder(outputDirectory)                   % Check whether the destination already exists.
    mkdir(outputDirectory);                     % Create it recursively when necessary.
end                                             % Finish the directory-existence check.
end                                             % Finish the output-path helper.
