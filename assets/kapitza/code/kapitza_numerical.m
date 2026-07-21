%KAPITZA_NUMERICAL Complete MATLAB study of the nonlinear Kapitza problem.
% This single script reproduces the teacher's fixed-point example, compares
% the three alpha regimes, draws cobweb paths, constructs the pitchfork, and
% unfolds it into a biased S-curve. The teacher's short code calls the material
% parameter d; this script uses the paper notation delta for the same quantity.

clearvars;                                      % Remove variables left by earlier MATLAB sessions.
close all;                                      % Close figures left by earlier MATLAB sessions.
clc;                                            % Clear the Command Window for readable output.

fprintf('Running the complete nonlinear Kapitza workflow.\n'); % Announce the batch calculation.
teacher_fixed_point_demo();                     % Reproduce the teacher's subcritical example.
fixed_point_regimes_demo();                     % Compare intersections below, at, and above alpha_c.
cobweb_regimes_demo();                          % Compare fixed-point iteration in all three regimes.
pitchfork_bifurcation_demo();                   % Draw the complete symmetric pitchfork diagram.
small_bias_alpha_bifurcation_demo();            % Unfold the pitchfork at small fixed h with alpha varying.
small_negative_bias_alpha_bifurcation_demo();   % Show the mirror unfolding for small negative h values.
critical_point_trajectory_demo();               % Trace the saddle-node path from the pitchfork limit.
s_curve_demo();                                 % Draw the biased S-curve and heat-flux intersections.
s_curve_fixed_point_cobweb_demo();              % Draw biased fixed points and their cobweb basins.
fixed_a_h_iteration_demo();                     % Recover h for one prescribed steady state a_star.
fprintf('All MATLAB outputs were written to output/matlab/.\n'); % Confirm successful completion.


function value = kapitza_beta(s, delta, alpha, h) % Evaluate the nonlinear Kapitza conductance.
%KAPITZA_BETA Use the teacher's alpha convention, whose threshold is alpha_c=2.
if nargin < 4                                  % Check whether the optional bias was omitted.
    h = 0;                                     % Use the symmetric conductance law by default.
end                                            % Finish the optional-argument check.
value = delta ./ (1 + delta) .* ...            % Apply the material prefactor element by element.
    exp(alpha .* (1 - 2 .* s) + 2 .* h);       % Apply the exponential jump dependence and bias.
end                                            % Finish the conductance function.


function value = kapitza_map(a, delta, alpha, h) % Evaluate the scalar fixed-point map F(a).
%KAPITZA_MAP Reconstruct the jump and apply the nonlinear interface law.
if nargin < 4                                  % Check whether the optional bias was omitted.
    h = 0;                                     % Use the symmetric interface law by default.
end                                            % Finish the optional-argument check.
s = 1 - (1 + delta) .* a;                      % Convert the affine state a to the jump s=[u].
b = kapitza_beta(s, delta, alpha, h);           % Evaluate the conductance at the reconstructed jump.
value = b ./ (delta + (1 + delta) .* b);        % Evaluate F(a)=beta/[delta+(1+delta)beta].
end                                            % Finish the fixed-point map.


function xRoots = kapitza_roots(alpha, h)      % Find every normalized steady-state root.
%KAPITZA_ROOTS Solve x=tanh(alpha*x/2+h), retaining double roots at folds.
if nargin < 2                                  % Check whether the optional bias was omitted.
    h = 0;                                     % Solve the symmetric problem by default.
end                                            % Finish the optional-argument check.
kappa = alpha / 2;                             % Convert alpha to the normalized sigmoid slope.
residual = @(x) tanh(kappa .* x + h) - x;      % Define the scalar root residual.
edge = 1e-12;                                  % Stay slightly inside the physical interval (-1,1).
tolerance = 2e-10;                             % Set the tolerance for roots located at folds.
if kappa > 1                                   % Check whether stationary points can exist.
    xFold = sqrt(1 - 1 / kappa);               % Solve residual'(x)=0 analytically.
    splitPoints = [-1 + edge, -xFold, xFold, 1 - edge]; % Split into monotone subintervals.
else                                           % Handle the monotone subcritical or critical case.
    splitPoints = [-1 + edge, 1 - edge];       % Search the whole physical interval at once.
end                                            % Finish construction of search intervals.
candidates = [];                               % Initialize the candidate-root list.
for j = 1:(numel(splitPoints) - 1)             % Inspect each monotone subinterval.
    left = splitPoints(j);                     % Read the left endpoint.
    right = splitPoints(j + 1);                % Read the right endpoint.
    fLeft = residual(left);                    % Evaluate the residual at the left endpoint.
    fRight = residual(right);                  % Evaluate the residual at the right endpoint.
    if abs(fLeft) < tolerance                  % Detect a root exactly at the left endpoint.
        candidates(end + 1) = left; %#ok<AGROW> % Retain a possible double root at a fold.
    end                                        % Finish the left-endpoint test.
    if fLeft * fRight < 0                      % Test whether the residual changes sign.
        candidates(end + 1) = fzero(residual, [left, right]); %#ok<AGROW> % Refine the root.
    end                                        % Finish the sign-change test.
    if abs(fRight) < tolerance                 % Detect a root exactly at the right endpoint.
        candidates(end + 1) = right; %#ok<AGROW> % Retain a possible double root at a fold.
    end                                        % Finish the right-endpoint test.
end                                            % Finish the interval loop.
candidates = sort(candidates);                 % Arrange candidates in increasing order.
xRoots = [];                                   % Initialize the duplicate-free root list.
for value = candidates                         % Examine each sorted candidate.
    if isempty(xRoots) || abs(value - xRoots(end)) > 1e-8 % Reject duplicate fold endpoints.
        xRoots(end + 1) = value; %#ok<AGROW>    % Keep one copy of each distinct root.
    end                                        % Finish the duplicate check.
end                                            % Finish the candidate loop.
end                                            % Finish the root solver.


function [pathX, pathY] = kapitza_cobweb_path(a0, delta, alpha, h, iterations)
%KAPITZA_COBWEB_PATH Build one vertical-horizontal cobweb polyline.
if nargin < 4                                  % Check whether the optional bias was omitted.
    h = 0;                                     % Use the symmetric map by default.
end                                            % Finish the optional-bias check.
if nargin < 5                                  % Check whether the iteration count was omitted.
    iterations = 30;                           % Use a readable default number of steps.
end                                            % Finish the optional-count check.
pathX = a0;                                    % Start the polyline at the initial state.
pathY = 0;                                     % Start from the horizontal axis.
aCurrent = a0;                                 % Store the current fixed-point iterate.
for j = 1:iterations                           % Generate every vertical-horizontal pair.
    aNext = kapitza_map(aCurrent, delta, alpha, h); % Apply the nonlinear map once.
    pathX = [pathX, aCurrent, aNext]; %#ok<AGROW> % Append vertical and horizontal x-coordinates.
    pathY = [pathY, aNext, aNext]; %#ok<AGROW> % Append the matching y-coordinates.
    aCurrent = aNext;                          % Advance to the next iterate.
end                                            % Finish the cobweb loop.
end                                            % Finish the cobweb helper.


function teacher_fixed_point_demo()            % Reproduce the teacher's fixed-point experiment.
%TEACHER_FIXED_POINT_DEMO Translate and extend the original Scilab code.
delta = 0.5;                                   % Use delta; the teacher's short code calls it d.
alpha = 1.5;                                   % Select the teacher's subcritical alpha value.
h = 0;                                         % Use the symmetric Kapitza law.
gamma = 1 / (1 + delta);                       % Compute the upper physical bound for a.
a = linspace(0, gamma, 2000);                  % Sample the complete physical interval.
F = kapitza_map(a, delta, alpha, h);            % Evaluate the fixed-point map on the grid.
figure('Color', 'w');                          % Open a white-background figure.
plot(a, F, 'b-', 'LineWidth', 2);              % Draw the nonlinear map in blue.
hold on;                                       % Keep the axes for the identity and roots.
plot(a, a, 'r--', 'LineWidth', 1.6);           % Draw the identity line in dashed red.
grid on;                                       % Add a grid for reading the intersection.
xlabel('\ita');                                % Label the physical state in mathematical italics.
ylabel('\itF(a)');                             % Label the mapped state in mathematical italics.
title(sprintf('Kapitza fixed points: \\delta = %.1f, \\alpha = %.1f', delta, alpha)); % Build title.
legend('\itF(a)', '\ita', 'Location', 'northwest', 'AutoUpdate', 'off'); % Identify both curves.
xRoots = kapitza_roots(alpha, h);               % Find every normalized root.
aRoots = (1 + xRoots) ./ (2 * (1 + delta));     % Convert roots to physical states a.
plot(aRoots, aRoots, 'ko', 'MarkerFaceColor', 'k'); % Mark each fixed-point intersection.
fprintf('\nTeacher example: delta=%.6g, alpha=%.6g, h=%.6g\n', delta, alpha, h); % Print inputs.
fprintf('Number of fixed points: %d\n', numel(aRoots)); % Print the root count.
fprintf('        a                 s=[u]              beta(s)             q=delta*a\n'); % Print headers.
for j = 1:numel(aRoots)                        % Visit every detected fixed point.
    s = 1 - (1 + delta) * aRoots(j);            % Reconstruct the temperature jump.
    b = kapitza_beta(s, delta, alpha, h);       % Evaluate the Kapitza conductance.
    fprintf('% .12f    % .12f    % .12f    % .12f\n', ... % Select numeric formatting.
        aRoots(j), s, b, delta * aRoots(j));    % Print state, jump, conductance, and flux.
end                                            % Finish the root-reporting loop.
kapitza_save_figure(gcf, 'teacher_fixed_point.png'); % Save the completed figure.
end                                            % Finish the teacher example.


function fixed_point_regimes_demo()            % Compare fixed-point intersections across regimes.
%FIXED_POINT_REGIMES_DEMO Show subcritical, critical, and supercritical maps.
delta = 0.5;                                   % Use delta; the teacher's short code calls it d.
alphaValues = [1.5, 2.0, 3.0];                 % Select cases below, at, and above alpha_c=2.
gamma = 1 / (1 + delta);                       % Compute the upper physical bound for a.
a = linspace(0, gamma, 2000);                  % Sample the full physical interval.
figure('Color', 'w');                          % Open a white-background comparison figure.
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact'); % Arrange three panels.
for j = 1:numel(alphaValues)                   % Visit all three bifurcation regimes.
    alpha = alphaValues(j);                    % Read the current alpha value.
    F = kapitza_map(a, delta, alpha, 0);        % Evaluate the symmetric fixed-point map.
    xRoots = kapitza_roots(alpha, 0);           % Find all normalized fixed points.
    aRoots = (1 + xRoots) ./ (2 * (1 + delta)); % Convert roots to the physical variable a.
    nexttile;                                  % Activate the current panel.
    plot(a, F, 'b-', 'LineWidth', 2);          % Draw the nonlinear map.
    hold on;                                   % Keep the axes for the identity and roots.
    plot(a, a, 'r--', 'LineWidth', 1.5);       % Draw the identity line.
    plot(aRoots, aRoots, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5); % Mark intersections.
    grid on;                                   % Add a comparison grid.
    xlabel('\ita');                            % Label the state variable in mathematical italics.
    ylabel('\itF(a)');                         % Label the map value in mathematical italics.
    title(sprintf('\\alpha = %.1f: %d fixed point(s)', alpha, numel(aRoots)), ... % Build title.
        'FontWeight', 'normal');               % Use the same weight in both export formats.
    xlim([0, gamma]);                          % Apply a common horizontal range.
    ylim([0, gamma]);                          % Apply a common vertical range.
    if j == 1                                 % Add the shared legend only once.
        legend('\itF(a)', '\ita', 'Fixed points', 'Location', 'northwest'); % Identify styles.
    end                                        % Finish the one-time legend block.
end                                            % Finish the regime loop.
sgtitle('Fixed-point geometry across the pitchfork threshold', ... % Summarize the figure.
    'FontName', 'Times New Roman', 'Interpreter', 'tex'); % Match all other figure text.
kapitza_save_figure(gcf, 'fixed_point_regimes.png'); % Save the comparison figure.
end                                            % Finish the intersection comparison.


function cobweb_regimes_demo()                 % Compare iterations across the three regimes.
%COBWEB_REGIMES_DEMO Show convergence and branch selection with cobweb paths.
delta = 0.5;                                   % Use delta; the teacher's short code calls it d.
alphaValues = [1.5, 2.0, 3.0];                 % Select subcritical, critical, and supercritical cases.
gamma = 1 / (1 + delta);                       % Compute the physical interval length.
a = linspace(0, gamma, 2000);                  % Sample all maps on the same interval.
defaultInitialValues = [0.12, 0.88] .* gamma;  % Keep wide starts for the one-root regimes.
pathColors = [0.10, 0.55, 0.30; 0.75, 0.20, 0.15]; % Assign colors to both initial states.
figure('Color', 'w');                          % Open a white-background comparison figure.
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact'); % Arrange three panels.
for j = 1:numel(alphaValues)                   % Visit each bifurcation regime.
    alpha = alphaValues(j);                    % Read the current alpha value.
    F = kapitza_map(a, delta, alpha, 0);        % Evaluate the symmetric fixed-point map.
    xRoots = kapitza_roots(alpha, 0);           % Find all normalized fixed points.
    aRoots = (1 + xRoots) ./ (2 * (1 + delta)); % Convert roots to physical states.
    if alpha > 2                               % Detect the three-root supercritical panel.
        middleFixedPoint = aRoots(2);           % Read the repelling central basin boundary.
        initialOffset = 0.08 * gamma;           % Stay close enough to reveal departure from it.
        panelInitialValues = middleFixedPoint + ... % Place one initial condition in each basin.
            [-initialOffset, initialOffset];    % Straddle the repelling central fixed point.
    else                                       % Preserve the one-root and critical comparisons.
        panelInitialValues = defaultInitialValues; % Use starts spanning most of the interval.
    end                                        % Finish the panel-specific start selection.
    nexttile;                                  % Activate the current cobweb panel.
    plot(a, F, 'b-', 'LineWidth', 1.8);        % Draw the nonlinear map.
    hold on;                                   % Keep the axes for all remaining elements.
    plot(a, a, 'k--', 'LineWidth', 1.2);       % Draw the identity line.
    for k = 1:numel(panelInitialValues)        % Draw a cobweb from each selected side.
        [pathX, pathY] = kapitza_cobweb_path(panelInitialValues(k), delta, alpha, 0, 35); % Build a path.
        plot(pathX, pathY, '-', 'Color', pathColors(k, :), 'LineWidth', 0.9); % Draw the path.
    end                                        % Finish the initial-condition loop.
    plot(aRoots, aRoots, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4); % Mark fixed points.
    grid on;                                   % Add a grid for reading iterations.
    xlabel('\ita_n');                          % Label the current iterate in mathematical italics.
    ylabel('\ita_{n+1}');                      % Label the next iterate in mathematical italics.
    title(sprintf('\\alpha = %.1f', alpha), ... % Build the current-parameter title.
        'FontWeight', 'normal');               % Match both export formats.
    xlim([0, gamma]);                          % Apply a common horizontal range.
    ylim([0, gamma]);                          % Apply a common vertical range.
    if j == 1                                 % Add the shared legend only once.
        legend('\itF(a)', '\ita', 'lower start', 'upper start', 'Fixed points', ... % Name elements.
            'Location', 'northwest');          % Place the legend inside the first panel.
    end                                        % Finish the one-time legend block.
end                                            % Finish the regime loop.
sgtitle('Cobweb iterations reveal stability and branch selection', ... % Summarize the figure.
    'FontName', 'Times New Roman', 'Interpreter', 'tex'); % Match all other figure text.
kapitza_save_figure(gcf, 'cobweb_regimes.png'); % Save the cobweb comparison.
end                                            % Finish the cobweb comparison.


function pitchfork_bifurcation_demo()          % Plot the analytic symmetric pitchfork.
%PITCHFORK_BIFURCATION_DEMO Use the teacher's universal threshold alpha_c=2.
delta = 0.5;                                   % Use delta; the teacher's short code calls it d.
gamma = 1 / (1 + delta);                       % Compute the physical interval length.
alphaMax = 5;                                  % Set the largest displayed alpha value.
xPositive = linspace(1e-4, 0.985, 1600);       % Parameterize the positive outer branch.
alphaOuter = 2 .* atanh(xPositive) ./ xPositive; % Evaluate the analytic continuation formula.
keep = alphaOuter <= alphaMax;                 % Retain points inside the displayed range.
figure('Color', 'w');                          % Open a white-background figure.
hCentralStable = plot([0.25, 2], [gamma / 2, gamma / 2], 'r-', 'LineWidth', 2); % Stable center.
hold on;                                       % Keep the axes for all remaining branches.
hCentralUnstable = plot([2, alphaMax], [gamma / 2, gamma / 2], ... % Select unstable center.
    'r--', 'LineWidth', 2);                    % Draw the unstable center as dashed.
hOuter = plot(alphaOuter(keep), gamma / 2 .* (1 + xPositive(keep)), ... % Select upper branch.
    'b-', 'LineWidth', 2);                     % Draw the stable upper branch as solid.
plot(alphaOuter(keep), gamma / 2 .* (1 - xPositive(keep)), ... % Select lower branch.
    'b-', 'LineWidth', 2, 'HandleVisibility', 'off'); % Draw it without a duplicate legend item.
hCritical = xline(2, 'k--', 'LineWidth', 1.3); % Mark the nonhyperbolic critical value.
grid on;                                       % Add a grid for reading branch values.
xlabel('\alpha');                             % Label the continuation parameter with TeX.
ylabel('\ita');                                % Label the physical state in mathematical italics.
title('Supercritical pitchfork: \delta = 0.5, \alpha_c = 2'); % State parameters.
legend([hCentralStable, hCentralUnstable, hOuter, hCritical], ... % Select style handles.
    {'Stable central branch', 'Unstable central branch', ... % Name stable and unstable branches.
     'Stable outer branches', 'Critical value'}, 'Location', 'east'); % Complete the legend.
xlim([0.25, alphaMax]);                        % Apply the selected alpha range.
ylim([0, gamma]);                              % Restrict a to its physical interval.
kapitza_save_figure(gcf, 'pitchfork_bifurcation.png'); % Save the pitchfork figure.
end                                            % Finish the pitchfork demo.


function small_bias_alpha_bifurcation_demo()  % Compare six fixed-h continuations in alpha.
%SMALL_BIAS_ALPHA_BIFURCATION_DEMO Display progressive imperfect unfolding.
delta = 0.5;                                   % Use delta; the teacher's short code calls it d.
hValues = [0, 0.002, 0.005, 0.010, 0.020, 0.040]; % Select six small biases including baseline.
gamma = 1 / (1 + delta);                       % Compute the physical interval length for a.
alphaMin = 0.25;                               % Start below the symmetric threshold.
alphaMax = 5;                                  % Continue far enough to display all branches.
stableColor = [0.09, 0.34, 0.79];             % Use blue for attracting fixed-point branches.
unstableColor = [0.76, 0.23, 0.14];           % Use red for repelling fixed-point branches.
figure('Color', 'w');                          % Open a white-background comparison figure.
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact'); % Arrange six matched panels.

for panel = 1:numel(hValues)                  % Visit every selected fixed bias.
    hSmall = hValues(panel);                   % Read the bias for the current continuation slice.
    nexttile;                                  % Activate the corresponding panel.
    if hSmall == 0                             % Draw the exact symmetric pitchfork as baseline.
        alphaCentralStable = linspace(alphaMin, 2, 400); % Sample the stable central segment.
        alphaCentralUnstable = linspace(2, alphaMax, 600); % Sample the unstable continuation.
        xOuter = linspace(1e-4, 0.995, 2200);  % Parameterize both symmetric outer branches.
        alphaOuter = 2 .* atanh(xOuter) ./ xOuter; % Recover alpha along the outer branches.
        outerVisible = alphaOuter <= alphaMax; % Keep only the displayed continuation interval.
        hStable = plot(alphaCentralStable, 0.5 * gamma .* ones(size(alphaCentralStable)), ... % Center.
            '-', 'Color', stableColor, 'LineWidth', 1.8); % Use a solid line for attraction.
        hold on;                               % Keep the axes for all other branches and points.
        hUnstable = plot(alphaCentralUnstable, 0.5 * gamma .* ones(size(alphaCentralUnstable)), ... % Center.
            '--', 'Color', unstableColor, 'LineWidth', 1.8); % Use a dashed line for repulsion.
        plot(alphaOuter(outerVisible), 0.5 * gamma .* (1 + xOuter(outerVisible)), ... % Upper branch.
            '-', 'Color', stableColor, 'LineWidth', 1.8, 'HandleVisibility', 'off'); % Stable.
        plot(alphaOuter(outerVisible), 0.5 * gamma .* (1 - xOuter(outerVisible)), ... % Lower branch.
            '-', 'Color', stableColor, 'LineWidth', 1.8, 'HandleVisibility', 'off'); % Stable.
        hCritical = plot(2, 0.5 * gamma, 'ko', 'MarkerFaceColor', 'k', ... % Mark pitchfork.
            'MarkerSize', 4.5);                % Keep the marker legible in the small panel.
        title('h = 0: perfect pitchfork', 'FontWeight', 'normal'); % Identify the baseline.
        legend([hStable, hUnstable, hCritical], ... % Explain the shared visual encoding once.
            {'Attracting', 'Repelling', 'Critical point'}, 'Location', 'east'); % Place legend.
        fprintf('\nSmall-bias alpha continuation: h=0, alpha_c=2\n'); % Report baseline threshold.
    else                                       % Draw the imperfect branches for nonzero h.
        xNegative = linspace(-0.995, -1e-4, 5000); % Parameterize both negative-state branches.
        xPositive = linspace(tanh(hSmall) + 1e-6, 0.995, 3000); % Parameterize favored branch.
        alphaNegative = 2 .* (atanh(xNegative) - hSmall) ./ xNegative; % Continue negative roots.
        alphaPositive = 2 .* (atanh(xPositive) - hSmall) ./ xPositive; % Continue positive roots.
        lambdaNegative = 0.5 .* alphaNegative .* (1 - xNegative .^ 2); % Evaluate multipliers.
        negativeVisible = alphaNegative >= alphaMin & alphaNegative <= alphaMax; % Clip alpha.
        negativeStable = negativeVisible & lambdaNegative < 1; % Select stable lower branch.
        negativeUnstable = negativeVisible & lambdaNegative >= 1; % Select repelling branch.
        positiveVisible = alphaPositive >= alphaMin & alphaPositive <= alphaMax; % Clip favored branch.
        plot(alphaPositive(positiveVisible), ... % Draw the smoothly continued favored branch.
            0.5 * gamma .* (1 + xPositive(positiveVisible)), '-', ... % Convert x to a.
            'Color', stableColor, 'LineWidth', 1.8); % Mark this branch as attracting.
        hold on;                               % Keep the axes for the saddle-node pair.
        plot(alphaNegative(negativeStable), ... % Draw the stable branch born at the fold.
            0.5 * gamma .* (1 + xNegative(negativeStable)), '-', ... % Convert x to a.
            'Color', stableColor, 'LineWidth', 1.8); % Mark the branch as attracting.
        plot(alphaNegative(negativeUnstable), ... % Draw the middle branch approaching x=0.
            0.5 * gamma .* (1 + xNegative(negativeUnstable)), '--', ... % Convert x to a.
            'Color', unstableColor, 'LineWidth', 1.8); % Mark this branch as repelling.
        foldResidual = @(x) atanh(x) - x ./ (1 - x .^ 2) - hSmall; % Combine fold equations.
        xFold = fzero(foldResidual, [-0.999999, -1e-12]); % Isolate the unique negative fold.
        alphaFold = 2 / (1 - xFold ^ 2);       % Recover the saddle-node alpha value.
        aFold = 0.5 * gamma * (1 + xFold);     % Convert the fold state to physical a.
        plot(alphaFold, aFold, 'ko', 'MarkerFaceColor', 'k', ... % Mark the saddle node.
            'MarkerSize', 4.5);                % Keep the marker clear but compact.
        title(sprintf('h = %.3f, \\alpha_{SN} = %.3f', hSmall, alphaFold), ... % State fold shift.
            'FontWeight', 'normal');           % Keep title weight consistent across panels.
        fprintf('Small-bias alpha continuation: h=%.6g, alpha_SN=%.12g, x_SN=%.12g\n', ... % Report.
            hSmall, alphaFold, xFold);         % Print the current fold coordinates.
    end                                        % Finish the symmetric/imperfect branch choice.
    grid on;                                   % Add the same reading grid to every panel.
    xlabel('\alpha');                        % Keep alpha identified as the bifurcation parameter.
    ylabel('\ita');                          % Label the physical steady-state variable.
    xlim([alphaMin, alphaMax]);                % Apply one common alpha scale.
    ylim([0, gamma]);                          % Apply one common physical state scale.
end                                            % Finish the six-bias comparison loop.

sgtitle('Progressive unfolding at six small fixed biases (\delta = 0.5; \alpha is the control parameter)', ... % Main result.
    'FontName', 'Times New Roman', 'Interpreter', 'tex'); % Use the common typography.
kapitza_save_figure(gcf, 'small_bias_alpha_bifurcation.png'); % Save PNG and vector PDF.
end                                            % Finish the six-bias alpha continuation.


function small_negative_bias_alpha_bifurcation_demo() % Compare six negative fixed-h continuations.
%SMALL_NEGATIVE_BIAS_ALPHA_BIFURCATION_DEMO Display the reflected imperfect unfolding.
delta = 0.5;                                   % Use delta; the teacher's short code calls it d.
hValues = [0, -0.002, -0.005, -0.010, -0.020, -0.040]; % Select the reflected bias sequence.
gamma = 1 / (1 + delta);                       % Compute the physical interval length for a.
alphaMin = 0.25;                               % Start below the symmetric threshold.
alphaMax = 5;                                  % Continue far enough to display all branches.
stableColor = [0.09, 0.34, 0.79];             % Use blue for attracting fixed-point branches.
unstableColor = [0.76, 0.23, 0.14];           % Use red for repelling fixed-point branches.
figure('Color', 'w');                          % Open a white-background comparison figure.
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact'); % Arrange six matched panels.

for panel = 1:numel(hValues)                  % Visit every selected fixed bias.
    hSmall = hValues(panel);                   % Read the nonpositive bias for this slice.
    nexttile;                                  % Activate the corresponding panel.
    if hSmall == 0                             % Draw the exact symmetric pitchfork as baseline.
        alphaCentralStable = linspace(alphaMin, 2, 400); % Sample the stable central segment.
        alphaCentralUnstable = linspace(2, alphaMax, 600); % Sample the unstable continuation.
        xOuter = linspace(1e-4, 0.995, 2200);  % Parameterize both symmetric outer branches.
        alphaOuter = 2 .* atanh(xOuter) ./ xOuter; % Recover alpha along the outer branches.
        outerVisible = alphaOuter <= alphaMax; % Keep only the displayed continuation interval.
        hStable = plot(alphaCentralStable, 0.5 * gamma .* ones(size(alphaCentralStable)), ... % Center.
            '-', 'Color', stableColor, 'LineWidth', 1.8); % Use a solid line for attraction.
        hold on;                               % Keep the axes for all other branches and points.
        hUnstable = plot(alphaCentralUnstable, 0.5 * gamma .* ones(size(alphaCentralUnstable)), ... % Center.
            '--', 'Color', unstableColor, 'LineWidth', 1.8); % Use a dashed line for repulsion.
        plot(alphaOuter(outerVisible), 0.5 * gamma .* (1 + xOuter(outerVisible)), ... % Upper branch.
            '-', 'Color', stableColor, 'LineWidth', 1.8, 'HandleVisibility', 'off'); % Stable.
        plot(alphaOuter(outerVisible), 0.5 * gamma .* (1 - xOuter(outerVisible)), ... % Lower branch.
            '-', 'Color', stableColor, 'LineWidth', 1.8, 'HandleVisibility', 'off'); % Stable.
        hCritical = plot(2, 0.5 * gamma, 'ko', 'MarkerFaceColor', 'k', ... % Mark pitchfork.
            'MarkerSize', 4.5);                % Keep the marker legible in the small panel.
        title('h = 0: perfect pitchfork', 'FontWeight', 'normal'); % Identify the baseline.
        legend([hStable, hUnstable, hCritical], ... % Explain the visual encoding once.
            {'Attracting', 'Repelling', 'Critical point'}, 'Location', 'east'); % Place legend.
        fprintf('\nSmall-negative-bias alpha continuation: h=0, alpha_c=2\n'); % Report baseline.
    else                                       % Draw the imperfect branches for negative h.
        xNegative = linspace(-0.995, tanh(hSmall) - 1e-6, 3000); % Parameterize favored branch.
        xPositive = linspace(1e-4, 0.995, 5000); % Parameterize the upper saddle-node pair.
        alphaNegative = 2 .* (atanh(xNegative) - hSmall) ./ xNegative; % Continue favored roots.
        alphaPositive = 2 .* (atanh(xPositive) - hSmall) ./ xPositive; % Continue upper roots.
        lambdaPositive = 0.5 .* alphaPositive .* (1 - xPositive .^ 2); % Evaluate multipliers.
        negativeVisible = alphaNegative >= alphaMin & alphaNegative <= alphaMax; % Clip alpha.
        positiveVisible = alphaPositive >= alphaMin & alphaPositive <= alphaMax; % Clip alpha.
        positiveStable = positiveVisible & lambdaPositive < 1; % Select stable upper branch.
        positiveUnstable = positiveVisible & lambdaPositive >= 1; % Select repelling branch.
        plot(alphaNegative(negativeVisible), ... % Draw the smoothly continued favored branch.
            0.5 * gamma .* (1 + xNegative(negativeVisible)), '-', ... % Convert x to a.
            'Color', stableColor, 'LineWidth', 1.8); % Mark this branch as attracting.
        hold on;                               % Keep the axes for the saddle-node pair.
        plot(alphaPositive(positiveStable), ... % Draw the stable branch born at the fold.
            0.5 * gamma .* (1 + xPositive(positiveStable)), '-', ... % Convert x to a.
            'Color', stableColor, 'LineWidth', 1.8); % Mark the branch as attracting.
        plot(alphaPositive(positiveUnstable), ... % Draw the middle branch approaching x=0.
            0.5 * gamma .* (1 + xPositive(positiveUnstable)), '--', ... % Convert x to a.
            'Color', unstableColor, 'LineWidth', 1.8); % Mark this branch as repelling.
        foldResidual = @(x) atanh(x) - x ./ (1 - x .^ 2) - hSmall; % Combine fold equations.
        xFold = fzero(foldResidual, [1e-12, 0.999999]); % Isolate the unique positive fold.
        alphaFold = 2 / (1 - xFold ^ 2);       % Recover the saddle-node alpha value.
        aFold = 0.5 * gamma * (1 + xFold);     % Convert the fold state to physical a.
        plot(alphaFold, aFold, 'ko', 'MarkerFaceColor', 'k', ... % Mark the saddle node.
            'MarkerSize', 4.5);                % Keep the marker clear but compact.
        title(sprintf('h = %.3f, \\alpha_{SN} = %.3f', hSmall, alphaFold), ... % State fold shift.
            'FontWeight', 'normal');           % Keep title weight consistent across panels.
        fprintf('Small-negative-bias alpha continuation: h=%.6g, alpha_SN=%.12g, x_SN=%.12g\n', ... % Report.
            hSmall, alphaFold, xFold);         % Print the current fold coordinates.
    end                                        % Finish the symmetric/imperfect branch choice.
    grid on;                                   % Add the same reading grid to every panel.
    xlabel('\alpha');                          % Keep alpha identified as the control parameter.
    ylabel('\ita');                            % Label the physical steady-state variable.
    xlim([alphaMin, alphaMax]);                % Apply one common alpha scale.
    ylim([0, gamma]);                          % Apply one common physical state scale.
end                                            % Finish the six-bias comparison loop.

sgtitle('Mirror unfolding at six small nonpositive biases (\delta = 0.5; \alpha is the control parameter)', ... % Main result.
    'FontName', 'Times New Roman', 'Interpreter', 'tex'); % Use the common typography.
kapitza_save_figure(gcf, 'small_negative_bias_alpha_bifurcation.png'); % Save PNG and vector PDF.
end                                            % Finish the negative-bias continuation.


function critical_point_trajectory_demo()     % Trace the fold point created by nonzero bias.
%CRITICAL_POINT_TRAJECTORY_DEMO Compare the exact cusp with its local scaling law.
delta = 0.5;                                   % Use delta; the teacher's short code calls it d.
gamma = 1 / (1 + delta);                       % Compute the physical interval length for a.
xFold = linspace(-0.45, 0.45, 2401);          % Parameterize both sides of the local cusp.
hExact = atanh(xFold) - xFold ./ (1 - xFold .^ 2); % Evaluate the exact fold bias.
alphaExact = 2 ./ (1 - xFold .^ 2);           % Recover alpha from the tangency condition.
aExact = 0.5 * gamma .* (1 + xFold);          % Convert the normalized fold state to physical a.
[hSorted, order] = sort(hExact);               % Order the parametric curve by increasing bias.
alphaSorted = alphaExact(order);               % Apply the same ordering to the alpha threshold.
aSorted = aExact(order);                       % Apply the same ordering to the physical state.
alphaAsymptotic = 2 + 2 .* (1.5 .* abs(hSorted)) .^ (2 / 3); % Evaluate the local cusp law.
hPositive = [0.002, 0.005, 0.010, 0.020, 0.040]; % Reuse the positive samples from Figure 5.
hSamples = [-fliplr(hPositive), 0, hPositive]; % Add their reflected negative-bias counterparts.
xSamples = zeros(size(hSamples));              % Allocate one exact fold state per sample.
alphaSamples = zeros(size(hSamples));          % Allocate one exact fold threshold per sample.
for j = 1:numel(hSamples)                      % Solve the fold equation at every sample bias.
    hValue = hSamples(j);                      % Read the current signed bias.
    if hValue > 0                              % A positive bias has a negative fold state.
        bracket = [-0.999999, -1e-12];         % Bracket that unique negative root.
    elseif hValue < 0                          % A negative bias has a positive fold state.
        bracket = [1e-12, 0.999999];           % Bracket the reflected positive root.
    else                                       % Handle the symmetric pitchfork limit directly.
        xSamples(j) = 0;                       % The limiting fold state is the center x=0.
        alphaSamples(j) = 2;                   % The limiting threshold is alpha=2.
        continue;                              % Skip the unnecessary numerical root solve.
    end                                        % Finish the sign-dependent bracket choice.
    foldResidual = @(x) atanh(x) - x ./ (1 - x .^ 2) - hValue; % Define the fold equation.
    xSamples(j) = fzero(foldResidual, bracket); % Isolate the exact saddle-node state.
    alphaSamples(j) = 2 / (1 - xSamples(j) ^ 2); % Recover the matching alpha threshold.
end                                            % Finish the sample-point loop.
aSamples = 0.5 * gamma .* (1 + xSamples);      % Convert all sample states to physical a.

figure('Color', 'w');                          % Open a white-background trajectory figure.
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact'); % Create two projections.
nexttile;                                      % Activate the parameter-space projection.
plot(hSorted, alphaSorted, 'b-', 'LineWidth', 2.2); % Draw the exact cusp boundary.
hold on;                                       % Keep the axes for approximation and samples.
plot(hSorted, alphaAsymptotic, 'r--', 'LineWidth', 1.8); % Draw the local scaling law.
plot(hSamples, alphaSamples, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4.5); % Mark samples.
plot(0, 2, 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 6); % Mark the pitchfork limit.
grid on;                                       % Add a light grid for quantitative comparison.
xlabel('\ith');                              % Label the fixed bias coordinate.
ylabel('\alpha_{SN}');                       % Label the moving multiplicity threshold.
title('Multiplicity threshold in parameter space'); % Describe the first trajectory projection.
legend('Exact fold trajectory', '2 + 2(3|h|/2)^{2/3}', ... % Identify exact and asymptotic curves.
    'Sampled biases', 'Pitchfork limit', 'Location', 'north'); % Complete the legend.

nexttile;                                      % Activate the physical-state projection.
plot(hSorted, aSorted, 'b-', 'LineWidth', 2.2); % Draw the exact physical fold state.
hold on;                                       % Keep the axes for samples and reference point.
plot(hSamples, aSamples, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4.5); % Mark samples.
plot(0, 0.5 * gamma, 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 6); % Mark center.
yline(0.5 * gamma, ':', 'Color', [0.45, 0.45, 0.45], 'LineWidth', 1); % Show symmetry level.
text(0.018, 0.285 * gamma, 'h > 0: lower saddle node', ... % Explain the positive-bias direction.
    'FontName', 'Times New Roman');            % Preserve the common text face.
text(-0.075, 0.725 * gamma, 'h < 0: upper saddle node', ... % Explain negative-bias direction.
    'FontName', 'Times New Roman');            % Preserve the common text face.
grid on;                                       % Add a light grid for reading the moving state.
xlabel('\ith');                              % Label the signed fixed bias coordinate.
ylabel('\ita_{SN}');                         % Label the moving physical critical state.
title('Critical state moves away from the center'); % Describe the second trajectory projection.
legend('Exact fold state', 'Sampled biases', 'Pitchfork limit', ... % Identify curves and markers.
    'Location', 'east');                      % Place the compact legend inside the panel.

sgtitle('Saddle-node trajectory from the pitchfork limit (\delta = 0.5)', ... % State origin.
    'FontName', 'Times New Roman', 'Interpreter', 'tex'); % Use the common typography.
kapitza_save_figure(gcf, 'critical_point_trajectory.png'); % Save PNG and vector PDF copies.
end                                            % Finish the critical-point trajectory demo.


function s_curve_demo()                        % Plot the imperfect pitchfork and flux roots.
%S_CURVE_DEMO Introduce a conductance-amplitude bias and compute both folds.
delta = 0.5;                                   % Use delta; the teacher's short code calls it d.
alpha = 3;                                     % Choose alpha>2 so two folds exist.
kappa = alpha / 2;                             % Convert alpha to the normalized slope.
xFold = sqrt(1 - 1 / kappa);                   % Evaluate the positive fold state.
hCritical = kappa * xFold - atanh(xFold);      % Evaluate the positive critical bias magnitude.
hSample = 0.5 * hCritical;                     % Select a bias inside the three-root region.
x = linspace(-0.995, 0.995, 2400);             % Parameterize the continuation curve.
hCurve = atanh(x) - kappa .* x;                % Evaluate the exact relation h(x).
xRoots = kapitza_roots(alpha, hSample);         % Find all states at the selected bias.
aRoots = (1 + xRoots) ./ (2 * (1 + delta));    % Convert roots to physical a values.
sRoots = (1 - xRoots) ./ 2;                    % Convert roots to jumps s=[u].
qRoots = delta .* aRoots;                      % Compute each steady heat flux.
figure('Color', 'w');                          % Open a white-background figure.
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact'); % Create two panels.
nexttile;                                      % Activate the continuation panel.
middle = abs(x) <= xFold;                      % Identify the repelling middle branch.
plot(hCurve(middle), x(middle), 'r--', 'LineWidth', 2); % Draw it as dashed.
hold on;                                       % Keep the axes for stable branches and points.
plot(hCurve(x < -xFold), x(x < -xFold), 'b-', 'LineWidth', 2); % Draw lower stable branch.
plot(hCurve(x > xFold), x(x > xFold), 'b-', 'LineWidth', 2); % Draw upper stable branch.
plot([hCritical, -hCritical], [-xFold, xFold], 'ko', 'MarkerFaceColor', 'k'); % Mark folds.
plot(hSample * ones(size(xRoots)), xRoots, 'mo', 'MarkerFaceColor', 'm'); % Mark sample roots.
grid on;                                       % Add a continuation grid.
xlabel('\ith');                                % Label the bias parameter in mathematical italics.
ylabel('\itx');                                % Label the normalized state in mathematical italics.
title(sprintf('S-curve: \\alpha = %.1f (\\kappa = %.2f)', alpha, kappa)); % Show both conventions.
legend('Unstable middle branch', 'Stable outer branches', '', 'Folds', ... % Name geometry.
    'Sample roots', 'Location', 'best');       % Complete the continuation legend.
xlim([-2.2 * hCritical, 2.2 * hCritical]);     % Match the focused Python continuation range.
nexttile;                                      % Activate the heat-flux panel.
s = linspace(0, 1, 2000);                     % Sample the complete jump interval.
qKapitza = kapitza_beta(s, delta, alpha, hSample) .* s; % Evaluate q_K=beta_h(s)s.
qBulk = delta / (1 + delta) .* (1 - s);       % Evaluate the bulk-required heat flux.
plot(s, qKapitza, 'r-', 'LineWidth', 2);       % Draw the nonlinear interface flux.
hold on;                                       % Keep the axes for the bulk flux and roots.
plot(s, qBulk, 'b-', 'LineWidth', 2);          % Draw the linear bulk flux.
plot(sRoots, qRoots, 'ko', 'MarkerFaceColor', 'k'); % Mark steady-state intersections.
grid on;                                       % Add a grid for reading flux values.
xlabel('\its = [u]');                          % Label the jump variable in mathematical italics.
ylabel('\itq');                                % Label the heat flux in mathematical italics.
title(sprintf('Three flux intersections: h = %.4f', hSample)); % Show the sample bias.
legend('\itq_K = \beta_h(s)s', '\itq_{bulk}', 'Roots', ... % Typeset flux formulas as mathematics.
    'Location', 'northeast');                  % Complete the heat-flux legend.
fprintf('\nS-curve: delta=%.6g, alpha=%.6g\n', delta, alpha); % Print selected parameters.
fprintf('x_f=%.12g, h_c=%.12g, h_sample=%.12g\n', xFold, hCritical, hSample); % Print folds.
fprintf('Number of sample roots: %d\n', numel(xRoots)); % Print the root count.
kapitza_save_figure(gcf, 's_curve.png');       % Save the S-curve figure.
end                                            % Finish the S-curve demo.


function s_curve_fixed_point_cobweb_demo()     % Show biased intersections and cobweb basins.
%S_CURVE_FIXED_POINT_COBWEB_DEMO Analyze one biased three-root state.
delta = 0.5;                                   % Use delta; the teacher's short code calls it d.
alpha = 3;                                     % Choose a supercritical alpha value.
kappa = alpha / 2;                             % Convert alpha to the normalized slope.
xFold = sqrt(1 - 1 / kappa);                   % Compute the positive fold state.
hCritical = kappa * xFold - atanh(xFold);      % Compute the critical bias magnitude.
hSample = 0.5 * hCritical;                     % Select a bias in the three-root region.
gamma = 1 / (1 + delta);                       % Compute the physical interval length.
a = linspace(0, gamma, 2000);                  % Sample the full physical interval.
F = kapitza_map(a, delta, alpha, hSample);     % Evaluate the biased fixed-point map.
xRoots = kapitza_roots(alpha, hSample);         % Find all three normalized roots.
aRoots = (1 + xRoots) ./ (2 * (1 + delta));    % Convert roots to physical states.
middleFixedPoint = aRoots(2);                   % Read the repelling root that separates both basins.
initialOffset = 0.08 * gamma;                   % Stay near the separator to reveal transient paths.
initialValues = middleFixedPoint + ...          % Place one initial condition in each basin.
    [-initialOffset, initialOffset];             % Start symmetrically around the repelling root.
pathColors = [0.10, 0.55, 0.30; 0.75, 0.20, 0.15]; % Assign colors to both paths.
figure('Color', 'w');                          % Open a white-background figure.
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact'); % Create two panels.
nexttile;                                      % Activate the intersection panel.
plot(a, F, 'b-', 'LineWidth', 2);              % Draw the biased nonlinear map.
hold on;                                       % Keep the axes for identity and roots.
plot(a, a, 'r--', 'LineWidth', 1.5);           % Draw the identity line.
plot(aRoots, aRoots, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6); % Mark roots.
grid on;                                       % Add a grid for reading intersections.
xlabel('\ita');                                % Label the physical state in mathematical italics.
ylabel('\itF_{\alpha,h}(a)');                 % Label the biased map in mathematical italics.
title('Three fixed-point intersections');      % State the panel conclusion.
legend('\itF_{\alpha,h}(a)', '\ita', 'Fixed points', ... % Typeset the map entries as mathematics.
    'Location', 'northwest');                  % Complete the intersection legend.
xlim([0, gamma]);                              % Restrict the horizontal interval.
ylim([0, gamma]);                              % Restrict the vertical interval.
nexttile;                                      % Activate the cobweb panel.
plot(a, F, 'b-', 'LineWidth', 1.8);            % Draw the same biased map.
hold on;                                       % Keep the axes for identity and paths.
plot(a, a, 'k--', 'LineWidth', 1.2);           % Draw the identity line.
for k = 1:numel(initialValues)                 % Draw paths from both outer basins.
    [pathX, pathY] = kapitza_cobweb_path(initialValues(k), delta, alpha, hSample, 35); % Build path.
    plot(pathX, pathY, '-', 'Color', pathColors(k, :), 'LineWidth', 0.9); % Draw path.
end                                            % Finish the two-basin loop.
plot(aRoots, aRoots, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5); % Mark roots.
grid on;                                       % Add a grid for reading iterations.
xlabel('\ita_n');                              % Label the current iterate in mathematical italics.
ylabel('\ita_{n+1}');                          % Label the next iterate in mathematical italics.
title('Two attracting basins');                % State the cobweb conclusion.
legend('\itF_{\alpha,h}(a)', '\ita', 'lower start', 'upper start', ... % Name curves and paths.
    'Fixed points', 'Location', 'northwest');  % Complete the legend.
xlim([0, gamma]);                              % Apply the physical horizontal range.
ylim([0, gamma]);                              % Apply the physical vertical range.
sgtitle(sprintf('Fixed-h slice: \\delta = %.1f, \\alpha = %.1f, h = %.4f', ... % Build title.
    delta, alpha, hSample), 'FontName', 'Times New Roman', ... % Insert parameters.
    'Interpreter', 'tex');                     % Match the rest of the figure typography.
kapitza_save_figure(gcf, 's_curve_fixed_point_cobweb.png'); % Save the completed figure.
end                                            % Finish the biased cobweb demo.


function fixed_a_h_iteration_demo()            % Recover h for one prescribed physical state a_star.
%FIXED_A_H_ITERATION_DEMO Show parameter inversion, not a physical bifurcation.
delta = 0.5;                                   % Use the same material parameter as all examples.
alpha = 3;                                     % Keep the supercritical sigmoid steepness fixed.
kappa = alpha / 2;                             % Convert alpha to the normalized slope.
xFold = sqrt(1 - 1 / kappa);                   % Compute the positive continuation fold state.
hCritical = kappa * xFold - atanh(xFold);      % Compute the positive fold-bias magnitude.
hStar = 0.5 * hCritical;                       % Reuse the sample bias from the S-curve figures.
xRoots = kapitza_roots(alpha, hStar);           % Find the three states at the sample bias.
xStar = xRoots(end);                           % Prescribe the upper attracting normalized state.
aStar = (1 + xStar) / (2 * (1 + delta));       % Convert the target state to the physical variable a.
h = linspace(-2.2 * hCritical, 2.2 * hCritical, 2000); % Use the focused continuation interval.
residual = tanh(kappa * xStar + h) - xStar;    % Evaluate the fixed-a residual R_{a_*}(h).
hMap = h + xStar - tanh(kappa * xStar + h);    % Define one relaxed parameter-recovery map.
initialValues = [-1.6, 1.8] .* hCritical;      % Start the h iteration on both sides of hStar.
pathColors = [0.10, 0.55, 0.30; 0.75, 0.20, 0.15]; % Assign colors to both starting guesses.
figure('Color', 'w');                          % Open a white-background comparison figure.
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact'); % Create two panels.
nexttile;                                      % Activate the scalar residual panel.
plot(h, residual, 'b-', 'LineWidth', 2);       % Draw the condition for the prescribed state.
hold on;                                       % Keep the axes for the zero line and recovered bias.
yline(0, 'r--', 'LineWidth', 1.5);             % Draw the residual root condition.
plot(hStar, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6); % Mark the unique recovered bias.
grid on;                                       % Add a grid for reading the root.
xlabel('\ith');                                % Label the bias variable in mathematical italics.
ylabel('\itR_{a_*}(h)');                       % Label the fixed-a residual in mathematical italics.
title('A fixed state determines one bias');    % State the parameter-inversion conclusion.
legend('\itR_{a_*}(h)', '\itR_{a_*}=0', '\ith_*', ... % Name the residual geometry.
    'Location', 'northwest');                  % Complete the residual legend.
xlim([-2.2 * hCritical, 2.2 * hCritical]);     % Apply the focused h interval.
nexttile;                                      % Activate the relaxed-iteration panel.
plot(h, hMap, 'b-', 'LineWidth', 1.8);         % Draw the parameter-recovery map.
hold on;                                       % Keep the axes for identity, paths, and fixed bias.
plot(h, h, 'k--', 'LineWidth', 1.2);           % Draw the identity line for h.
for k = 1:numel(initialValues)                 % Build one cobweb from each side of hStar.
    pathX = zeros(1, 71);                      % Allocate 35 vertical-horizontal cobweb steps.
    pathY = zeros(1, 71);                      % Allocate the matching ordinate coordinates.
    hCurrent = initialValues(k);               % Read the current initial bias guess.
    pathX(1) = hCurrent;                       % Start the path on the identity line.
    pathY(1) = hCurrent;                       % Use the same initial ordinate.
    for m = 1:35                               % Generate the relaxed parameter-recovery iteration.
        hNext = hCurrent + xStar - tanh(kappa * xStar + hCurrent); % Update the bias guess.
        pathX(2 * m) = hCurrent;               % Add the vertical segment abscissa.
        pathY(2 * m) = hNext;                  % Move vertically to the recovery map.
        pathX(2 * m + 1) = hNext;              % Move horizontally to the identity line.
        pathY(2 * m + 1) = hNext;              % Complete this cobweb step.
        hCurrent = hNext;                      % Advance to the next bias iterate.
    end                                        % Finish the parameter-recovery loop.
    plot(pathX, pathY, '-', 'Color', pathColors(k, :), 'LineWidth', 0.9); % Draw this path.
end                                            % Finish the two-start loop.
plot(hStar, hStar, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5); % Mark the recovered fixed bias.
grid on;                                       % Add a grid for reading the h iteration.
xlabel('\ith_n');                              % Label the current bias iterate.
ylabel('\ith_{n+1}');                          % Label the next bias iterate.
title('Relaxed iteration recovers the bias');  % Explain the numerical role of this cobweb.
legend('\itH_{a_*}(h)', '\ith', 'lower start', 'upper start', ... % Name curves and paths.
    '\ith_*', 'Location', 'northwest');        % Complete the parameter-recovery legend.
xlim([-2.2 * hCritical, 2.2 * hCritical]);     % Apply the focused horizontal interval.
ylim([-2.2 * hCritical, 2.2 * hCritical]);     % Apply the same vertical interval.
sgtitle(sprintf('Fixed-a parameter recovery: \\delta = %.1f, \\alpha = %.1f, a_* = %.6f, h_* = %.6f', ... % Build title.
    delta, alpha, aStar, hStar), 'FontName', 'Times New Roman', ... % Insert selected values.
    'Interpreter', 'tex');                     % Render symbols and indices consistently.
kapitza_save_figure(gcf, 'fixed_a_h_iteration.png'); % Save PNG and vector-PDF copies.
end                                            % Finish the fixed-a parameter-recovery demo.


function kapitza_save_figure(fig, filename)    % Save a figure in output/matlab/.
%KAPITZA_SAVE_FIGURE Locate the project root from this consolidated script.
scriptFile = which('kapitza_numerical');       % Obtain the absolute path of this script.
matlabDir = fileparts(scriptFile);             % Extract the code/matlab directory.
projectRoot = fileparts(fileparts(matlabDir)); % Move upward twice to the project root.
outputDir = fullfile(projectRoot, 'output', 'matlab'); % Construct the MATLAB output folder.
if ~exist(outputDir, 'dir')                    % Test whether the output folder exists.
    mkdir(outputDir);                          % Create the missing output folder.
end                                            % Finish the folder check.
target = fullfile(outputDir, filename);        % Construct the complete PNG output filename.
[~, baseName, ~] = fileparts(filename);        % Remove the PNG extension for the PDF copy.
pdfTarget = fullfile(outputDir, [baseName, '.pdf']); % Construct the matching PDF output filename.
fontObjects = findall(fig, '-property', 'FontName'); % Include axes, legends, and layout titles.
for j = 1:numel(fontObjects)                   % Visit every object that exposes a font face.
    fontObjects(j).FontName = 'Times New Roman'; % Apply Times New Roman without exceptions.
end                                            % Finish the all-object font loop.
interpreterObjects = findall(fig, '-property', 'Interpreter'); % Find all markup-aware objects.
for j = 1:numel(interpreterObjects)            % Visit titles, labels, annotations, and legends.
    interpreterObjects(j).Interpreter = 'tex'; % Preserve Greek letters, subscripts, and superscripts.
end                                            % Finish the all-object interpreter loop.
axesObjects = findall(fig, 'Type', 'axes');    % Find every ordinary axes object in the figure.
for j = 1:numel(axesObjects)                   % Visit each axes before exporting the figure.
    axesObjects(j).Toolbar.Visible = 'off';    % Hide interactive axes controls in the saved image.
    axesObjects(j).FontName = 'Times New Roman'; % Use Times New Roman for every tick label.
    axesObjects(j).TickLabelInterpreter = 'tex'; % Interpret tick-label markup consistently.
end                                            % Finish the axes-toolbar cleanup loop.
textObjects = findall(fig, 'Type', 'text');    % Find titles and all x/y-axis labels.
for j = 1:numel(textObjects)                   % Visit every textual annotation in the figure.
    textObjects(j).FontName = 'Times New Roman'; % Apply the requested Times New Roman face.
    textObjects(j).Interpreter = 'tex';        % Retain subscripts and Greek mathematical symbols.
end                                            % Finish the text-formatting loop.
legendObjects = findall(fig, 'Type', 'legend'); % Find every legend in the current figure.
for j = 1:numel(legendObjects)                 % Visit all legends before either export.
    legendObjects(j).FontName = 'Times New Roman'; % Apply Times New Roman to legend entries.
    legendObjects(j).Interpreter = 'tex';      % Interpret formulas and subscripts in the legend.
end                                            % Finish the legend-formatting loop.
originalUnits = fig.Units;                     % Remember the figure's on-screen unit system.
fig.Units = 'inches';                          % Express the current figure size in physical units.
figurePosition = fig.Position;                 % Read the width and height in inches.
fig.PaperUnits = 'inches';                     % Use the same units for the PDF page.
fig.PaperPosition = [0, 0, figurePosition(3), figurePosition(4)]; % Fill the PDF page.
fig.PaperSize = figurePosition(3:4);            % Match the PDF page to the figure aspect ratio.
print(fig, target, '-dpng', '-r220', '-painters'); % Use the TeX-aware print path for the PNG.
print(fig, pdfTarget, '-dpdf', '-painters');    % Save a tightly sized vector PDF.
fig.Units = originalUnits;                     % Restore the original on-screen unit system.
end                                            % Finish the output helper.
