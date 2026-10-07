%% Plot lambda sweep results from simulation_choose_params_lambda_EMD.m
%
% Loads FlexMF_choose_lambda_<data_type>_lambdaM=<lambda_M>.mat and creates
% a 4-panel figure (vs. lambda, log x-axis):
%   1. Normalized regularization cost (reg_costs) and L1 norm of R
%      (transport residual), with the per-sim crossing
%   2. EMD of W
%   3. EMD of H
%   4. Number of significant factors
% and, for a selected lambda/simulation, a ground-truth-vs-reconstruction
% figure over a plot_range covering ~5 occurrences of each sequence.
%
% Configure data_type and lambda_M below, then run this script from the
% FlexMF directory.

clear all
close all
clc

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir)
    this_dir = pwd;
end
addpath(genpath(this_dir))

%% Configuration
data_type = 'warpnoise';
lambda_M = 1e-1;
output_dir = 'Simulation_EMD';
save_figure = true;
selected_lambda_ids = [5,8,9,12];  % [] selects the lambda with the lowest EMD of W
sim_idx = 1;  % simulation index used for the reconstruction figure
n_occurrences = 5;  % number of sequence occurrences to show in the reconstruction figure

results_file = fullfile(this_dir, output_dir, ...
    sprintf('FlexMF_choose_lambda_%s_lambdaM=%0.3e.mat', data_type, lambda_M));
assert(exist(results_file, 'file') == 2, ...
    'Results file not found: %s', results_file)

S = load(results_file);
required_fields = {'lambdas', 'reg_costs', 'Rs', 'emds_W', 'emds_H', 'num_significant', ...
    'W_hats', 'H_hats', 'Ws', 'Hs', 'Xs', 'num_detected', 'ids_match'};
for i = 1:numel(required_fields)
    assert(isfield(S, required_fields{i}), ...
        'Missing variable ''%s'' in %s', required_fields{i}, results_file)
end

lambdas = S.lambdas(:);
reg_costs = S.reg_costs;
L1_R = cellfun(@(r) norm(r(:), 1), S.Rs);
emds_W = cellfun(@mean_emd, S.emds_W);
emds_H = cellfun(@mean_emd, S.emds_H);
num_significant = cell2mat(S.num_significant);
nSim = size(reg_costs, 2);

%% Crossing of normalized reg cost and ||R||_1
reg_norm = normalize_prctile(reg_costs);
R_norm = normalize_prctile(L1_R);

loglam = log10(lambdas);
loglam_cross = nan(1, nSim);
for ii = 1:nSim
    loglam_cross(ii) = find_crossing(loglam, reg_norm(:, ii) - R_norm(:, ii));
end
lambda_cross = 10.^loglam_cross;
lambda_cross_med = 10^median(loglam_cross, 'omitnan');

fprintf('\nCrossing of normalized reg cost and ||R||_1 per simulation:\n');
for ii = 1:nSim
    fprintf('  sim %2d: lambda = %.3e\n', ii, lambda_cross(ii));
end
fprintf('  median (in log lambda): %.3e  (%d/%d sims cross)\n', ...
    lambda_cross_med, sum(~isnan(loglam_cross)), nSim);

true_K = 3;  % ground-truth number of motifs (K in simulation_choose_params_lambda_EMD.m)

%% Figure
reg_color = [0.85 0.30 0.20];
r_color = [0.20 0.45 0.85];
k_color = [0.15 0.60 0.35];
figure('Color', 'w', 'Position', [80 40 500 900]);
layout = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(layout, 1);
hold(ax1, 'on')
h1 = errorbar_stats(ax1, lambdas, reg_norm, reg_color);
h2 = errorbar_stats(ax1, lambdas, R_norm, r_color);
if ~isnan(lambda_cross_med)
    xline(ax1, lambda_cross_med, '--k');
end
hold(ax1, 'off')
set(ax1, 'XScale', 'log', 'XTickLabel', [])
box(ax1, 'off')
ylabel(ax1, 'Normalized cost')
legend(ax1, [h1, h2], {'Reg cost', '$\|R\|_1$'}, 'Location', 'best', ...
    'Box', 'off', 'Interpreter', 'latex')
title(ax1, sprintf(['%s: FlexMF lambda sweep (lambda_M=%0.3e)\n', ...
    'median crossing lambda = %.3e'], data_type, lambda_M, lambda_cross_med), ...
    'Interpreter', 'none')

ax2 = nexttile(layout, 2);
errorbar_stats(ax2, lambdas, emds_W, 'k');
set(ax2, 'XScale', 'log', 'XTickLabel', [])
box(ax2, 'off')
ylabel(ax2, 'EMD of W')

ax3 = nexttile(layout, 3);
errorbar_stats(ax3, lambdas, emds_H, 'k');
set(ax3, 'XScale', 'log', 'XTickLabel', [])
box(ax3, 'off')
ylabel(ax3, 'EMD of H')

ax4 = nexttile(layout, 4);
hold(ax4, 'on')
errorbar_stats(ax4, lambdas, num_significant, 'k');
yline(ax4, true_K, '--', 'Color', k_color, 'Label', sprintf('true K=%d', true_K), ...
    'LabelHorizontalAlignment', 'left');
hold(ax4, 'off')
set(ax4, 'XScale', 'log')
box(ax4, 'off')
ylabel(ax4, 'Number significant')
ylim(ax4, [0, max(true_K, max(num_significant(:))) + 0.5])

xlabel(layout, '\lambda')
linkaxes([ax1, ax2, ax3, ax4], 'x')

if save_figure
    if ~exist(fullfile(this_dir, output_dir), 'dir')
        mkdir(fullfile(this_dir, output_dir))
    end
    output_file = fullfile(this_dir, output_dir, ...
        sprintf('Choose_lambda_%s_lambdaM=%0.3e_plot.pdf', data_type, lambda_M));
    export_vector_pdf(output_file, gcf)
    fprintf('Saved figure: %s\n', output_file)
end

%% Ground truth vs. reconstruction at a selected lambda
if isempty(selected_lambda_ids)
    [~, selected_lambda_ids] = min(emds_W(:, sim_idx));
end
for li = 1:numel(selected_lambda_ids)
    level = selected_lambda_ids(li);
    lambda = lambdas(selected_lambda_ids(li)); 

    W_true = S.Ws{selected_lambda_ids(li), sim_idx};
    H_true = S.Hs{selected_lambda_ids(li), sim_idx};
    X_true = S.Xs{selected_lambda_ids(li), sim_idx};
    W_hat = S.W_hats{selected_lambda_ids(li), sim_idx};
    H_hat = S.H_hats{selected_lambda_ids(li), sim_idx};
    is_sig = S.is_significant{selected_lambda_ids(li), sim_idx};

    % Reorder estimated factors to align with ground-truth motif order
    [W_hat_plot, H_hat_plot, order] = helper.sort_matched_factors(W_hat, H_hat, S.ids_match{selected_lambda_ids(li), sim_idx});
    is_sig_plot = is_sig(order);

    T = size(H_true, 2);
    Htrain = H_true(:, 1:round(T/2));
    Xtrain = X_true(:, 1:round(T/2));

    plot_range = helper.estimate_plot_range(Htrain, n_occurrences);
    lambda_label = sprintf('%s: lambda=%0.3e, lambda_M=%0.3e', data_type, lambda, lambda_M);

    figure('Color', 'w', 'Name', sprintf('%s - ground truth', lambda_label));
    SimpleWHPlot_patch(W_true, Htrain, 'Data', Xtrain, 'plot_range', plot_range, ...
        'is_significant', ones(1, size(W_true, 2)), 'center', true);
    sgtitle(sprintf('Ground truth (%s)', lambda_label), 'Interpreter', 'none')
    save_lambda_figure(save_figure, this_dir, output_dir, data_type, lambda_M, selected_lambda_ids(li), 'groundtruth')

    figure('Color', 'w', 'Name', sprintf('%s - FlexMF', lambda_label));
    SimpleWHPlot_patch(W_hat_plot, H_hat_plot, 'Data', Xtrain, 'plot_range', plot_range, ...
        'is_significant', is_sig_plot, 'center', true);
    sgtitle(sprintf('FlexMF reconstruction (%s, %d/%d significant)', lambda_label, ...
        num_significant(selected_lambda_ids(li), sim_idx), size(W_hat, 2)), 'Interpreter', 'none')
    save_lambda_figure(save_figure, this_dir, output_dir, data_type, lambda_M, selected_lambda_ids(li), 'FlexMF')
end

%% Local functions
function value = mean_emd(x)
if isempty(x)
    value = NaN;
else
    value = mean(x(:), 'omitnan');
end
end

function save_lambda_figure(save_figure, this_dir, output_dir, data_type, lambda_M, lambda_idx, tag)
if ~save_figure
    return
end
if ~exist(fullfile(this_dir, output_dir), 'dir')
    mkdir(fullfile(this_dir, output_dir))
end
output_file = fullfile(this_dir, output_dir, ...
    sprintf('Choose_lambda_%s_lambdaM=%0.3e_idx%d_%s.pdf', data_type, lambda_M, lambda_idx, tag));
exportgraphics(gcf, output_file, 'ContentType', 'vector')
fprintf('Saved figure: %s\n', output_file)
end

function y = normalize_prctile(x)
% Rescale so the 10th / 90th percentiles over the whole grid map to 0 / 1
lo = prctile(x(:), 10);
hi = prctile(x(:), 90);
y = (x - lo) / max(hi - lo, eps);
end

function x0 = find_crossing(x, d)
% First sign change of d(x), linearly interpolated in x; NaN if none
x0 = NaN;
for k = 1:numel(d) - 1
    if d(k) == 0
        x0 = x(k);
        return
    end
    if d(k) * d(k + 1) < 0
        x0 = x(k) + d(k) / (d(k) - d(k + 1)) * (x(k + 1) - x(k));
        return
    end
end
end

function h = errorbar_stats(ax, x, values, color)
% values: [nLambdas x nSim] matrix; plots median with IQR error bars
med = median(values, 2, 'omitnan');
lo = prctile(values, 25, 2);
hi = prctile(values, 75, 2);
h = errorbar(ax, x, med, med - lo, hi - med, ...
    '-', 'Marker', '.', 'MarkerSize', 12, 'Color', color);
end
