%% Plot lambda sweep results from simulation_choose_params_lambda_EMD.m
%
% Loads FlexMF_choose_lambda_<data_type>_lambdaM=<lambda_M>.mat and creates
% a 4-panel figure (vs. lambda, log x-axis):
%   1. Regularization cost (reg_costs) and L1 norm of R (transport residual)
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
data_type = 'jitternoise';
lambda_M = 1e-1;
output_dir = 'Simulation_EMD';
save_figure = true;
selected_lambda_ids = [6,7,8,9,10];  % [] selects the lambda with the lowest EMD of W
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

lambdas = S.lambdas;
reg_costs = S.reg_costs;
L1_R = cellfun(@(r) norm(r(:), 1), S.Rs);
emds_W = cellfun(@mean_emd, S.emds_W);
emds_H = cellfun(@mean_emd, S.emds_H);
num_significant = cell2mat(S.num_significant);

%% Figure
reg_color = [0.85 0.30 0.20];
r_color = [0.20 0.45 0.85];
figure('Color', 'w', 'Position', [80 40 500 900]);
layout = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(layout, 1);
hold(ax1, 'on')
yyaxis(ax1, 'left')
errorbar_stats(ax1, lambdas, reg_costs, reg_color)
ylabel(ax1, 'Reg cost')
ax1.YColor = reg_color;
yyaxis(ax1, 'right')
errorbar_stats(ax1, lambdas, L1_R, r_color)
ylabel(ax1, '||R||_1')
ax1.YColor = r_color;
set(ax1, 'XScale', 'log', 'XTickLabel', [])
title(ax1, sprintf('%s: FlexMF lambda sweep (lambda_M=%0.3e)', data_type, lambda_M), ...
    'Interpreter', 'none')

ax2 = nexttile(layout, 2);
errorbar_stats(ax2, lambdas, emds_W, 'k')
set(ax2, 'XScale', 'log', 'XTickLabel', [])
ylabel(ax2, 'EMD of W')

ax3 = nexttile(layout, 3);
errorbar_stats(ax3, lambdas, emds_H, 'k')
set(ax3, 'XScale', 'log', 'XTickLabel', [])
ylabel(ax3, 'EMD of H')

ax4 = nexttile(layout, 4);
errorbar_stats(ax4, lambdas, num_significant, 'k')
set(ax4, 'XScale', 'log')
ylabel(ax4, 'Number significant')

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
        'is_significant', ones(1, size(W_true, 2)));
    sgtitle(sprintf('Ground truth (%s)', lambda_label), 'Interpreter', 'none')
    save_lambda_figure(save_figure, this_dir, output_dir, data_type, lambda_M, selected_lambda_ids(li), 'groundtruth')

    figure('Color', 'w', 'Name', sprintf('%s - FlexMF', lambda_label));
    SimpleWHPlot_patch(W_hat_plot, H_hat_plot, 'Data', Xtrain, 'plot_range', plot_range, ...
        'is_significant', is_sig_plot);
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

function errorbar_stats(ax, x, values, color)
% values: [nLambdas x nSim] matrix; plots median with IQR error bars
med = median(values, 2, 'omitnan');
lo = prctile(values, 25, 2);
hi = prctile(values, 75, 2);
errorbar(ax, x, med, med - lo, hi - med, ...
    '-', 'Marker', '.', 'MarkerSize', 12, 'Color', color)
end
