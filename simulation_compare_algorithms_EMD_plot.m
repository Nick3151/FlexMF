%% Plot SeqNMF and FlexMF EMD comparison results
%
% Loads Compare_FlexMF_SeqNMF_<data_type>.mat and creates:
%   - A summary figure with EMD of W, EMD of H, and number of significant
%     factors across simulations
%   - For each selected level, separate figures showing the ground truth,
%     SeqNMF reconstruction, and FlexMF reconstruction over the same
%     plot_range (covering ~5 occurrences of each sequence)
%
% Configure data_type and selected_levels below, then run this script from
% the FlexMF directory.

clear all
close all
clc

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir)
    this_dir = pwd;
end
addpath(genpath(this_dir))

%% Configuration
data_type = 'jitter';
selected_levels = 5;  % [] selects the last level
output_dir = 'Simulation_EMD';
save_figure = true;
true_K = 3;  % ground-truth number of motifs (K in simulation_compare_algorithms_EMD.m)

results_file = fullfile(this_dir, output_dir, ...
    sprintf('Compare_FlexMF_SeqNMF_%s.mat', data_type));
assert(exist(results_file, 'file') == 2, ...
    'Results file not found: %s', results_file)

S = load(results_file);
required_fields = {'emds_W_SeqNMF', 'emds_W_FlexMF', ...
    'emds_H_SeqNMF', 'emds_H_FlexMF', ...
    'is_significant_SeqNMF', 'is_significant_FlexMF', ...
    'Whats_SeqNMF', 'Hhats_train_SeqNMF', ...
    'Whats_FlexMF', 'Hhats_train_FlexMF', ...
    'ids_SeqNMF', 'ids_FlexMF', ...
    'Ws', 'Hs_train', 'Xs_train'};
for i = 1:numel(required_fields)
    assert(isfield(S, required_fields{i}), ...
        'Missing variable ''%s'' in %s', required_fields{i}, results_file)
end

emds_W_seq = cellfun(@mean_emd, S.emds_W_SeqNMF);
emds_W_flex = cellfun(@mean_emd, S.emds_W_FlexMF);
emds_H_seq = cellfun(@mean_emd, S.emds_H_SeqNMF);
emds_H_flex = cellfun(@mean_emd, S.emds_H_FlexMF);
n_sig_seq = cellfun(@count_significant, S.is_significant_SeqNMF);
n_sig_flex = cellfun(@count_significant, S.is_significant_FlexMF);

% EMD summaries only reflect runs that recovered exactly true_K significant motifs
emds_W_seq(n_sig_seq ~= true_K) = NaN;
emds_H_seq(n_sig_seq ~= true_K) = NaN;
emds_W_flex(n_sig_flex ~= true_K) = NaN;
emds_H_flex(n_sig_flex ~= true_K) = NaN;

[n_levels, n_sim] = size(emds_W_seq);
assert(isequal(size(emds_W_flex), [n_levels, n_sim]), ...
    'SeqNMF and FlexMF W results have different dimensions')

level_values = comparison_levels(data_type, n_levels);
if isempty(selected_levels)
    selected_levels = n_levels;
end
selected_levels = selected_levels(:)';
assert(all(selected_levels >= 1 & selected_levels <= n_levels) && ...
    all(selected_levels == floor(selected_levels)), ...
    'selected_levels must contain valid level indices')

%% Summary figure
seq_color = [0.20 0.45 0.85];
flex_color = [0.85 0.30 0.20];
offset = 0.16;
figure('Color', 'w', 'Position', [80 40 1000 700]);
layout = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax_w = nexttile(layout, 1);
plot_summary(ax_w, emds_W_seq, emds_W_flex, seq_color, flex_color, ...
    'EMD of W', false, offset, level_values)

title(ax_w, sprintf('%s: SeqNMF vs FlexMF (%d sims/level; EMD from runs with exactly %d significant motifs)', ...
    data_type, n_sim, true_K), 'Interpreter', 'none')

ax_h = nexttile(layout, 2);
plot_summary(ax_h, emds_H_seq, emds_H_flex, seq_color, flex_color, ...
    'EMD of H', false, offset, level_values)

ax_sig = nexttile(layout, 3);
plot_summary(ax_sig, n_sig_seq, n_sig_flex, seq_color, flex_color, ...
    'Number of significant factors', true, offset, level_values)
yline(ax_sig, true_K, 'k--', 'true K', 'LabelHorizontalAlignment', 'left')

xlabel(layout, level_axis_label(data_type), 'Interpreter', 'none')

if save_figure
    if ~exist(fullfile(this_dir, output_dir), 'dir')
        mkdir(fullfile(this_dir, output_dir))
    end
    output_file = fullfile(this_dir, output_dir, ...
        sprintf('Compare_FlexMF_SeqNMF_%s_summary.pdf', data_type));
    export_vector_pdf(output_file, gcf)
    fprintf('Saved figure: %s\n', output_file)
end

%% Example reconstructions: ground truth vs. SeqNMF vs. FlexMF
n_occurrences = 5; % number of sequence occurrences to show in each example plot
sim_idx = 1; % simulation index used for example plots

for li = 1:numel(selected_levels)
    level = selected_levels(li);
    W_true = S.Ws{level, sim_idx};
    H_true = S.Hs_train{level, sim_idx};
    X_true = S.Xs_train{level, sim_idx};
    What_seq = S.Whats_SeqNMF{level, sim_idx};
    Hhat_seq = S.Hhats_train_SeqNMF{level, sim_idx};
    is_sig_seq = S.is_significant_SeqNMF{level, sim_idx};
    What_flex = S.Whats_FlexMF{level, sim_idx};
    Hhat_flex = S.Hhats_train_FlexMF{level, sim_idx};
    is_sig_flex = S.is_significant_FlexMF{level, sim_idx};

    if isempty(W_true) || isempty(What_seq) || isempty(What_flex)
        continue
    end

    % Reorder estimated factors to align with ground-truth motif order
    [What_seq_plot, Hhat_seq_plot, order_seq] = helper.sort_matched_factors(What_seq, Hhat_seq, S.ids_SeqNMF{level, sim_idx});
    is_sig_seq_plot = is_sig_seq(order_seq);
    [What_flex_plot, Hhat_flex_plot, order_flex] = helper.sort_matched_factors(What_flex, Hhat_flex, S.ids_FlexMF{level, sim_idx});
    is_sig_flex_plot = is_sig_flex(order_flex);

    plot_range = helper.estimate_plot_range(H_true, n_occurrences);
    level_label = level_title_label(data_type, level_values(level));

    figure('Color', 'w', 'Name', sprintf('%s - ground truth', level_label));
    SimpleWHPlot_patch(W_true, H_true, 'Data', X_true, 'plot_range', plot_range, ...
        'is_significant', ones(1, size(W_true, 2)));
    sgtitle(sprintf('Ground truth (%s)', level_label), 'Interpreter', 'none')
    save_example_figure(save_figure, this_dir, output_dir, data_type, level, 'groundtruth')

    figure('Color', 'w', 'Name', sprintf('%s - SeqNMF', level_label));
    SimpleWHPlot_patch(What_seq_plot, Hhat_seq_plot, 'plot_range', plot_range, ...
        'is_significant', is_sig_seq_plot);
    sgtitle(sprintf('SeqNMF reconstruction (%s)', level_label), 'Interpreter', 'none')
    save_example_figure(save_figure, this_dir, output_dir, data_type, level, 'SeqNMF')

    figure('Color', 'w', 'Name', sprintf('%s - FlexMF', level_label));
    SimpleWHPlot_patch(What_flex_plot, Hhat_flex_plot, 'plot_range', plot_range, ...
        'is_significant', is_sig_flex_plot);
    sgtitle(sprintf('FlexMF reconstruction (%s)', level_label), 'Interpreter', 'none')
    save_example_figure(save_figure, this_dir, output_dir, data_type, level, 'FlexMF')
end

%% Local functions
function value = mean_emd(x)
if isempty(x)
    value = NaN;
else
    value = mean(x(:), 'omitnan');
end
end

function value = count_significant(x)
if isempty(x)
    value = NaN;
else
    value = sum(logical(x(:)));
end
end

function plot_summary(ax, seq_values, flex_values, seq_color, flex_color, ...
    y_label, show_x_labels, offset, level_values)
axes(ax)
hold(ax, 'on')
[n_levels, n_sim] = size(seq_values);
x = 1:n_levels;

for level = 1:n_levels
    seq = seq_values(level, :);
    flex = flex_values(level, :);
    seq_stats = quartile_stats(seq);
    flex_stats = quartile_stats(flex);
    errorbar(ax, level - offset, seq_stats(1), ...
        seq_stats(1) - seq_stats(2), seq_stats(3) - seq_stats(1), ...
        'o', 'Color', seq_color, 'MarkerFaceColor', seq_color, ...
        'CapSize', 5, 'LineWidth', 1)
    errorbar(ax, level + offset, flex_stats(1), ...
        flex_stats(1) - flex_stats(2), flex_stats(3) - flex_stats(1), ...
        'o', 'Color', flex_color, 'MarkerFaceColor', flex_color, ...
        'CapSize', 5, 'LineWidth', 1)
end

xlim(ax, [0.5 n_levels + 0.5])
ylabel(ax, y_label)
box(ax, 'off')
set(ax, 'XTick', x)
if show_x_labels
    set(ax, 'XTickLabel', string(level_values))
else
    set(ax, 'XTickLabel', [])
end
if n_sim > 1
    legend(ax, {'SeqNMF', 'FlexMF'}, 'Location', 'best', 'Box', 'off')
end
hold(ax, 'off')
end

function stats = quartile_stats(values)
values = values(isfinite(values));
if isempty(values)
    stats = [NaN NaN NaN];
else
    stats = [median(values), prctile(values, 25), prctile(values, 75)];
end
end

function save_example_figure(save_figure, this_dir, output_dir, data_type, level, tag)
if ~save_figure
    return
end
if ~exist(fullfile(this_dir, output_dir), 'dir')
    mkdir(fullfile(this_dir, output_dir))
end
output_file = fullfile(this_dir, output_dir, ...
    sprintf('Compare_FlexMF_SeqNMF_%s_level%d_%s.pdf', data_type, level, tag));
exportgraphics(gcf, output_file, 'ContentType', 'vector')
fprintf('Saved figure: %s\n', output_file)
end

function label = level_title_label(data_type, level_value)
% Fixed noise level used for 'warpnoise'/'jitternoise' (noise_levels(3) in
% simulation_compare_algorithms_EMD.m)
fixed_noise = 0.002;
switch lower(data_type)
    case 'warp'
        label = sprintf('warp=%g, noise=0', level_value);
    case 'warpnoise'
        label = sprintf('warp=%g, noise=%g', level_value, fixed_noise);
    case 'jitter'
        label = sprintf('jitter=%g, noise=0', level_value);
    case 'jitternoise'
        label = sprintf('jitter=%g, noise=%g', level_value, fixed_noise);
    case 'participation'
        label = sprintf('participation=%g, noise=0', level_value);
    case 'noise'
        label = sprintf('noise=%g', level_value);
    otherwise
        label = sprintf('level=%g', level_value);
end
end

function levels = comparison_levels(data_type, n_levels)
switch lower(data_type)
    case {'noise'}
        levels = 0:0.001:0.01;
    case {'warp', 'warpnoise'}
        levels = 0:9;
    case {'jitter', 'jitternoise'}
        levels = 0:9;
    case 'participation'
        levels = 1:-0.1:0.1;
    otherwise
        levels = 1:n_levels;
end
levels = levels(1:min(n_levels, numel(levels)));
if numel(levels) < n_levels
    levels = 1:n_levels;
end
end

function label = level_axis_label(data_type)
switch lower(data_type)
    case 'noise'
        label = 'Additive noise level';
    case 'warp'
        label = 'Warp level';
    case 'warpnoise'
        label = 'Warp level (fixed noise)';
    case 'jitter'
        label = 'Jitter level';
    case 'jitternoise'
        label = 'Jitter level (fixed noise)';
    case 'participation'
        label = 'Participation Rate';
    otherwise
        label = 'Level';
end
end
