%% Plot results from compare_init_robustness
%
% Loads Simulation_Results/compare_init_robustness_<dataType>.mat files
% produced by compare_init_robustness(dataType), and writes figures to
% Simulation_Results/.
%
% Configure which data types to plot via dataTypes below. Missing mats are
% skipped with a warning so a partial job array can still be visualized.

clear all
close all
clc

thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir)
    thisDir = pwd;
end
root = fileparts(thisDir);
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))
addpath(genpath(thisDir));

%% Configuration
outDir = 'Simulation_Results';
dataTypes = {'clean', 'noise', 'warp', 'jitter', 'participation'};
% dataTypes = {'noise'};   % plot a subset

constraintTolDefault = 0.05;
plotAll = 1;

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% Load available results
results = struct('name', {}, 'label', {}, 'emds_W', {}, 'emds_H', {}, ...
    'nDetected', {}, 'constraints_rel', {}, 'objs', {}, 'nIters', {}, ...
    'time_flex', {}, 'time_match', {}, 'time_match_runs', {}, 'time_total', {}, ...
    'info_warm', {}, 'info_rand', {});
data = struct('name', {}, 'label', {}, 'X', {}, 'Wtrue', {}, 'Htrue', {}, 'L', {});
methodLabels = {};
nSim = [];
K = [];
Khat = [];
maxiter = [];
constraintTol = constraintTolDefault;
loaded = cell(numel(dataTypes), 1);
nLoaded = 0;

for i = 1:numel(dataTypes)
    matFile = fullfile(outDir, sprintf('compare_init_robustness_%s.mat', dataTypes{i}));
    if ~exist(matFile, 'file')
        warning('Missing results for ''%s'': %s', dataTypes{i}, matFile);
        continue
    end
    S = load(matFile);
    d = numel(results) + 1;
    results(d) = S.results;
    data(d) = S.data;
    nLoaded = nLoaded + 1;
    loaded{nLoaded} = dataTypes{i};

    if isempty(methodLabels)
        methodLabels = S.methodLabels;
        nSim = S.nSim;
        K = S.K;
        Khat = S.Khat;
        maxiter = S.maxiter;
        if isfield(S, 'constraintTol')
            constraintTol = S.constraintTol;
        end
    else
        assert(isequal(methodLabels, S.methodLabels), ...
            'methodLabels mismatch in %s', matFile);
        assert(nSim == S.nSim, 'nSim mismatch in %s', matFile);
        assert(K == S.K, 'K mismatch in %s', matFile);
        assert(Khat == S.Khat, 'Khat mismatch in %s', matFile);
    end
    fprintf('Loaded %s\n', matFile);
end

assert(~isempty(results), 'No result mats found in %s for: %s', ...
    outDir, strjoin(dataTypes, ', '));
loaded = loaded(1:nLoaded);
nData = numel(results);
nMethod = numel(methodLabels);
fprintf('Plotting %d data type(s): %s\n', nData, strjoin(loaded, ', '));

%% ------------------------------------------------------------------------
%  Running time / iterations / constraint / summary tables (log)
%  -------------------------------------------------------------------------
fprintf('\n===== Running time (seconds) =====\n');
fprintf('%-12s %-34s %10s %12s\n', 'data', 'stage', 'total', 'per fit');
for d = 1:nData
    fprintf('%-12s %-34s %10.1f %12.1f\n', results(d).label, ...
        'FlexMF fit, SeqNMF init (+SeqNMF)', results(d).time_flex(1), ...
        results(d).time_flex(1)/nSim);
    fprintf('%-12s %-34s %10.1f %12.1f\n', results(d).label, ...
        'FlexMF fit, rand init', results(d).time_flex(2), ...
        results(d).time_flex(2)/nSim);
    for m = 1:nMethod
        fprintf('%-12s %-34s %10.1f %12.1f\n', results(d).label, ...
            sprintf('matching, %s', methodLabels{m}), results(d).time_match(m), ...
            mean(results(d).time_match_runs(:,m)));
    end
    fprintf('%-12s %-34s %10.1f %12s\n', results(d).label, 'data type total', ...
        results(d).time_total, '-');
end

fprintf('\n===== FlexMF iterations run (maxiter = %d) =====\n', maxiter);
fprintf('%-12s %-22s %8s %6s %6s %14s\n', ...
    'data', 'method', 'median', 'min', 'max', '#stopped early');
for d = 1:nData
    for m = 2:3
        it = results(d).nIters(:, m-1);
        fprintf('%-12s %-22s %8.1f %6d %6d %14d\n', results(d).label, ...
            methodLabels{m}, median(it), min(it), max(it), sum(it < maxiter));
    end
end

fprintf('\n===== Constraint validation: ||Xcorr - R - W*H||_1 / ||X||_1 =====\n');
fprintf('%-12s %-22s %10s %10s %10s %8s\n', ...
    'data', 'method', 'median', 'max', 'min', '#flagged');
nFlagged = 0;
for d = 1:nData
    for m = 2:3
        c = results(d).constraints_rel(:,m);
        flagged = sum(c > constraintTol);
        nFlagged = nFlagged + flagged;
        fprintf('%-12s %-22s %10.3e %10.3e %10.3e %8d\n', ...
            results(d).label, methodLabels{m}, median(c), max(c), min(c), flagged);
    end
end
if nFlagged == 0
    fprintf('All FlexMF runs satisfy the constraint below %.3g of ||X||_1.\n', constraintTol);
else
    warning('%d FlexMF run(s) exceeded the constraint tolerance %.3g.', ...
        nFlagged, constraintTol);
end

fprintf('\n===== Median over %d restarts =====\n', nSim);
fprintf('%-12s %-22s %10s %10s %10s\n', 'data', 'method', 'EMD W', 'EMD H', '#seq');
for d = 1:nData
    for m = 1:nMethod
        fprintf('%-12s %-22s %10.4g %10.4g %10.1f\n', results(d).label, methodLabels{m}, ...
            median(results(d).emds_W(:,m), 'omitnan'), ...
            median(results(d).emds_H(:,m), 'omitnan'), ...
            median(results(d).nDetected(:,m)));
    end
end

%% ------------------------------------------------------------------------
%  Constraint residual across data types
%  -------------------------------------------------------------------------
figure;
hold on
offsets = linspace(-0.15, 0.15, 2);
for m = 2:3
    xs = [];
    ys = [];
    for d = 1:nData
        c = results(d).constraints_rel(:,m);
        xs = [xs; d + offsets(m-1)*ones(numel(c),1)]; %#ok<AGROW>
        ys = [ys; c]; %#ok<AGROW>
    end
    swarmchart(xs, ys, 25, 'filled', 'DisplayName', methodLabels{m});
end
yline(constraintTol, 'k--', 'tolerance', 'HandleVisibility', 'off');
set(gca, 'XTick', 1:nData, 'XTickLabel', {results.label}, 'YScale', 'log')
ylabel('||Xcorr - R - W*H||_1 / ||X||_1')
title('EMD constraint residual per run', 'FontSize', 14)
legend('Location', 'best')
set(gcf, 'Position', [100 100 700 420])
save2pdf(fullfile(outDir, 'compare_init_constraint_validation.pdf'), gcf)

%% ------------------------------------------------------------------------
%  One comparison figure per data type
%  -------------------------------------------------------------------------
for d = 1:nData
    emds_W = results(d).emds_W;
    emds_H = results(d).emds_H;
    nDetected = results(d).nDetected;

    figure;
    ax1 = subplot('Position', [0.13 0.68 0.8 0.26]);
    hold on
    for m = 1:nMethod
        swarmchart(m*ones(nSim,1), emds_W(:,m), 30, 'filled')
    end
    set(gca, 'XTickLabel', [])
    xlim([0.5, nMethod+0.5])
    ylabel(ax1, 'EMD of W')
    title(sprintf('%s  (K=%d, Khat=%d, %d restarts)', ...
        results(d).label, K, Khat, nSim), 'FontSize', 14)

    ax2 = subplot('Position', [0.13 0.38 0.8 0.26]);
    hold on
    for m = 1:nMethod
        swarmchart(m*ones(nSim,1), emds_H(:,m), 30, 'filled')
    end
    set(gca, 'XTickLabel', [])
    xlim([0.5, nMethod+0.5])
    ylabel(ax2, 'EMD of H')

    ax3 = subplot('Position', [0.13 0.12 0.8 0.18]);
    hold on
    errorbar(1:nMethod, median(nDetected), ...
        median(nDetected)-prctile(nDetected,25), ...
        prctile(nDetected,75)-median(nDetected), ...
        '-', 'Marker', '.', 'MarkerSize', 14, 'Color', 'k');
    yline(K, 'r--', 'true K', 'LabelHorizontalAlignment', 'left');
    set(gca, 'XTick', 1:nMethod, 'XTickLabel', methodLabels)
    xlim([0.5, nMethod+0.5])
    ylim([0, Khat])
    ylabel('#Sequences')

    linkaxes([ax1, ax2, ax3], 'x')
    set(gcf, 'Position', [100, 100, 600, 800])
    save2pdf(fullfile(outDir, sprintf('compare_init_%s.pdf', results(d).name)), gcf)
end

%% ------------------------------------------------------------------------
%  Example run visualization (every loaded data type)
%  -------------------------------------------------------------------------
for dEx = 1:nData
    Xex = data(dEx).X;
    info_warm = results(dEx).info_warm;
    info_rand = results(dEx).info_rand;
    nWarm = info_warm.best_idx;
    nRand = info_rand.best_idx;
    exW = {info_warm.seqNMF_W{nWarm}, info_warm.W_all{nWarm}, info_rand.W_all{nRand}};
    exH = {info_warm.seqNMF_H{nWarm}, info_warm.H_all{nWarm}, info_rand.H_all{nRand}};

    fprintf('\nExample: %s, best warm-start restart %d, best random restart %d\n', ...
        data(dEx).label, nWarm, nRand);

    figure;
    SimpleWHPlot(data(dEx).Wtrue, data(dEx).Htrue, 'Data', Xex, 'plotAll', plotAll);
    title(sprintf('%s: ground truth', data(dEx).label), 'FontSize', 16)
    set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8])
    save2pdf(fullfile(outDir, sprintf('compare_init_example_%s_truth.pdf', data(dEx).name)), gcf)

    for m = 1:nMethod
        figure;
        SimpleWHPlot_patch(exW{m}, exH{m}, 'plotAll', plotAll);
        title(sprintf('%s: %s reconstruction', data(dEx).label, methodLabels{m}), 'FontSize', 16)
        set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8])
        save2pdf(fullfile(outDir, sprintf('compare_init_example_%s_method%d_recon.pdf', ...
            data(dEx).name, m)), gcf)

        figure;
        SimpleWHPlot_patch(exW{m}, exH{m}, 'Data', Xex, 'plotAll', plotAll);
        title(sprintf('%s: %s factors, with raw data', data(dEx).label, methodLabels{m}), 'FontSize', 16)
        set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8])
        save2pdf(fullfile(outDir, sprintf('compare_init_example_%s_method%d_data.pdf', ...
            data(dEx).name, m)), gcf)
    end

    figure;
    plot_MR(info_warm.M_all{nWarm}, info_warm.R_all{nWarm})
    sgtitle(sprintf('%s: FlexMF SeqNMF init, M and R', data(dEx).label), 'FontSize', 16)
    save2pdf(fullfile(outDir, sprintf('compare_init_example_%s_MR_seqnmfinit.pdf', ...
        data(dEx).name)), gcf)

    figure;
    plot_MR(info_rand.M_all{nRand}, info_rand.R_all{nRand})
    sgtitle(sprintf('%s: FlexMF rand init, M and R', data(dEx).label), 'FontSize', 16)
    save2pdf(fullfile(outDir, sprintf('compare_init_example_%s_MR_randinit.pdf', ...
        data(dEx).name)), gcf)
end

fprintf('\nFigures written to %s/\n', outDir);
