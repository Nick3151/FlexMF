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
outDir = 'Simulation_EMD';
dataTypes = {'clean', 'noise', 'warp', 'jitter', 'participation', ...
             'warpnoise', 'jitternoise'};
% dataTypes = {'noise'};   % plot a subset

constraintTolDefault = 0.05;
plotAll = 1;

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% Load available results
results = struct([]);
data = struct([]);
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
    r = S.results;
    if ~isfield(r, 'nSignificant'), r.nSignificant = nan(size(r.nDetected)); end
    if ~isfield(r, 'info_warm_reseed'), r.info_warm_reseed = struct(); end
    if ~isfield(r, 'time_sig'), r.time_sig = []; end
    if ~isfield(r, 'time_sig_runs'), r.time_sig_runs = []; end
    if d == 1
        results = r;
    else
        results(d) = r; %#ok<AGROW>
    end
    dd = S.data;
    if ~isfield(dd, 'Xtrain'), dd.Xtrain = dd.X; end
    if ~isfield(dd, 'Xtest'), dd.Xtest = []; end
    if ~isfield(dd, 'Hfull'), dd.Hfull = dd.Htrue; end
    if ~isfield(dd, 'tSplit'), dd.tSplit = []; end
    if d == 1
        data = dd;
    else
        data(d) = dd; %#ok<AGROW>
    end
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
flexMethodIdx = 2:nMethod;
fprintf('Plotting %d data type(s): %s\n', nData, strjoin(loaded, ', '));

flexFitLabels = { ...
    'FlexMF fit, SeqNMF init no reseed (+SeqNMF)', ...
    'FlexMF fit, SeqNMF init + reseed', ...
    'FlexMF fit, rand init no reseed'};

%% ------------------------------------------------------------------------
%  Running time / iterations / constraint / summary tables (log)
%  -------------------------------------------------------------------------
fprintf('\n===== Running time (seconds) =====\n');
fprintf('%-12s %-40s %10s %12s\n', 'data', 'stage', 'total', 'per fit');
for d = 1:nData
    tf = results(d).time_flex;
    for fi = 1:min(numel(tf), numel(flexFitLabels))
        fprintf('%-12s %-40s %10.1f %12.1f\n', results(d).label, ...
            flexFitLabels{fi}, tf(fi), tf(fi)/nSim);
    end
    for m = 1:nMethod
        fprintf('%-12s %-40s %10.1f %12.1f\n', results(d).label, ...
            sprintf('matching, %s', methodLabels{m}), results(d).time_match(m), ...
            mean(results(d).time_match_runs(:,m)));
    end
    if ~isempty(results(d).time_sig)
        for m = 1:nMethod
            fprintf('%-12s %-40s %10.1f %12.1f\n', results(d).label, ...
                sprintf('significance, %s', methodLabels{m}), results(d).time_sig(m), ...
                mean(results(d).time_sig_runs(:,m)));
        end
    end
    fprintf('%-12s %-40s %10.1f %12s\n', results(d).label, 'data type total', ...
        results(d).time_total, '-');
end

fprintf('\n===== FlexMF iterations run (maxiter = %d) =====\n', maxiter);
fprintf('%-12s %-36s %8s %6s %6s %14s\n', ...
    'data', 'method', 'median', 'min', 'max', '#stopped early');
for d = 1:nData
    for mi = 1:numel(flexMethodIdx)
        m = flexMethodIdx(mi);
        if mi > size(results(d).nIters, 2)
            continue
        end
        it = results(d).nIters(:, mi);
        fprintf('%-12s %-36s %8.1f %6d %6d %14d\n', results(d).label, ...
            methodLabels{m}, median(it), min(it), max(it), sum(it < maxiter));
    end
end

fprintf('\n===== Constraint validation: ||Xcorr - R - W*H||_1 / ||X||_1 =====\n');
fprintf('%-12s %-36s %10s %10s %10s %8s\n', ...
    'data', 'method', 'median', 'max', 'min', '#flagged');
nFlagged = 0;
for d = 1:nData
    for m = flexMethodIdx
        c = results(d).constraints_rel(:,m);
        flagged = sum(c > constraintTol);
        nFlagged = nFlagged + flagged;
        fprintf('%-12s %-36s %10.3e %10.3e %10.3e %8d\n', ...
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
fprintf('%-12s %-36s %10s %10s %10s %10s\n', ...
    'data', 'method', 'EMD W', 'EMD H', '#seq', '#sig');
for d = 1:nData
    for m = 1:nMethod
        fprintf('%-12s %-36s %10.4g %10.4g %10.1f %10.1f\n', ...
            results(d).label, methodLabels{m}, ...
            median(results(d).emds_W(:,m), 'omitnan'), ...
            median(results(d).emds_H(:,m), 'omitnan'), ...
            median(results(d).nDetected(:,m)), ...
            median(results(d).nSignificant(:,m)));
    end
end

%% ------------------------------------------------------------------------
%  Constraint residual across data types
%  -------------------------------------------------------------------------
nFlex = numel(flexMethodIdx);
figure;
hold on
offsets = linspace(-0.2, 0.2, nFlex);
for mi = 1:nFlex
    m = flexMethodIdx(mi);
    xs = [];
    ys = [];
    for d = 1:nData
        c = results(d).constraints_rel(:,m);
        xs = [xs; d + offsets(mi)*ones(numel(c),1)]; %#ok<AGROW>
        ys = [ys; c]; %#ok<AGROW>
    end
    swarmchart(xs, ys, 25, 'filled', 'DisplayName', methodLabels{m});
end
yline(constraintTol, 'k--', 'tolerance', 'HandleVisibility', 'off');
set(gca, 'XTick', 1:nData, 'XTickLabel', {results.label}, 'YScale', 'log')
ylabel('||Xcorr - R - W*H||_1 / ||X||_1')
title('EMD constraint residual per run', 'FontSize', 14)
legend('Location', 'best')
set(gcf, 'Position', [100 100 800 420])
export_vector_pdf(fullfile(outDir, 'compare_init_constraint_validation.pdf'), gcf);

%% ------------------------------------------------------------------------
%  One comparison figure per data type
%  -------------------------------------------------------------------------
for d = 1:nData
    emds_W = results(d).emds_W;
    emds_H = results(d).emds_H;
    nDetected = results(d).nDetected;
    nSignificant = results(d).nSignificant;

    figure;
    ax1 = subplot('Position', [0.13 0.74 0.8 0.20]);
    hold on
    for m = 1:nMethod
        swarmchart(m*ones(nSim,1), emds_W(:,m), 30, 'filled')
    end
    set(gca, 'XTickLabel', [])
    xlim([0.5, nMethod+0.5])
    ylabel(ax1, 'EMD of W')
    title(sprintf('%s  (K=%d, Khat=%d, %d restarts)', ...
        results(d).label, K, Khat, nSim), 'FontSize', 14)

    ax2 = subplot('Position', [0.13 0.52 0.8 0.20]);
    hold on
    for m = 1:nMethod
        swarmchart(m*ones(nSim,1), emds_H(:,m), 30, 'filled')
    end
    set(gca, 'XTickLabel', [])
    xlim([0.5, nMethod+0.5])
    ylabel(ax2, 'EMD of H')

    ax3 = subplot('Position', [0.13 0.30 0.8 0.18]);
    hold on
    errorbar(1:nMethod, median(nDetected), ...
        median(nDetected)-prctile(nDetected,25), ...
        prctile(nDetected,75)-median(nDetected), ...
        '-', 'Marker', '.', 'MarkerSize', 14, 'Color', 'k');
    yline(K, 'r--', 'true K', 'LabelHorizontalAlignment', 'left');
    set(gca, 'XTickLabel', [])
    xlim([0.5, nMethod+0.5])
    ylim([0, Khat])
    ylabel('#Sequences')

    ax4 = subplot('Position', [0.13 0.08 0.8 0.18]);
    hold on
    errorbar(1:nMethod, median(nSignificant), ...
        median(nSignificant)-prctile(nSignificant,25), ...
        prctile(nSignificant,75)-median(nSignificant), ...
        '-', 'Marker', '.', 'MarkerSize', 14, 'Color', 'k');
    yline(K, 'r--', 'true K', 'LabelHorizontalAlignment', 'left');
    set(gca, 'XTick', 1:nMethod, 'XTickLabel', methodLabels)
    xlim([0.5, nMethod+0.5])
    ylim([0, Khat])
    ylabel('#Significant')

    linkaxes([ax1, ax2, ax3, ax4], 'x')
    set(gcf, 'Position', [100, 100, 700, 900])
    export_vector_pdf(fullfile(outDir, sprintf('compare_init_%s.pdf', results(d).name)), gcf);
end

%% ------------------------------------------------------------------------
%  Example run visualization (one selected data type)
%  -------------------------------------------------------------------------
exampleDataType = 'participation';   % must be among loaded results
exampleRunId = 2;           % [] => each method's best_idx; else fixed restart index
saveExampleFigs = false;       % false: plot only, skip PDF export

[exKnown, dEx] = ismember(exampleDataType, {results.name});
assert(exKnown, 'exampleDataType ''%s'' not among loaded results: %s', ...
    exampleDataType, strjoin({results.name}, ', '));

if isfield(data(dEx), 'Xtrain') && ~isempty(data(dEx).Xtrain)
    Xex = data(dEx).Xtrain;
else
    Xex = data(dEx).X;
end
info_warm = results(dEx).info_warm;
info_warm_reseed = results(dEx).info_warm_reseed;
info_rand = results(dEx).info_rand;

if isempty(exampleRunId)
    nWarm = info_warm.best_idx;
    nRand = info_rand.best_idx;
    if isfield(info_warm_reseed, 'best_idx')
        nWarmR = info_warm_reseed.best_idx;
    else
        nWarmR = nWarm;
    end
else
    assert(isscalar(exampleRunId) && exampleRunId >= 1 && exampleRunId <= nSim, ...
        'exampleRunId must be an integer in 1:%d', nSim);
    nWarm = exampleRunId;
    nRand = exampleRunId;
    nWarmR = exampleRunId;
end

% Method order: SeqNMF, rand, warm no-reseed, warm+reseed
exW = {info_warm.seqNMF_W{nWarm}, info_rand.W_all{nRand}, ...
    info_warm.W_all{nWarm}, info_warm_reseed.W_all{nWarmR}};
exH = {info_warm.seqNMF_H{nWarm}, info_rand.H_all{nRand}, ...
    info_warm.H_all{nWarm}, info_warm_reseed.H_all{nWarmR}};
exRunIds = [nWarm, nRand, nWarm, nWarmR];

exSig = cell(1, nMethod);
hasSig = isfield(results(dEx), 'is_significant') && ...
    ~isempty(results(dEx).is_significant);
for m = 1:nMethod
    if hasSig
        exSig{m} = results(dEx).is_significant{exRunIds(m), m};
    else
        exSig{m} = [];
    end
end

methodFileTags = cell(1, nMethod);
for m = 1:nMethod
    tag = regexprep(methodLabels{m}, '[^\w]+', '_');
    tag = regexprep(tag, '_+', '_');
    tag = regexprep(tag, '^_|_$', '');
    methodFileTags{m} = tag;
end

fprintf('\nExample: %s, run ids [SeqNMF=%d, rand=%d, warm=%d, warm+reseed=%d]\n', ...
    data(dEx).label, nWarm, nRand, nWarm, nWarmR);

figure;
SimpleWHPlot_patch(data(dEx).Wtrue, data(dEx).Htrue, 'Data', Xex, ...
    'is_significant', ones(1, size(data(dEx).Wtrue, 2)), 'plotAll', plotAll);
title(sprintf('%s: ground truth (train)', data(dEx).label), 'FontSize', 16)
set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8])
if saveExampleFigs
    export_vector_pdf(fullfile(outDir, sprintf('compare_init_example_%s_truth.pdf', data(dEx).name)), gcf);
end

for m = 1:nMethod
    if isfield(results(dEx), 'ids_match') && ~isempty(results(dEx).ids_match) ...
            && ~isempty(results(dEx).ids_match{exRunIds(m), m})
        ids = results(dEx).ids_match{exRunIds(m), m};
    else
        [~, ~, ids] = helper.similarity_WH_EMD( ...
            data(dEx).Wtrue, data(dEx).Htrue, exW{m}, exH{m});
    end
    [Wp, Hp, factorOrder] = helper.sort_matched_factors(exW{m}, exH{m}, ids);
    if ~isempty(exSig{m})
        sigPlot = exSig{m}(factorOrder);
    else
        sigPlot = [];
    end

    figure;
    SimpleWHPlot_patch(Wp, Hp, 'is_significant', sigPlot, 'plotAll', plotAll);
    title(sprintf('%s: %s reconstruction (run %d, GT order)', ...
        data(dEx).label, methodLabels{m}, exRunIds(m)), 'FontSize', 16)
    set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8])
    if saveExampleFigs
        export_vector_pdf(fullfile(outDir, sprintf('compare_init_example_%s_%s_recon.pdf', ...
            data(dEx).name, methodFileTags{m})), gcf);
    end

end

figure;
plot_MR(info_rand.M_all{nRand}, info_rand.R_all{nRand})
sgtitle(sprintf('%s: FlexMF rand (no reseed), M and R (run %d)', ...
    data(dEx).label, nRand), 'FontSize', 16)
if saveExampleFigs
    export_vector_pdf(fullfile(outDir, sprintf('compare_init_example_%s_%s_MR.pdf', ...
        data(dEx).name, methodFileTags{2})), gcf);
end

figure;
plot_MR(info_warm.M_all{nWarm}, info_warm.R_all{nWarm})
sgtitle(sprintf('%s: FlexMF SeqNMF init (no reseed), M and R (run %d)', ...
    data(dEx).label, nWarm), 'FontSize', 16)
if saveExampleFigs
    export_vector_pdf(fullfile(outDir, sprintf('compare_init_example_%s_%s_MR.pdf', ...
        data(dEx).name, methodFileTags{3})), gcf);
end

if isfield(info_warm_reseed, 'M_all') && ~isempty(info_warm_reseed.M_all)
    figure;
    plot_MR(info_warm_reseed.M_all{nWarmR}, info_warm_reseed.R_all{nWarmR})
    sgtitle(sprintf('%s: FlexMF SeqNMF init (reseed), M and R (run %d)', ...
        data(dEx).label, nWarmR), 'FontSize', 16)
    if saveExampleFigs
        export_vector_pdf(fullfile(outDir, sprintf('compare_init_example_%s_%s_MR.pdf', ...
            data(dEx).name, methodFileTags{4})), gcf);
    end
end

if saveExampleFigs
    fprintf('\nFigures written to %s/\n', outDir);
else
    fprintf('\nExample figures plotted (not saved; set saveExampleFigs=true to export)\n');
end
