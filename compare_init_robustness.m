function results = compare_init_robustness(dataType, varargin)
%COMPARE_INIT_ROBUSTNESS  Init robustness of SeqNMF vs FlexMF for one data type
%
%   results = compare_init_robustness(dataType)
%   results = compare_init_robustness(dataType, 'Name', value, ...)
%
% Fits four methods on the training half of one synthetic corruption:
%   1. SeqNMF
%   2. FlexMF (EMD) from random initialization, no reseed
%   3. FlexMF (EMD) warm-started from SeqNMF, no reseed
%   4. FlexMF (EMD) warm-started from SeqNMF, with reseed
% Methods 3 and 4 share the same SeqNMF init per restart. Restarts run in a
% single parfor over a local pool when UseParallel is true.
%
% Data are split train/test (half/half), each Frobenius-normalized to Khat.
% Training fits are scored against ground truth with helper.similarity_WH_EMD.
% Held-out significance uses test_significance (SeqNMF) or
% test_significance_EMD after a W-fixed FlexMF refit on Xtest (FlexMF).
% Results are written to
%   Simulation_Results/compare_init_robustness_<dataType>.mat
% Plotting is separate: see compare_init_robustness_plot.m.
%
% dataType is one of:
%   'clean' | 'noise' | 'warp' | 'jitter' | 'participation' |
%   'warpnoise' | 'jitternoise'
%
% Optional name/value pairs:
%   'quick'           false   Reduced T / nSim / maxiter; disables parallel
%   'nSim'            10      Restarts per method
%   'Khat'            5       Overcomplete number of factors
%   'maxiter'         50      FlexMF iteration budget
%   'seqNMF_maxiter'  50      SeqNMF warm-start budget
%   'lambda'          0.1
%   'lambda_M'        0.1
%   'lambda_R'        1
%   'constraintTol'   0.05    Relative constraint residual flag threshold
%   'simOpts'         {}      Extra args for helper.similarity_WH_EMD
%   'UseParallel'     true    Start a local parpool and parfor over restarts
%   'nWorkers'        10      Pool size, capped by nSim and allocated cores
%   'outDir'          'Simulation_Results'
%   'saveResults'     true
%
% Rockfish job array (see compare_init_robustness.sh):
%   dataTypes = {'clean','noise','warp','jitter','participation',...
%                'warpnoise','jitternoise'};
%   id = str2double(getenv('SLURM_ARRAY_TASK_ID'));
%   compare_init_robustness(dataTypes{id});
% Each array task must request --cpus-per-task >= nWorkers; the local pool
% runs on those cores, so no second SLURM job is queued.

%% Paths
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir)
    thisDir = pwd;
end
root = fileparts(thisDir);
addpath(fullfile(root, 'TFOCS'))
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))
if exist(fullfile(root, 'seqNMF-master'), 'dir')
    rmpath(genpath(fullfile(root, 'seqNMF-master')));
end
addpath(genpath(thisDir));

%% Parse options
p = inputParser;
p.FunctionName = 'compare_init_robustness';
addRequired(p, 'dataType', @(s) ischar(s) || isstring(s));
addParameter(p, 'quick', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'nSim', 10, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'Khat', 5, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'maxiter', 50, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'seqNMF_maxiter', 50, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'lambda', 0.1, @isnumeric);
addParameter(p, 'lambda_M', 0.1, @isnumeric);
addParameter(p, 'lambda_R', 1, @isnumeric);
addParameter(p, 'lambdaL1H', 0, @isnumeric);
addParameter(p, 'tolerance', 1e-3, @isnumeric);
addParameter(p, 'constraintTol', 0.05, @isnumeric);
addParameter(p, 'simOpts', {}, @iscell);
addParameter(p, 'UseParallel', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'nWorkers', 10, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'outDir', 'Simulation_Results', @(s) ischar(s) || isstring(s));
addParameter(p, 'saveResults', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'seed', 1, @isnumeric);
addParameter(p, 'T', 4000, @isnumeric);
addParameter(p, 'K', 3, @isnumeric);
parse(p, dataType, varargin{:});
opt = p.Results;
dataType = char(opt.dataType);
opt.quick = logical(opt.quick);
opt.saveResults = logical(opt.saveResults);
opt.UseParallel = logical(opt.UseParallel);
opt.outDir = char(opt.outDir);

%% Ground-truth / corruption defaults
K = opt.K;
T = opt.T;
Nneurons = 5*ones(K,1);
Dt = 3.*ones(K,1);
noise = .005;
jitter = 5*ones(K,1);
warp = 5;
participation = .8*ones(K,1);

nSim = opt.nSim;
Khat = opt.Khat;
maxiter = opt.maxiter;
seqNMF_maxiter = opt.seqNMF_maxiter;
lambda = opt.lambda;
lambda_M = opt.lambda_M;
lambda_R = opt.lambda_R;
lambdaL1H = opt.lambdaL1H;
tolerance = opt.tolerance;
constraintTol = opt.constraintTol;
simOpts = opt.simOpts;
seed = opt.seed;
reseedEmpty = 5;  % FlexMF default for method 4
nWorkers = opt.nWorkers;
useParallel = opt.UseParallel;

if opt.quick
    T = 800;
    nSim = min(nSim, 2);
    maxiter = min(maxiter, 5);
    seqNMF_maxiter = min(seqNMF_maxiter, 10);
    Khat = min(Khat, 4);
    useParallel = false;
    if isempty(simOpts)
        simOpts = {'MaxShift', 5};
    end
end

%% Resolve data type
allNames  = {'clean', 'noise', 'warp', 'jitter', 'participation', ...
             'warpnoise', 'jitternoise'};
allLabels = {'X (clean)', 'Xnoise', 'Xwarp', 'Xjit', 'Xpart', ...
             'Xwarpnoise', 'Xjitternoise'};
allArgs   = { {'noise', 0}, ...
              {'noise', noise}, ...
              {'noise', 0, 'warp', warp}, ...
              {'noise', 0, 'jitter', jitter}, ...
              {'noise', 0, 'participation', participation}, ...
              {'noise', noise, 'warp', warp}, ...
              {'noise', noise, 'jitter', jitter} };
[known, loc] = ismember(dataType, allNames);
assert(known, 'Unknown dataType ''%s''. Choose one of: %s', ...
    dataType, strjoin(allNames, ', '));
dataLabel = allLabels{loc};
dataArgs = allArgs{loc};

methodLabels = { ...
    'SeqNMF', ...
    'FlexMF rand (no reseed)', ...
    'FlexMF SeqNMF init (no reseed)', ...
    'FlexMF SeqNMF init (reseed)'};
nMethod = numel(methodLabels);
flexMethodIdx = 2:nMethod;  % columns that are FlexMF

fprintf('compare_init_robustness: dataType=%s (%s)\n', dataType, dataLabel);
fprintf('  T=%d  K=%d  Khat=%d  nSim=%d  maxiter=%d  UseParallel=%d\n', ...
    T, K, Khat, nSim, maxiter, useParallel);

%% Parallel pool (local: workers run on the cores this job already holds)
if useParallel
    nCores = str2double(getenv('SLURM_CPUS_PER_TASK'));
    if isnan(nCores) || nCores < 1
        nCores = feature('numcores');
    end
    nWorkers = min([nWorkers, nSim, nCores]);
    pool = gcp('nocreate');
    if isempty(pool)
        % Per-session storage keeps concurrent array tasks from colliding
        arrayJob = getenv('SLURM_ARRAY_JOB_ID');
        if isempty(arrayJob)
            runTag = sprintf('pid%d', feature('getpid'));
        else
            runTag = sprintf('%s_%s', arrayJob, getenv('SLURM_ARRAY_TASK_ID'));
        end
        jobDir = fullfile(tempdir, sprintf('flexmf_%s_%s', dataType, runTag));
        if ~exist(jobDir, 'dir')
            mkdir(jobDir);
        end
        lc = parcluster('local');
        lc.JobStorageLocation = jobDir;
        lc.NumWorkers = max(nWorkers, lc.NumWorkers);
        fprintf('Starting local parallel pool (%d workers of %d cores)...\n', ...
            nWorkers, nCores);
        parpool(lc, nWorkers);
    else
        nWorkers = pool.NumWorkers;
        fprintf('Using existing parallel pool (%d workers)\n', nWorkers);
    end
end

%% Generate data, split train/test, normalize
[Xd, Wd, Hd, ~] = generate_data(T, Nneurons, Dt, dataArgs{:}, ...
    'seed', seed, 'len_burst', 1, 'dynamic', 0);
tSplit = round(T/2);
Xtrain = Xd(:, 1:tSplit);
Xtest  = Xd(:, tSplit+1:end);
Htrain = Hd(:, 1:tSplit);

frob_train = norm(Xtrain(:));
Xtrain = Xtrain / frob_train * Khat;
Wd = Wd / frob_train * Khat;
frob_test = norm(Xtest(:));
Xtest = Xtest / frob_test * Khat;
Ld = size(Wd, 3);
normX1 = norm(Xtrain(:), 1);

data = struct();
data.name = dataType;
data.label = dataLabel;
data.X = Xtrain;           % training data (for plots / backward compat)
data.Xtrain = Xtrain;
data.Xtest = Xtest;
data.Wtrue = Wd;
data.Htrue = Htrain;
data.Hfull = Hd;
data.L = Ld;
data.tSplit = tSplit;
fprintf('%-8s : %d neurons x %d bins (train %d / test %d), L=%d\n', ...
    dataType, size(Xtrain,1), size(Xd,2), size(Xtrain,2), size(Xtest,2), Ld);

%% Allocate sliced outputs for parfor
seqNMF_W = cell(nSim, 1);
seqNMF_H = cell(nSim, 1);
W_rand = cell(nSim, 1);
H_rand = cell(nSim, 1);
M_rand = cell(nSim, 1);
R_rand = cell(nSim, 1);
W_warm = cell(nSim, 1);
H_warm = cell(nSim, 1);
M_warm = cell(nSim, 1);
R_warm = cell(nSim, 1);
W_warmR = cell(nSim, 1);
H_warmR = cell(nSim, 1);
M_warmR = cell(nSim, 1);
R_warmR = cell(nSim, 1);

objs_rand = nan(nSim, 1);
objs_warm = nan(nSim, 1);
objs_warmR = nan(nSim, 1);
nIters_rand = nan(nSim, 1);
nIters_warm = nan(nSim, 1);
nIters_warmR = nan(nSim, 1);

emds_W = nan(nSim, nMethod);
emds_H = nan(nSim, nMethod);
nDetected = nan(nSim, nMethod);
nSignificant = nan(nSim, nMethod);
pvals = cell(nSim, nMethod);
is_significant = cell(nSim, nMethod);
constraints_rel = nan(nSim, nMethod);

%% Fit / match / significance (one parfor over restarts)
fprintf('\n================ %s ================\n', dataLabel);
fprintf('-- Fitting %d restarts (4 methods + match + significance)\n', nSim);
tDataset = tic;
tParfor = tic;

if useParallel
    parfor n = 1:nSim
        [seqNMF_W{n}, seqNMF_H{n}, W_rand{n}, H_rand{n}, M_rand{n}, R_rand{n}, ...
            W_warm{n}, H_warm{n}, M_warm{n}, R_warm{n}, ...
            W_warmR{n}, H_warmR{n}, M_warmR{n}, R_warmR{n}, ...
            objs_rand(n), objs_warm(n), objs_warmR(n), ...
            nIters_rand(n), nIters_warm(n), nIters_warmR(n), ...
            emds_W(n,:), emds_H(n,:), nDetected(n,:), ...
            nSignificant(n,:), pvals(n,:), is_significant(n,:), ...
            constraints_rel(n,:)] = compare_init_one_restart( ...
            Xtrain, Xtest, Wd, Htrain, ...
            Khat, Ld, lambda, lambda_M, lambda_R, lambdaL1H, ...
            maxiter, seqNMF_maxiter, tolerance, reseedEmpty, ...
            normX1, simOpts, nMethod);
        fprintf('compare_init_robustness: restart %d/%d done\n', n, nSim);
    end
else
    for n = 1:nSim
        fprintf('compare_init_robustness: restart %d/%d\n', n, nSim);
        [seqNMF_W{n}, seqNMF_H{n}, W_rand{n}, H_rand{n}, M_rand{n}, R_rand{n}, ...
            W_warm{n}, H_warm{n}, M_warm{n}, R_warm{n}, ...
            W_warmR{n}, H_warmR{n}, M_warmR{n}, R_warmR{n}, ...
            objs_rand(n), objs_warm(n), objs_warmR(n), ...
            nIters_rand(n), nIters_warm(n), nIters_warmR(n), ...
            emds_W(n,:), emds_H(n,:), nDetected(n,:), ...
            nSignificant(n,:), pvals(n,:), is_significant(n,:), ...
            constraints_rel(n,:)] = compare_init_one_restart( ...
            Xtrain, Xtest, Wd, Htrain, ...
            Khat, Ld, lambda, lambda_M, lambda_R, lambdaL1H, ...
            maxiter, seqNMF_maxiter, tolerance, reseedEmpty, ...
            normX1, simOpts, nMethod);
    end
end

time_parfor = toc(tParfor);
% Stage timings collapsed: wall-clock for the whole parallel block
time_flex = [time_parfor, nan, nan];
time_match = nan(1, nMethod);
time_match_runs = nan(nSim, nMethod);
time_sig = nan(1, nMethod);
time_sig_runs = nan(nSim, nMethod);
fprintf('   parfor wall time: %.1f s (%.1f s/restart avg)\n', ...
    time_parfor, time_parfor/nSim);

%% Rebuild info structs for plot compatibility
info_warm = struct();
info_warm.W_all = W_warm;
info_warm.H_all = H_warm;
info_warm.M_all = M_warm;
info_warm.R_all = R_warm;
info_warm.seqNMF_W = seqNMF_W;
info_warm.seqNMF_H = seqNMF_H;
info_warm.objs = objs_warm;
info_warm.nIters = nIters_warm;
info_warm.constraints_rel = constraints_rel(:, 3);
[~, info_warm.best_idx] = min(objs_warm);

info_warm_reseed = struct();
info_warm_reseed.W_all = W_warmR;
info_warm_reseed.H_all = H_warmR;
info_warm_reseed.M_all = M_warmR;
info_warm_reseed.R_all = R_warmR;
info_warm_reseed.objs = objs_warmR;
info_warm_reseed.nIters = nIters_warmR;
info_warm_reseed.constraints_rel = constraints_rel(:, 4);
[~, info_warm_reseed.best_idx] = min(objs_warmR);

info_rand = struct();
info_rand.W_all = W_rand;
info_rand.H_all = H_rand;
info_rand.M_all = M_rand;
info_rand.R_all = R_rand;
info_rand.objs = objs_rand;
info_rand.nIters = nIters_rand;
info_rand.constraints_rel = constraints_rel(:, 2);
[~, info_rand.best_idx] = min(objs_rand);

%% Pack results (scalar struct; plot script concatenates across data types)
% objs / nIters columns align with method order 2,3,4 (rand, warm0, warm+reseed)
results = struct();
results.name = dataType;
results.label = dataLabel;
results.emds_W = emds_W;
results.emds_H = emds_H;
results.nDetected = nDetected;
results.nSignificant = nSignificant;
results.pvals = pvals;
results.is_significant = is_significant;
results.constraints_rel = constraints_rel;
results.objs = [objs_rand, objs_warm, objs_warmR];
results.nIters = [nIters_rand, nIters_warm, nIters_warmR];
results.time_flex = time_flex;
results.time_match = time_match;
results.time_match_runs = time_match_runs;
results.time_sig = time_sig;
results.time_sig_runs = time_sig_runs;
results.time_total = toc(tDataset);
results.info_warm = info_warm;
results.info_warm_reseed = info_warm_reseed;
results.info_rand = info_rand;
fprintf('-- %s done in %.1f s\n', dataLabel, results.time_total);

%% Text summary for the cluster log
fprintf('\n===== Running time (seconds) =====\n');
fprintf('%-40s %10s %12s\n', 'stage', 'total', 'per restart');
fprintf('%-40s %10.1f %12.1f\n', 'parfor (fit+match+significance)', ...
    time_parfor, time_parfor/nSim);
fprintf('%-40s %10.1f\n', 'data type total', results.time_total);

fprintf('\n===== FlexMF iterations run (maxiter = %d) =====\n', maxiter);
fprintf('%-36s %8s %6s %6s %14s\n', 'method', 'median', 'min', 'max', '#stopped early');
for mi = 1:numel(flexMethodIdx)
    m = flexMethodIdx(mi);
    it = results.nIters(:, mi);
    fprintf('%-36s %8.1f %6d %6d %14d\n', methodLabels{m}, ...
        median(it), min(it), max(it), sum(it < maxiter));
end

fprintf('\n===== Constraint validation: ||Xcorr - R - W*H||_1 / ||X||_1 =====\n');
fprintf('%-36s %10s %10s %10s %8s\n', 'method', 'median', 'max', 'min', '#flagged');
nFlagged = 0;
for m = flexMethodIdx
    c = constraints_rel(:,m);
    flagged = sum(c > constraintTol);
    nFlagged = nFlagged + flagged;
    fprintf('%-36s %10.3e %10.3e %10.3e %8d\n', methodLabels{m}, ...
        median(c), max(c), min(c), flagged);
end
if nFlagged == 0
    fprintf('All FlexMF runs satisfy the constraint below %.3g of ||X||_1.\n', constraintTol);
else
    warning('%d FlexMF run(s) exceeded the constraint tolerance %.3g.', ...
        nFlagged, constraintTol);
end

fprintf('\n===== Median over %d restarts =====\n', nSim);
fprintf('%-36s %10s %10s %10s %10s\n', 'method', 'EMD W', 'EMD H', '#seq', '#sig');
for m = 1:nMethod
    fprintf('%-36s %10.4g %10.4g %10.1f %10.1f\n', methodLabels{m}, ...
        median(emds_W(:,m), 'omitnan'), ...
        median(emds_H(:,m), 'omitnan'), ...
        median(nDetected(:,m)), ...
        median(nSignificant(:,m)));
end

%% Save
if opt.saveResults
    if ~exist(opt.outDir, 'dir')
        mkdir(opt.outDir);
    end
    resultsFile = fullfile(opt.outDir, sprintf('compare_init_robustness_%s.mat', dataType));
    save(resultsFile, 'results', 'data', 'methodLabels', 'nSim', 'K', 'Khat', ...
        'lambda', 'lambda_M', 'lambda_R', 'maxiter', 'tolerance', ...
        'constraintTol', 'seed', 'dataType', 'simOpts', 'T', 'reseedEmpty', ...
        'useParallel', 'nWorkers', '-v7.3');
    fprintf('\nSaved results to %s\n', resultsFile);
end

end

%% ------------------------------------------------------------------------
function [W0, H0, Wr, Hr, Mr, Rr, Ww, Hw, Mw, Rw, Wwr, Hwr, Mwr, Rwr, ...
    obj_r, obj_w, obj_wr, nIt_r, nIt_w, nIt_wr, ...
    eWrow, eHrow, nDetRow, nSigRow, pvalsRow, isigRow, cRelRow] = ...
    compare_init_one_restart(Xtrain, Xtest, Wd, Htrain, ...
    Khat, Ld, lambda, lambda_M, lambda_R, lambdaL1H, ...
    maxiter, seqNMF_maxiter, tolerance, reseedEmpty, ...
    normX1, simOpts, nMethod)

    eWrow = nan(1, nMethod);
    eHrow = nan(1, nMethod);
    nDetRow = nan(1, nMethod);
    nSigRow = nan(1, nMethod);
    pvalsRow = cell(1, nMethod);
    isigRow = cell(1, nMethod);
    cRelRow = nan(1, nMethod);

    % 1) SeqNMF
    [W0, H0] = seqNMF(Xtrain, 'K', Khat, 'L', Ld, ...
        'lambda', lambda, 'maxiter', seqNMF_maxiter, 'showPlot', 0);

    % 3) FlexMF warm, no reseed
    [Ww, Hw, costw, errorsw, ~, ~, Mw, Rw] = FlexMF(Xtrain, ...
        'K', Khat, 'L', Ld, 'EMD', 1, ...
        'lambda', lambda, 'lambda_M', lambda_M, 'lambda_R', lambda_R, ...
        'lambdaL1H', lambdaL1H, 'maxiter', maxiter, 'tolerance', tolerance, ...
        'neg_prop', 0, 'Reweight', 0, 'reseedEmpty', 0, ...
        'W_init', W0, 'H_init', H0, 'showPlot', 0, 'verbal', 0);
    nIt_w = numel(costw) - 1;
    obj_w = lambda * errorsw(end, 2) + ...
        lambda_M * norm(Mw(:), 1) + lambda_R * norm(Rw(:), 1);

    % 4) FlexMF warm, with reseed (same SeqNMF init)
    [Wwr, Hwr, costwr, errorswr, ~, ~, Mwr, Rwr] = FlexMF(Xtrain, ...
        'K', Khat, 'L', Ld, 'EMD', 1, ...
        'lambda', lambda, 'lambda_M', lambda_M, 'lambda_R', lambda_R, ...
        'lambdaL1H', lambdaL1H, 'maxiter', maxiter, 'tolerance', tolerance, ...
        'neg_prop', 0, 'Reweight', 0, 'reseedEmpty', reseedEmpty, ...
        'W_init', W0, 'H_init', H0, 'showPlot', 0, 'verbal', 0);
    nIt_wr = numel(costwr) - 1;
    obj_wr = lambda * errorswr(end, 2) + ...
        lambda_M * norm(Mwr(:), 1) + lambda_R * norm(Rwr(:), 1);

    % 2) FlexMF random, no reseed
    [Wr, Hr, costr, errorsr, ~, ~, Mr, Rr] = FlexMF(Xtrain, ...
        'K', Khat, 'L', Ld, 'EMD', 1, ...
        'lambda', lambda, 'lambda_M', lambda_M, 'lambda_R', lambda_R, ...
        'lambdaL1H', lambdaL1H, 'maxiter', maxiter, 'tolerance', tolerance, ...
        'neg_prop', 0, 'Reweight', 0, 'reseedEmpty', 0, ...
        'showPlot', 0, 'verbal', 0);
    nIt_r = numel(costr) - 1;
    obj_r = lambda * errorsr(end, 2) + ...
        lambda_M * norm(Mr(:), 1) + lambda_R * norm(Rr(:), 1);

    % Method order: 1 SeqNMF, 2 rand, 3 warm, 4 warm+reseed
    Wcell = {W0, Wr, Ww, Wwr};
    Hcell = {H0, Hr, Hw, Hwr};
    Mcell = {[], Mr, Mw, Mwr};
    Rcell = {[], Rr, Rw, Rwr};

    for m = 1:nMethod
        [eW, eH, ids] = helper.similarity_WH_EMD(Wd, Htrain, Wcell{m}, Hcell{m}, simOpts{:});
        eWrow(m) = mean(eW, 'omitnan');
        eHrow(m) = mean(eH, 'omitnan');
        nDetRow(m) = nnz(ids);

        if m == 1
            [pv, isig] = test_significance(Xtest, Wcell{m});
        else
            [Wtest, ~, ~, ~, ~, ~, M_test, ~] = FlexMF(Xtest, ...
                'K', Khat, 'L', Ld, 'W_fixed', 1, 'W_init', Wcell{m}, ...
                'EMD', 1, 'lambda', lambda, 'lambda_M', lambda_M, ...
                'lambda_R', lambda_R, 'maxiter', maxiter, ...
                'showPlot', 0, 'verbal', 0);
            [pv, isig, ~] = test_significance_EMD(Xtest, Wtest, M_test, 'plot', 0);
            Xcorr = helper.correct_warp(Xtrain, Mcell{m});
            constraint = Xcorr - Rcell{m} - helper.reconstruct(Wcell{m}, Hcell{m});
            cRelRow(m) = norm(constraint(:), 1) / normX1;
        end
        pvalsRow{m} = pv;
        isigRow{m} = isig;
        nSigRow(m) = sum(isig);
    end
end
