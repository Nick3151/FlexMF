function results = compare_init_robustness(dataType, varargin)
%COMPARE_INIT_ROBUSTNESS  Init robustness of SeqNMF vs FlexMF for one data type
%
%   results = compare_init_robustness(dataType)
%   results = compare_init_robustness(dataType, 'Name', value, ...)
%
% Fits three methods on one synthetic corruption of the data:
%   1. SeqNMF
%   2. FlexMF (EMD) warm-started from SeqNMF
%   3. FlexMF (EMD) from random initialization
% and scores them against ground truth with helper.similarity_WH_EMD.
% Methods 1-2 come from one FlexMF_multistart warm-start call; method 3 from
% a second random-init call. Results are written to
%   Simulation_Results/compare_init_robustness_<dataType>.mat
% Plotting is separate: see compare_init_robustness_plot.m.
%
% dataType is one of:
%   'clean' | 'noise' | 'warp' | 'jitter' | 'participation'
%
% Optional name/value pairs:
%   'quick'           false   Reduced T / nSim / maxiter for a smoke test
%   'nSim'            10      Restarts per method
%   'Khat'            5       Overcomplete number of factors
%   'maxiter'         50      FlexMF iteration budget
%   'seqNMF_maxiter'  50      SeqNMF warm-start budget
%   'lambda'          0.1
%   'lambda_M'        0.1
%   'lambda_R'        1
%   'constraintTol'   0.05    Relative constraint residual flag threshold
%   'simOpts'         {}      Extra args for helper.similarity_WH_EMD
%   'outDir'          'Simulation_Results'
%   'saveResults'     true
%
% Rockfish job array (see compare_init_robustness.sh):
%   dataTypes = {'clean','noise','warp','jitter','participation'};
%   id = str2double(getenv('SLURM_ARRAY_TASK_ID'));
%   compare_init_robustness(dataTypes{id});

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
addParameter(p, 'outDir', 'Simulation_Results', @(s) ischar(s) || isstring(s));
addParameter(p, 'saveResults', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'seed', 1, @isnumeric);
addParameter(p, 'T', 2000, @isnumeric);
addParameter(p, 'K', 3, @isnumeric);
parse(p, dataType, varargin{:});
opt = p.Results;
dataType = char(opt.dataType);
opt.quick = logical(opt.quick);
opt.saveResults = logical(opt.saveResults);
opt.outDir = char(opt.outDir);

%% Ground-truth / corruption defaults
K = opt.K;
T = opt.T;
Nneurons = 5*ones(K,1);
Dt = 3.*ones(K,1);
noise = .001;
jitter = 2*ones(K,1);
warp = 2;
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

if opt.quick
    T = 400;
    nSim = min(nSim, 2);
    maxiter = min(maxiter, 5);
    seqNMF_maxiter = min(seqNMF_maxiter, 10);
    Khat = min(Khat, 4);
    if isempty(simOpts)
        simOpts = {'MaxShift', 5};
    end
end

%% Resolve data type
allNames  = {'clean', 'noise', 'warp', 'jitter', 'participation'};
allLabels = {'X (clean)', 'Xnoise', 'Xwarp', 'Xjit', 'Xpart'};
allArgs   = { {'noise', 0}, ...
              {'noise', noise}, ...
              {'noise', noise, 'warp', warp}, ...
              {'noise', noise, 'jitter', jitter}, ...
              {'noise', 0, 'participation', participation} };
[known, loc] = ismember(dataType, allNames);
assert(known, 'Unknown dataType ''%s''. Choose one of: %s', ...
    dataType, strjoin(allNames, ', '));
dataLabel = allLabels{loc};
dataArgs = allArgs{loc};

methodLabels = {'SeqNMF', 'FlexMF SeqNMF init', 'FlexMF rand init'};
nMethod = numel(methodLabels);
meanOmitNan = @(v) mean(v, 'omitnan');

flexOpts = {'EMD', 1, 'lambda', lambda, 'lambda_M', lambda_M, ...
    'lambda_R', lambda_R, 'lambdaL1H', lambdaL1H, ...
    'maxiter', maxiter, 'tolerance', tolerance, ...
    'neg_prop', 0, 'Reweight', 1};

fprintf('compare_init_robustness: dataType=%s (%s)\n', dataType, dataLabel);
fprintf('  T=%d  K=%d  Khat=%d  nSim=%d  maxiter=%d\n', T, K, Khat, nSim, maxiter);

%% Generate data
[Xd, Wd, Hd, ~] = generate_data(T, Nneurons, Dt, dataArgs{:}, ...
    'seed', seed, 'len_burst', 1, 'dynamic', 0);
frob_norm = norm(Xd(:));
Xd = Xd/frob_norm*Khat;
Wd = Wd/frob_norm*Khat;
Ld = size(Wd, 3);

data = struct();
data.name = dataType;
data.label = dataLabel;
data.X = Xd;
data.Wtrue = Wd;
data.Htrue = Hd;
data.L = Ld;
fprintf('%-8s : %d neurons x %d bins, L=%d\n', dataType, size(Xd,1), size(Xd,2), Ld);

%% Fit
fprintf('\n================ %s ================\n', dataLabel);
tDataset = tic;
time_flex = nan(1, 2);

fprintf('-- FlexMF with SeqNMF warm-start (%d restarts)\n', nSim);
tFit = tic;
[~, ~, ~, ~, ~, ~, ~, ~, info_warm] = FlexMF_multistart(Xd, ...
    'K', Khat, 'L', Ld, flexOpts{:}, ...
    'nRestarts', nSim, 'warmStart', true, ...
    'seqNMF_maxiter', seqNMF_maxiter, 'verbal', 1);
time_flex(1) = toc(tFit);
fprintf('   FlexMF (SeqNMF init, incl. SeqNMF): %.1f s total, %.1f s/restart\n', ...
    time_flex(1), time_flex(1)/nSim);

fprintf('-- FlexMF with random init (%d restarts)\n', nSim);
tFit = tic;
[~, ~, ~, ~, ~, ~, ~, ~, info_rand] = FlexMF_multistart(Xd, ...
    'K', Khat, 'L', Ld, flexOpts{:}, ...
    'nRestarts', nSim, 'warmStart', false, 'verbal', 1);
time_flex(2) = toc(tFit);
fprintf('   FlexMF (rand init): %.1f s total, %.1f s/restart\n', ...
    time_flex(2), time_flex(2)/nSim);

W_runs = cell(nSim, nMethod);
H_runs = cell(nSim, nMethod);
for n = 1:nSim
    W_runs{n,1} = info_warm.seqNMF_W{n};
    H_runs{n,1} = info_warm.seqNMF_H{n};
    W_runs{n,2} = info_warm.W_all{n};
    H_runs{n,2} = info_warm.H_all{n};
    W_runs{n,3} = info_rand.W_all{n};
    H_runs{n,3} = info_rand.H_all{n};
end

%% Match
emds_W = nan(nSim, nMethod);
emds_H = nan(nSim, nMethod);
nDetected = nan(nSim, nMethod);
time_match_runs = nan(nSim, nMethod);
fprintf('-- Matching %d fits with similarity_WH_EMD\n', nSim*nMethod);
tMatchAll = tic;
for m = 1:nMethod
    for n = 1:nSim
        tMatch = tic;
        [eW, eH, ids] = helper.similarity_WH_EMD(Wd, Hd, W_runs{n,m}, H_runs{n,m}, simOpts{:});
        time_match_runs(n,m) = toc(tMatch);
        emds_W(n,m) = meanOmitNan(eW);
        emds_H(n,m) = meanOmitNan(eH);
        nDetected(n,m) = nnz(ids);
    end
    fprintf('   %-22s matching: %.1f s total, %.1f s/fit\n', methodLabels{m}, ...
        sum(time_match_runs(:,m)), mean(time_match_runs(:,m)));
end
time_match = sum(time_match_runs, 1);
fprintf('   matching total: %.1f s\n', toc(tMatchAll));

%% Constraint validation
constraints_rel = nan(nSim, nMethod);
normX1 = norm(Xd(:), 1);
for n = 1:nSim
    for m = 2:3
        if m == 2
            Mn = info_warm.M_all{n};
            Rn = info_warm.R_all{n};
        else
            Mn = info_rand.M_all{n};
            Rn = info_rand.R_all{n};
        end
        Xcorr = helper.correct_warp(Xd, Mn);
        constraint = Xcorr - Rn - helper.reconstruct(W_runs{n,m}, H_runs{n,m});
        constraints_rel(n,m) = norm(constraint(:), 1) / normX1;
    end
end

%% Pack results (scalar struct; plot script concatenates across data types)
results = struct();
results.name = dataType;
results.label = dataLabel;
results.emds_W = emds_W;
results.emds_H = emds_H;
results.nDetected = nDetected;
results.constraints_rel = constraints_rel;
results.objs = [info_warm.objs, info_rand.objs];
results.nIters = [info_warm.nIters, info_rand.nIters];
results.time_flex = time_flex;
results.time_match = time_match;
results.time_match_runs = time_match_runs;
results.time_total = toc(tDataset);
results.info_warm = info_warm;
results.info_rand = info_rand;
fprintf('-- %s done in %.1f s\n', dataLabel, results.time_total);

%% Text summary for the cluster log
fprintf('\n===== Running time (seconds) =====\n');
fprintf('%-34s %10s %12s\n', 'stage', 'total', 'per fit');
fprintf('%-34s %10.1f %12.1f\n', 'FlexMF fit, SeqNMF init (+SeqNMF)', ...
    time_flex(1), time_flex(1)/nSim);
fprintf('%-34s %10.1f %12.1f\n', 'FlexMF fit, rand init', ...
    time_flex(2), time_flex(2)/nSim);
for m = 1:nMethod
    fprintf('%-34s %10.1f %12.1f\n', sprintf('matching, %s', methodLabels{m}), ...
        time_match(m), mean(time_match_runs(:,m)));
end
fprintf('%-34s %10.1f\n', 'data type total', results.time_total);

fprintf('\n===== FlexMF iterations run (maxiter = %d) =====\n', maxiter);
fprintf('%-22s %8s %6s %6s %14s\n', 'method', 'median', 'min', 'max', '#stopped early');
for m = 2:3
    it = results.nIters(:, m-1);
    fprintf('%-22s %8.1f %6d %6d %14d\n', methodLabels{m}, ...
        median(it), min(it), max(it), sum(it < maxiter));
end

fprintf('\n===== Constraint validation: ||Xcorr - R - W*H||_1 / ||X||_1 =====\n');
fprintf('%-22s %10s %10s %10s %8s\n', 'method', 'median', 'max', 'min', '#flagged');
nFlagged = 0;
for m = 2:3
    c = constraints_rel(:,m);
    flagged = sum(c > constraintTol);
    nFlagged = nFlagged + flagged;
    fprintf('%-22s %10.3e %10.3e %10.3e %8d\n', methodLabels{m}, ...
        median(c), max(c), min(c), flagged);
end
if nFlagged == 0
    fprintf('All FlexMF runs satisfy the constraint below %.3g of ||X||_1.\n', constraintTol);
else
    warning('%d FlexMF run(s) exceeded the constraint tolerance %.3g.', ...
        nFlagged, constraintTol);
end

fprintf('\n===== Median over %d restarts =====\n', nSim);
fprintf('%-22s %10s %10s %10s\n', 'method', 'EMD W', 'EMD H', '#seq');
for m = 1:nMethod
    fprintf('%-22s %10.4g %10.4g %10.1f\n', methodLabels{m}, ...
        median(emds_W(:,m), 'omitnan'), ...
        median(emds_H(:,m), 'omitnan'), ...
        median(nDetected(:,m)));
end

%% Save
if opt.saveResults
    if ~exist(opt.outDir, 'dir')
        mkdir(opt.outDir);
    end
    resultsFile = fullfile(opt.outDir, sprintf('compare_init_robustness_%s.mat', dataType));
    save(resultsFile, 'results', 'data', 'methodLabels', 'nSim', 'K', 'Khat', ...
        'lambda', 'lambda_M', 'lambda_R', 'maxiter', 'tolerance', ...
        'constraintTol', 'seed', 'dataType', 'simOpts', '-v7.3');
    fprintf('\nSaved results to %s\n', resultsFile);
end

end
