%% Test FlexMF over a lambda x lambda_M grid (lambda_R fixed at 1)
% FlexMF is warm-started from SeqNMF factors (as in EMD_demo_compare.m).
% Xtrain/Xtest are Frobenius-normalized; significance is tested on held-out data.
clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))
rmpath(genpath(fullfile(root, 'seqNMF-master')));
addpath(genpath(fullfile(root, 'FlexMF')));

data_type = 'warp+noise';

%% Generate some synthetic data
nSim = 1;
rng(1)
seeds = randperm(1000, nSim);

K = 3;
Khat = 5;
T = 4000; % length of data to generate
Nneurons = 5*ones(K,1); % number of neurons in each sequence
Dt = 3.*ones(K,1); % gap between each member of the sequence
noise = .005;
participation = .7.*ones(K,1);
jitter = 5*ones(K,1);
warp = 5;
maxiter = 50;
lambda_SeqNMF = .01;
lambda_R = 1;

%% Run simulation on different combinations of lambda, lambda_M
nLambdas = 5;
nMs = 5;
nGrid = nLambdas * nMs;

lambdas = logspace(-1, 1, nLambdas);
lambda_Ms = logspace(-2, 0, nMs);

%% Parallel pool (local: workers run on the cores this job already holds)
nCores = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isnan(nCores) || nCores < 1
    nCores = feature('numcores');
end
nWorkers = min([nGrid, nCores]);
pool = gcp('nocreate');
if isempty(pool)
    arrayJob = getenv('SLURM_ARRAY_JOB_ID');
    if isempty(arrayJob)
        runTag = sprintf('pid%d', feature('getpid'));
    else
        runTag = sprintf('%s_%s', arrayJob, getenv('SLURM_ARRAY_TASK_ID'));
    end
    jobDir = fullfile(tempdir, sprintf('flexmf_%s_%s', ...
        helper.sanitize_name(data_type), runTag));
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

times = cell(nSim, nLambdas, nMs);
times_emd = cell(nSim, nLambdas, nMs);
times_test = cell(nSim, nLambdas, nMs);
emds_W = cell(nSim, nLambdas, nMs);
emds_H = cell(nSim, nLambdas, nMs);
num_detected = cell(nSim, nLambdas, nMs);
num_significant = cell(nSim, nLambdas, nMs);
ids_match = cell(nSim, nLambdas, nMs);
pvals_all = cell(nSim, nLambdas, nMs);
is_significant_all = cell(nSim, nLambdas, nMs);
W_hats = cell(nSim, nLambdas, nMs);
H_hats = cell(nSim, nLambdas, nMs);
Ms = cell(nSim, nLambdas, nMs);
Rs = cell(nSim, nLambdas, nMs);
constraints_rel = nan(nSim, nLambdas, nMs);
reg_costs = nan(nSim, nLambdas, nMs);
What_SeqNMFs = cell(nSim, 1);
Hhat_SeqNMFs = cell(nSim, 1);
ids_match_SeqNMF = cell(nSim, 1);
L1_Xtrain = nan(nSim, 1);
Ws = cell(nSim, 1);
Hs = cell(nSim, 1);
Xs = cell(nSim, 1);

for n = 1:nSim
    fprintf('Simulation %d/%d\n', n, nSim);
    switch data_type
        case 'jitter+noise'
            [X, W, H, ~] = generate_data(T, Nneurons, Dt, 'noise', noise, ...
                'jitter', jitter, 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
        case 'warp+noise'
            [X, W, H, ~] = generate_data(T, Nneurons, Dt, 'noise', noise, ...
                'warp', warp, 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
        otherwise
            error('Unknown data_type ''%s''.', data_type);
    end
    Xtrain = X(:, 1:round(T/2));
    Xtest = X(:, 1+round(T/2):end);
    Htrain = H(:, 1:round(T/2));

    % Normalize train/test (and W with the train scale for GT matching)
    frob_norm_train = norm(Xtrain(:));
    Xtrain = Xtrain / frob_norm_train * Khat;
    W = W / frob_norm_train * Khat;
    frob_norm_test = norm(Xtest(:));
    Xtest = Xtest / frob_norm_test * Khat;
    L1_n = norm(Xtrain(:), 1);
    L1_Xtrain(n) = L1_n;

    Ws{n} = W;
    Hs{n} = H;
    Xs{n} = X;
    L = size(W, 3);

    % SeqNMF warm-start (once per simulation)
    fprintf('  SeqNMF warm-start (lambda=%g)\n', lambda_SeqNMF);
    [What_SeqNMF, Hhat_SeqNMF] = seqNMF(Xtrain, 'K', Khat, 'L', L, ...
        'lambda', lambda_SeqNMF, 'maxiter', maxiter, 'showPlot', 0);
    What_SeqNMFs{n} = What_SeqNMF;
    Hhat_SeqNMFs{n} = Hhat_SeqNMF;
    [~, ~, ids_Seq] = helper.similarity_WH_EMD(W, Htrain, What_SeqNMF, Hhat_SeqNMF);
    ids_match_SeqNMF{n} = ids_Seq;

    % Flat outputs for valid parfor slicing, then reshape into 3D cells
    times_n = cell(nGrid, 1);
    times_emd_n = cell(nGrid, 1);
    times_test_n = cell(nGrid, 1);
    emds_W_n = cell(nGrid, 1);
    emds_H_n = cell(nGrid, 1);
    num_detected_n = cell(nGrid, 1);
    num_significant_n = cell(nGrid, 1);
    ids_match_n = cell(nGrid, 1);
    pvals_n = cell(nGrid, 1);
    is_significant_n = cell(nGrid, 1);
    W_hats_n = cell(nGrid, 1);
    H_hats_n = cell(nGrid, 1);
    Ms_n = cell(nGrid, 1);
    Rs_n = cell(nGrid, 1);
    constraints_rel_n = nan(nGrid, 1);
    reg_costs_n = nan(nGrid, 1);

    parfor g = 1:nGrid
        [Li, Mi] = ind2sub([nLambdas, nMs], g);
        lambda = lambdas(Li);
        lambda_M = lambda_Ms(Mi);
        fprintf('  FlexMF lambda=%g, lambda_M=%g, lambda_R=%g\n', ...
            lambda, lambda_M, lambda_R);

        tic
        [W_hat, H_hat, ~, ~, ~, ~, M, R] = FlexMF(Xtrain, ...
            'K', Khat, 'L', L, 'EMD', 1, 'maxiter', maxiter, ...
            'lambda', lambda, 'lambda_M', lambda_M, 'lambda_R', lambda_R, ...
            'W_init', What_SeqNMF, 'H_init', Hhat_SeqNMF, ...
            'showPlot', 0, 'verbal', 0);
        times_n{g} = toc;
        W_hats_n{g} = W_hat;
        H_hats_n{g} = H_hat;
        Ms_n{g} = M;
        Rs_n{g} = R;
        Xcorr = helper.correct_warp(Xtrain, M);
        constraint = Xcorr - R - helper.reconstruct(W_hat, H_hat);
        constraints_rel_n(g) = norm(constraint(:), 1) / L1_n;
        [~, reg_cross, ~, ~] = helper.get_FlexMF_cost(Xtrain, W_hat, H_hat);
        reg_costs_n(g) = reg_cross;

        tic
        [emds_W_n{g}, emds_H_n{g}, ids] = ...
            helper.similarity_WH_EMD(W, Htrain, W_hat, H_hat);
        times_emd_n{g} = toc;
        num_detected_n{g} = sum(ids > 0);
        ids_match_n{g} = ids;

        % Fit H/M/R on test with W fixed, then significance test
        tic
        [W_hat, ~, ~, ~, ~, ~, M_test, ~] = FlexMF(Xtest, ...
            'K', Khat, 'L', L, 'W_fixed', 1, 'W_init', W_hat, ...
            'EMD', 1, 'lambda', lambda, 'lambda_R', lambda_R, ...
            'lambda_M', lambda_M, 'maxiter', maxiter, 'showPlot', 0, 'verbal', 0);
        times_test_n{g} = toc;
        [pvals, is_significant, ~] = test_significance_EMD(Xtest, W_hat, M_test, 'plot', 0);
        pvals_n{g} = pvals;
        is_significant_n{g} = is_significant;
        num_significant_n{g} = sum(is_significant);
    end

    sz = [1, nLambdas, nMs];
    times(n, :, :) = reshape(times_n, sz);
    times_emd(n, :, :) = reshape(times_emd_n, sz);
    times_test(n, :, :) = reshape(times_test_n, sz);
    emds_W(n, :, :) = reshape(emds_W_n, sz);
    emds_H(n, :, :) = reshape(emds_H_n, sz);
    num_detected(n, :, :) = reshape(num_detected_n, sz);
    num_significant(n, :, :) = reshape(num_significant_n, sz);
    ids_match(n, :, :) = reshape(ids_match_n, sz);
    pvals_all(n, :, :) = reshape(pvals_n, sz);
    is_significant_all(n, :, :) = reshape(is_significant_n, sz);
    W_hats(n, :, :) = reshape(W_hats_n, sz);
    H_hats(n, :, :) = reshape(H_hats_n, sz);
    Ms(n, :, :) = reshape(Ms_n, sz);
    Rs(n, :, :) = reshape(Rs_n, sz);
    constraints_rel(n, :, :) = reshape(constraints_rel_n, sz);
    reg_costs(n, :, :) = reshape(reg_costs_n, sz);
end

%% Save results
out_dir = 'Simulation_EMD';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
results_file = fullfile(out_dir, sprintf( ...
    'FlexMF_choose_lambda_lambdaM_%s_lambdaR=%0.3e.mat', ...
    helper.sanitize_name(data_type), lambda_R));
settings = struct('data_type', data_type, 'K', K, 'Khat', Khat, 'T', T, ...
    'Nneurons', Nneurons, 'Dt', Dt, 'noise', noise, 'jitter', jitter, ...
    'participation', participation, 'warp', warp, 'seeds', seeds, ...
    'maxiter', maxiter, 'lambda_SeqNMF', lambda_SeqNMF, ...
    'lambda_R', lambda_R, 'lambdas', lambdas, 'lambda_Ms', lambda_Ms);
save(results_file, 'times', 'times_emd', 'times_test', ...
    'lambdas', 'lambda_Ms', 'lambda_R', ...
    'W_hats', 'H_hats', 'What_SeqNMFs', 'Hhat_SeqNMFs', 'ids_match_SeqNMF', ...
    'Ws', 'Hs', 'Xs', 'L1_Xtrain', 'Ms', 'Rs', 'constraints_rel', 'reg_costs', 'ids_match', 'emds_W', 'emds_H', ...
    'num_detected', 'num_significant', 'pvals_all', 'is_significant_all', ...
    'settings');
fprintf('Saved results: %s\n', results_file);

% Shut down the parallel pool
delete(gcp('nocreate'));
