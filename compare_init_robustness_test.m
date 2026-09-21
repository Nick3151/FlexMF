%% compare_init_robustness_test
% One-simulation smoke test of the four methods used in
% compare_init_robustness, focused on whether FlexMF reseed runs.
%
% Methods (same order as compare_init_robustness):
%   1. SeqNMF
%   2. FlexMF random init, reseedEmpty=0
%   3. FlexMF SeqNMF warm-start, reseedEmpty=0
%   4. FlexMF SeqNMF warm-start, reseedEmpty=5
%
% Evidence that reseed is active: FlexMF prints
%   "Reseeded <k> unused factors at iter <t>"
% when verbal=1. This script diaries that output and counts matches.
% Overcomplete Khat > K makes unused SeqNMF factors likely so method 4
% has something to reseed.

clear all
close all
clc

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

%% Settings (small for a quick local run)
dataType = 'noise';
K = 3;
Khat = 5;
T = 800;
maxiter = 10;
seqNMF_maxiter = 20;
lambda = 0.1;
lambda_M = 0.1;
lambda_R = 1;
reseedEmpty = 5;
seed = 1;
noise = 0.005;

Nneurons = 5*ones(K,1);
Dt = 3.*ones(K,1);

methodLabels = { ...
    'SeqNMF', ...
    'FlexMF rand (no reseed)', ...
    'FlexMF SeqNMF init (no reseed)', ...
    'FlexMF SeqNMF init (reseed)'};

fprintf('compare_init_robustness_test: dataType=%s\n', dataType);
fprintf('  T=%d  K=%d  Khat=%d  maxiter=%d  reseedEmpty=%d\n', ...
    T, K, Khat, maxiter, reseedEmpty);

%% Data: train half only (no significance / test path)
[Xd, Wd, Hd, ~] = generate_data(T, Nneurons, Dt, 'noise', noise, ...
    'seed', seed, 'len_burst', 1, 'dynamic', 0);
tSplit = round(T/2);
Xtrain = Xd(:, 1:tSplit);
Htrain = Hd(:, 1:tSplit);
frob_train = norm(Xtrain(:));
Xtrain = Xtrain / frob_train * Khat;
Wd = Wd / frob_train * Khat;
Ld = size(Wd, 3);
fprintf('  train: %d neurons x %d bins, L=%d\n', size(Xtrain,1), size(Xtrain,2), Ld);

%% 1) SeqNMF (shared warm-start for methods 3 and 4)
fprintf('\n== 1) SeqNMF ==\n');
[W0, H0] = seqNMF(Xtrain, 'K', Khat, 'L', Ld, ...
    'lambda', lambda, 'maxiter', seqNMF_maxiter, 'showPlot', 0);
nEmpty0 = count_near_empty_factors(W0, H0, Xtrain);
recon0 = zeros(1, Khat);
for kk = 1:Khat
    Xk = helper.reconstruct(W0(:, kk, :), H0(kk, :));
    recon0(kk) = norm(Xk(:));
end
fprintf('  SeqNMF ||W_k(*)H_k||: '); fprintf('%.3e ', recon0); fprintf('\n');
fprintf('  SeqNMF near-empty factors: %d / %d\n', nEmpty0, Khat);

diaryFile = fullfile(tempdir, 'compare_init_robustness_test_diary.txt');
if exist(diaryFile, 'file')
    delete(diaryFile);
end

%% 2) FlexMF random, no reseed
fprintf('\n== 2) FlexMF rand, reseedEmpty=0 ==\n');
diary(diaryFile);
[Wr, Hr] = FlexMF(Xtrain, 'K', Khat, 'L', Ld, 'EMD', 1, ...
    'lambda', lambda, 'lambda_M', lambda_M, 'lambda_R', lambda_R, ...
    'maxiter', maxiter, 'neg_prop', 0, 'Reweight', 0, ...
    'reseedEmpty', 0, 'showPlot', 0, 'verbal', 1);
diary off
nReseedRand = count_reseed_msgs(diaryFile);
nEmptyRand = count_near_empty_factors(Wr, Hr, Xtrain);
fprintf('  reseed messages: %d\n', nReseedRand);
fprintf('  near-empty factors after fit: %d / %d\n', nEmptyRand, Khat);

%% 3) FlexMF warm, no reseed
fprintf('\n== 3) FlexMF SeqNMF init, reseedEmpty=0 ==\n');
if exist(diaryFile, 'file'), delete(diaryFile); end
diary(diaryFile);
[Ww, Hw] = FlexMF(Xtrain, 'K', Khat, 'L', Ld, 'EMD', 1, ...
    'lambda', lambda, 'lambda_M', lambda_M, 'lambda_R', lambda_R, ...
    'maxiter', maxiter, 'neg_prop', 0, 'Reweight', 0, ...
    'reseedEmpty', 0, 'W_init', W0, 'H_init', H0, ...
    'showPlot', 0, 'verbal', 1);
diary off
nReseedWarm0 = count_reseed_msgs(diaryFile);
nEmptyWarm0 = count_near_empty_factors(Ww, Hw, Xtrain);
fprintf('  reseed messages: %d\n', nReseedWarm0);
fprintf('  near-empty factors after fit: %d / %d\n', nEmptyWarm0, Khat);

%% 4) FlexMF warm, with reseed (same SeqNMF init)
fprintf('\n== 4) FlexMF SeqNMF init, reseedEmpty=%d ==\n', reseedEmpty);
if exist(diaryFile, 'file'), delete(diaryFile); end
diary(diaryFile);
[Wwr, Hwr] = FlexMF(Xtrain, 'K', Khat, 'L', Ld, 'EMD', 1, ...
    'lambda', lambda, 'lambda_M', lambda_M, 'lambda_R', lambda_R, ...
    'maxiter', maxiter, 'neg_prop', 0, 'Reweight', 0, ...
    'reseedEmpty', reseedEmpty, 'W_init', W0, 'H_init', H0, ...
    'showPlot', 0, 'verbal', 1);
diary off
nReseedWarmR = count_reseed_msgs(diaryFile);
nEmptyWarmR = count_near_empty_factors(Wwr, Hwr, Xtrain);
fprintf('  reseed messages: %d\n', nReseedWarmR);
fprintf('  near-empty factors after fit: %d / %d\n', nEmptyWarmR, Khat);

%% Summary / pass criteria
fprintf('\n===== Reseed test summary =====\n');
fprintf('%-40s %10s %12s\n', 'method', '#reseed msgs', '#near-empty');
fprintf('%-40s %10s %12d\n', methodLabels{1}, 'n/a', nEmpty0);
fprintf('%-40s %10d %12d\n', methodLabels{2}, nReseedRand, nEmptyRand);
fprintf('%-40s %10d %12d\n', methodLabels{3}, nReseedWarm0, nEmptyWarm0);
fprintf('%-40s %10d %12d\n', methodLabels{4}, nReseedWarmR, nEmptyWarmR);

okNoReseed = (nReseedRand == 0) && (nReseedWarm0 == 0);
okReseed = (nReseedWarmR > 0);
if okNoReseed && okReseed
    fprintf('\nPASS: methods 2-3 emitted no reseed messages; method 4 did (%d).\n', ...
        nReseedWarmR);
elseif okNoReseed && ~okReseed
    warning(['FAIL: method 4 (reseedEmpty=%d) emitted no "Reseeded" messages. ', ...
        'Either no unused factors met reseedThresh, or reseed did not run. ', ...
        'SeqNMF near-empty at init was %d / %d.'], ...
        reseedEmpty, nEmpty0, Khat);
else
    warning('FAIL: unexpected reseed messages with reseedEmpty=0 (rand=%d, warm=%d).', ...
        nReseedRand, nReseedWarm0);
end

if exist(diaryFile, 'file')
    delete(diaryFile);
end

%% ------------------------------------------------------------------------
function n = count_reseed_msgs(diaryFile)
    n = 0;
    if ~exist(diaryFile, 'file')
        return
    end
    txt = fileread(diaryFile);
    n = numel(regexp(txt, 'Reseeded \d+ unused factors at iter'));
end

function nEmpty = count_near_empty_factors(W, H, X)
    % Match FlexMF reseed criterion: ||W_k (*) H_k|| vs max(factor recon, ||X||)
    thresh = 1e-3;  % matches FlexMF reseedThresh default
    K = size(W, 2);
    reconNorms = zeros(1, K);
    for kk = 1:K
        Xk = helper.reconstruct(W(:, kk, :), H(kk, :));
        reconNorms(kk) = norm(Xk(:));
    end
    nEmpty = nnz(reconNorms < thresh * max([max(reconNorms), norm(X(:))]));
end
