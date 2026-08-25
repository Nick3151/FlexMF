clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))

data_types = {'warp', 'jitter', 'participation', 'noise', 'warpnoise', 'jitternoise'};

% Load the cluster profile
rf = parcluster('rockfish');

% Set SLURM resource parameters
rf.AdditionalProperties.Partition = 'parallel';  % Specify partition
rf.AdditionalProperties.WallTime = '72:00:00';   % Wall time must match SLURM script
rf.AdditionalProperties.AdditionalSubmitArgs = '--nodes=1 --ntasks-per-node=32 --cpus-per-task=1';

% Display properties (optional)
disp(rf.AdditionalProperties);

% Start parallel pool
disp('Starting a parallel pool...');
parpool(rf, 32);  % 1 nodes × 32 tasks

% parpool('local', 32);

%% Generate some synthetic data, with warping
Li = str2double(getenv("SLURM_ARRAY_TASK_ID"));
data_type = data_types{Li}

nSim = 10;
rng(1)
seeds = randperm(1000, nSim);

K = 3;
T = 4000; % length of data to generate
Nneurons = 5*ones(K,1); % number of neurons in each sequence
Dt = 3.*ones(K,1); % gap between each member of the sequence
neg = 0;
noise = .001;
participation = .8.*ones(K,1); 
jitter = 5*ones(K,1);
warp = 2;
gap = 100;

%% Run simulation on different combinations of lambda_M, lambda_R
lambda = 1e-1;
nMs = 9;

lambda_R = 1;
% lambdas = logspace(-6,-2,nLambdas);
lambda_Ms = logspace(-4,0,nMs);

times = cell(nMs, nSim);
emds_W = cell(nMs, nSim);
emds_H = cell(nMs, nSim);
recon_costs = zeros(nMs, nSim);   
reg_costs = zeros(nMs, nSim);
constraints = zeros(nMs, nSim);
num_detected = cell(nMs, nSim);
num_significant = cell(nMs, nSim);
W_hats = cell(nMs, nSim);
H_hats = cell(nMs, nSim);
Ws = cell(nSim, 1);
Hs = cell(nSim, 1);
Xs = cell(nSim, 1);
Rs = cell(nMs, nSim);
Ms = cell(nMs, nSim);
ids_match = cell(nMs, nSim);

opts = tfocs_SCD;
opts.continuation = 1;
opts.tol = 1e-4;
opts.stopCrit = 4;
opts.maxIts = 500;
opts.printEvery = 0;
opts.alg = 'N83';
continue_opts = continuation();
continue_opts.verbose = 0;

for n=1:nSim
    display(['n=' num2str(n)])
    switch data_type
        case 'warp'
            [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',0, 'warp', warp, 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
        case 'jitter'
            [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',0, 'jitter', jitter, 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
        case 'participation'
            [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',0, 'participation', participation, 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
        case 'noise'
            [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise', noise, 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
        case 'warpnoise'
            [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',noise, 'warp', warp, 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
        case 'jitternoise'
            [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',noise, 'jitter', jitter, 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
    end
    Xtrain = X(:,1:round(T/2));
    Xtest = X(:,1+round(T/2):end);
    Htrain = H(:,1:round(T/2));

    Ws{n} = W;
    Hs{n} = H;
    Xs{n} = X;
    L = size(W,3);
    parfor Mi=1:nMs
        lambda_M = lambda_Ms(Mi);
        display(['Testing lambda_M=' num2str(lambda_M)])
        tic
        [W_hat, H_hat, cost, errors, loadings, power, M, R]= FlexMF(Xtrain,'K',K,'L',L, 'EMD',1, 'maxiter', 50,...
            'lambda', lambda, 'lambda_M', lambda_M, 'lambda_R', lambda_R, 'showPlot', 0, 'verbal', 0); 
        [recon_err, reg_cross, reg_W, reg_H] = helper.get_FlexMF_cost(Xtrain,W_hat,H_hat);
        Xhat = helper.reconstruct(W_hat, H_hat);
        Xtrain_corr = helper.correct_warp(Xtrain, M);
        constraint = Xtrain_corr-R-Xhat;
        recon_costs(Mi,n) = sum((Xtrain_corr(:)-Xhat(:)).^2)/2;
        reg_costs(Mi,n) = reg_cross;
        constraints(Mi,n) = norm(constraint(:),1);

        time = toc
        times{Mi,n} = time;
        W_hats{Mi,n} = W_hat;
        H_hats{Mi,n} = H_hat;
        Ms{Mi,n} = M;
        Rs{Mi,n} = R;
        
        disp('Evaluate EMDs of W and H to ground truth')
        tic
        [emds_W{Mi,n}, emds_H{Mi,n}, ids] = helper.similarity_WH_EMD(W, Htrain, W_hat, H_hat);
        num_detected{Mi,n} = length(ids);
        ids_match{Mi,n} = ids;
        toc

        disp('Fit on Test data')
        [W_hat, Hhat_test, cost_test, errors_test, ~, ~, M_test, R_test] = FlexMF(Xtest, 'K', K, 'L', L, 'W_fixed', 1, 'W_init', W_hat,...
        'EMD',1, 'lambda', lambda, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'maxiter', 1, 'showPlot', 0, 'verbal', 0);
        [pvals,is_significant,is_single] = test_significance_EMD(Xtest, W_hat, M_test, 'plot', 0);
        num_significant{Mi,n} = sum(is_significant);
    end
end

save(fullfile('Simulation_Results', sprintf('FlexMF_EMD_results_%s_lambda=%0.3e.mat', data_type, lambda)), ...
    'times', 'lambda_Ms', 'lambda_R', 'lambda', 'W_hats', 'H_hats', 'Ws', 'Hs', 'Xs', 'Ms', 'Rs', 'ids_match', ...
    'emds_W', 'emds_H', 'recon_costs', 'reg_costs', 'constraints', 'num_detected', 'num_significant')

% Shut down the parallel pool
delete(gcp('nocreate'));
