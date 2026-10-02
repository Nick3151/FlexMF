clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))

data_types = {'warpnoise', 'jitternoise'};

% Start a local pool inside the CPUs allocated by the SLURM job.
nWorkers = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isnan(nWorkers) || nWorkers < 1
    nWorkers = 1;
end
pool = gcp('nocreate');
if isempty(pool)
    fprintf('Starting local parallel pool (%d workers)...\n', nWorkers);
    parpool('local', nWorkers);
else
    fprintf('Using existing parallel pool (%d workers).\n', pool.NumWorkers);
end

%% Generate some synthetic data, with warping
id = str2double(getenv("SLURM_ARRAY_TASK_ID"));
data_type = data_types{id}

nSim = 10;
rng(1)
seeds = randperm(1000, nSim);

K = 3;
Khat = 5;
T = 4000; % length of data to generate
Nneurons = 10*ones(K,1); % number of neurons in each sequence
Dt = 3.*ones(K,1); % gap between each member of the sequence
neg = 0;
noise = .005;
participation = .8.*ones(K,1); 
jitter = 5*ones(K,1);
warp = 5;
gap = 100;

%% Run simulation on different combinations of lambda_M, lambda_R
lambda_M = 1e-1;
nlambdas = 17;

lambda_R = 1;
lambda_SeqNMF = .05;
lambdas = logspace(-2, 2, nlambdas);

nTotal = nlambdas*nSim;

times_flat = cell(nTotal,1);
emds_W_flat = cell(nTotal,1);
emds_H_flat = cell(nTotal,1);
recon_costs_flat = zeros(nTotal,1);
reg_costs_flat = zeros(nTotal,1);
constraints_flat = zeros(nTotal,1);
num_detected_flat = cell(nTotal,1);
num_significant_flat = cell(nTotal,1);
is_significant_flat = cell(nTotal,1);
W_hats_flat = cell(nTotal,1);
H_hats_flat = cell(nTotal,1);
Ws_flat = cell(nTotal,1);
Hs_flat = cell(nTotal,1);
Xs_flat = cell(nTotal,1);
Rs_flat = cell(nTotal,1);
Ms_flat = cell(nTotal,1);
ids_match_flat = cell(nTotal,1);

parfor g = 1:nTotal
    [Li, n] = ind2sub([nlambdas, nSim], g);
    lambda = lambdas(Li);
    display(['n=' num2str(n) ', lambda=' num2str(lambda)])
    switch data_type
        case 'warpnoise'
            [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',noise, 'warp', warp, 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
        case 'jitternoise'
            [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',noise, 'jitter', jitter, 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
    end
    Xtrain = X(:,1:round(T/2));
    Xtest = X(:,1+round(T/2):end);
    Htrain = H(:,1:round(T/2));
    L = size(W,3);

    % Normalize training and test data
    frob_norm = norm(Xtrain(:));
    Xtrain = Xtrain/frob_norm*Khat;
    W = W/frob_norm*Khat;
    frob_norm = norm(Xtest(:));
    Xtest = Xtest/frob_norm*Khat;

    % Warm-start FlexMF with SeqNMF results
    [What_SeqNMF, Hhat_SeqNMF] = seqNMF(Xtrain,'K',Khat,'L',L,...
        'lambda', lambda_SeqNMF, 'maxiter', 50, 'showPlot', 0);

    tic
    [W_hat, H_hat, cost, errors, loadings, power, M, R]= FlexMF(Xtrain,'K',Khat,'L',L, 'EMD',1, 'maxiter', 50,...
        'lambda', lambda, 'lambda_M', lambda_M, 'lambda_R', lambda_R, ...
        'W_init', What_SeqNMF, 'H_init', Hhat_SeqNMF, 'showPlot', 0, 'verbal', 0);
    [recon_err, reg_cross, reg_W, reg_H] = helper.get_FlexMF_cost(Xtrain,W_hat,H_hat);
    Xhat = helper.reconstruct(W_hat, H_hat);

    Xtrain_corr = helper.correct_warp(Xtrain, M);
    constraint = Xtrain_corr-R-Xhat;
    recon_costs_flat(g) = sum((Xtrain_corr(:)-Xhat(:)).^2)/2;
    reg_costs_flat(g) = reg_cross;
    constraints_flat(g) = norm(constraint(:),1);

    times_flat{g} = toc;
    W_hats_flat{g} = W_hat;
    H_hats_flat{g} = H_hat;
    Ms_flat{g} = M;
    Rs_flat{g} = R;
    Ws_flat{g} = W;
    Hs_flat{g} = H;
    Xs_flat{g} = X;

    disp('Evaluate EMDs of W and H to ground truth')
    [emds_W_flat{g}, emds_H_flat{g}, ids] = helper.similarity_WH_EMD(W, Htrain, W_hat, H_hat);
    num_detected_flat{g} = length(ids);
    ids_match_flat{g} = ids;

    disp('Test Significance')
    [W_hat, Hhat_test, cost_test, errors_test, ~, ~, M_test, R_test] = FlexMF(Xtest, 'K', Khat, 'L', L, 'W_fixed', 1, 'W_init', W_hat,...
    'EMD',1, 'lambda', lambda, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'maxiter', 50, 'showPlot', 0, 'verbal', 0);
    [pvals,is_significant,is_single] = test_significance_EMD(Xtest, W_hat, M_test, 'plot', 0);
    num_significant_flat{g} = sum(is_significant);
    is_significant_flat{g} = is_significant;
end

times = reshape(times_flat, [nlambdas, nSim]);
emds_W = reshape(emds_W_flat, [nlambdas, nSim]);
emds_H = reshape(emds_H_flat, [nlambdas, nSim]);
recon_costs = reshape(recon_costs_flat, [nlambdas, nSim]);
reg_costs = reshape(reg_costs_flat, [nlambdas, nSim]);
constraints = reshape(constraints_flat, [nlambdas, nSim]);
num_detected = reshape(num_detected_flat, [nlambdas, nSim]);
num_significant = reshape(num_significant_flat, [nlambdas, nSim]);
is_significant = reshape(is_significant_flat, [nlambdas, nSim]);
W_hats = reshape(W_hats_flat, [nlambdas, nSim]);
H_hats = reshape(H_hats_flat, [nlambdas, nSim]);
Ws = reshape(Ws_flat, [nlambdas, nSim]);
Hs = reshape(Hs_flat, [nlambdas, nSim]);
Xs = reshape(Xs_flat, [nlambdas, nSim]);
Rs = reshape(Rs_flat, [nlambdas, nSim]);
Ms = reshape(Ms_flat, [nlambdas, nSim]);
ids_match = reshape(ids_match_flat, [nlambdas, nSim]);

save(fullfile('Simulation_Results', sprintf('FlexMF_choose_lambda_%s_lambdaM=%0.3e.mat', data_type, lambda_M)), ...
    'times', 'lambdas', 'lambda_R', 'lambda_M', 'W_hats', 'H_hats', 'Ws', 'Hs', 'Xs', 'Ms', 'Rs', 'ids_match', ...
    'emds_W', 'emds_H', 'recon_costs', 'reg_costs', 'constraints', 'num_detected', 'num_significant', 'is_significant')

% Shut down the parallel pool
delete(gcp('nocreate'));
