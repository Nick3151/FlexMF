clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))

data_types = {'warp', 'jitter', 'participation', 'noise'};

% Start parallel pool
disp('Starting a parallel pool...');
% parpool(rf, 24);  % 2 nodes × 12 tasks
parpool('local', 32);

%% Generate some synthetic data, with warping
id = str2double(getenv("SLURM_ARRAY_TASK_ID"));
data_type = data_types{id}

nSim = 10;
rng(1)
seeds = randperm(1000, nSim);

K = 3;
T = 4000; % length of data to generate
Nneurons = 5*ones(K,1); % number of neurons in each sequence
Dt = 3.*ones(K,1); % gap between each member of the sequence
neg = 0;
noise = .01;
participation = .7.*ones(K,1); 
jitter = 5*ones(K,1);
warp = 2;
gap = 100;

%% Run simulation on different combinations of lambda_M, lambda_R
lambda = 1e-2;
lambda_M = 1e-2;
lambda_R = 1;
Ks = 1:10;
nKs = length(Ks);

times = cell(nKs, nSim);
emds_W = cell(nKs, nSim);
emds_H = cell(nKs, nSim);
emds_X = zeros(nKs, nSim);   % EMDs between X and Xhat
reg_costs = zeros(nKs, nSim);
constraints = zeros(nKs, nSim);
num_detected = cell(nKs, nSim);
num_significant = cell(nKs, nSim);
W_hats = cell(nKs, nSim);
H_hats = cell(nKs, nSim);
Ws = cell(nSim, 1);
Hs = cell(nSim, 1);
Xs = cell(nSim, 1);
Rs = cell(nKs, nSim);
Ms = cell(nKs, nSim);
ids_match = cell(nKs, nSim);

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
    end
    Xtrain = X(:,1:round(T/2));
    Xtest = X(:,1+round(T/2):end);
    Htrain = H(:,1:round(T/2));

    Ws{n} = W;
    Hs{n} = H;
    Xs{n} = X;
    L = size(W,3);
    parfor Ki=1:nKs
        K_tmp = Ks(Ki);
        display(['Testing K=' num2str(K_tmp)])
        tic
        [W_hat, H_hat, cost, errors, loadings, power, M, R]= FlexMF(Xtrain,'K',K_tmp,'L',L, 'EMD',1, 'maxiter', 50,...
            'lambda', lambda, 'lambda_M', lambda_M, 'lambda_R', lambda_R, 'showPlot', 0, 'verbal', 0); 
        [recon_err, reg_cross, reg_W, reg_H] = helper.get_FlexMF_cost(Xtrain,W_hat,H_hat);
        Xhat = helper.reconstruct(W_hat, H_hat);
        Ttrain = size(Xtrain,2);
        D = eye(Ttrain) - diag(ones(Ttrain-1,1),-1);
        D(Ttrain,Ttrain) = 0;
        constraint = M*D'-R-Xhat+Xtrain;
        disp('Evaluate EMDs between Xtrain and Xhat')
        emds_X(Ki,n) = compute_EMD(Xtrain,Xhat,opts, 'continuationOptions', continue_opts, 'lambdaR', 1e2);
        reg_costs(Ki,n) = reg_cross;
        constraints(Ki,n) = norm(constraint(:),1);

        time = toc
        times{Ki,n} = time;
        W_hats{Ki,n} = W_hat;
        H_hats{Ki,n} = H_hat;
        Ms{Ki,n} = M;
        Rs{Ki,n} = R;
        
        disp('Evaluate EMDs of W and H to ground truth')
        tic
        [emds_W{Ki,n}, emds_H{Ki,n}, ids] = helper.similarity_WH_EMD(W, Htrain, W_hat, H_hat);
        num_detected{Ki,n} = length(ids);
        ids_match{Ki,n} = ids;
        toc

        disp('Fit on Test data')
        [W_hat, Hhat_test, cost_test, errors_test, ~, ~, M_test, R_test] = FlexMF(Xtest, 'K', K_tmp, 'L', L, 'W_fixed', 1, 'W_init', W_hat,...
        'EMD',1, 'lambda', lambda, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'maxiter', 1, 'showPlot', 0, 'verbal', 0);
        [pvals,is_significant,is_single] = test_significance_EMD(Xtest, W_hat, M_test, 'plot', 0);
        num_significant{Ki,n} = sum(is_significant);
    end
end

save(fullfile('Simulation_Results', sprintf('FlexMF_choose_lambda_%s_lambdaM=%0.3e.mat', data_type, lambda_M)), ...
    'times', 'Ks', 'lambda_R', 'lambda_M', 'lambda', 'W_hats', 'H_hats', 'Ws', 'Hs', 'Xs', 'Ms', 'Rs', 'ids_match', ...
    'emds_W', 'emds_H', 'emds_X', 'reg_costs', 'constraints', 'num_detected', 'num_significant')

% Shut down the parallel pool
delete(gcp('nocreate'));
