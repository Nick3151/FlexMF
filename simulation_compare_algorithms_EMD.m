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

%% Generate some synthetic data
id = str2double(getenv("SLURM_ARRAY_TASK_ID"));
data_type = data_types{id}

K = 3;
T = 4000; % length of data to generate
Nneurons = 5*ones(K,1); % number of neurons in each sequence
Dt = 5.*ones(K,1); % gap between each member of the sequence
neg = 0;
noise_levels = 0:.0005:.002;
participation_levels = 1:-0.1:.6; 
jitter_levels = 0:4;
warp_levels = 0:4;
gap = 100;

nSim = 10;
rng(1)
seeds = randperm(1000, nSim);

emds_W_SeqNMF = cell(5,nSim);
emds_H_SeqNMF = cell(5,nSim);
emds_W_FlexMF = cell(5,nSim);
emds_H_FlexMF = cell(5,nSim);
ids_SeqNMF = cell(5,nSim);
ids_FlexMF = cell(5,nSim);
pvals_SeqNMF = cell(5,nSim);
pvals_FlexMF = cell(5,nSim);
is_significant_SeqNMF = cell(5,nSim);
is_significant_FlexMF = cell(5,nSim);
Ws = cell(5,nSim);
Hs_train = cell(5,nSim);
Whats_SeqNMF = cell(5,nSim);
Hhats_train_SeqNMF = cell(5,nSim);
Whats_FlexMF = cell(5,nSim);
Hhats_train_FlexMF = cell(5,nSim);
Xs_train = cell(5,nSim);
Xs_test = cell(5,nSim);

parfor n = 1:nSim
    display(['n=' num2str(n)])
    for l=1:5
        switch data_type
            case 'warp'
                [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',0, 'warp', warp_levels(l), 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
                display(['warp level = ' num2str(warp_levels(l))])
            case 'jitter'
                [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',0, 'jitter', jitter_levels(l).*ones(K,1), 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
                display(['jitter level = ' num2str(jitter_levels(l))])
            case 'participation'
                [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',0, 'participation', participation_levels(l).*ones(K,1), 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
                display(['participation level = ' num2str(participation_levels(l))])
            case 'noise'
                [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise', noise_levels(l), 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
                display(['noise level = ' num2str(noise_levels(l))])
            case 'warpnoise'
                [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',noise_levels(3), 'warp', warp_levels(l), 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
                display(['warp level = ' num2str(warp_levels(l))])
            case 'jitternoise'
                [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',noise_levels(3), 'jitter', jitter_levels(l).*ones(K,1), 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
                display(['jitter level = ' num2str(jitter_levels(l))])
        end    
    
        L = size(W,3);
        X_train = X(:,1:round(T/2));
        X_test = X(:,1+round(T/2):end);
        H_train = H(:,1:round(T/2));
        H_test = H(:,1+round(T/2):end);
        
        % Split into training and test set, normalize data
        frob_norm = norm(X_train(:));
        X_train = X_train/frob_norm*K;
        W = W/frob_norm*K;
        frob_norm = norm(X_test(:));
        X_test = X_test/frob_norm*K;

        % Run SeqNMF
        lambda = .05;
        [What_SeqNMF, Hhat_train_SeqNMF]= seqNMF(X_train,'K',K,'L',L,...
                    'lambda', lambda, 'maxiter', 50, 'showPlot', 0); 

        % Run FlexMF
        lambda = .05;
        lambda_M = .01;
        lambda_R = 1;
        [What_FlexMF, Hhat_train_FlexMF, ~, ~, loadings, power, M_train, R_train] = FlexMF(X_train, 'K', K, 'L', L, ...
        'EMD',1, 'lambda', lambda, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'maxiter', 50, 'showPlot', 0, 'verbal', 0);

        % Compare algorithms
        disp('Evaluate EMDs of results')
        tic
        [emds_W_SeqNMF{l,n}, emds_H_SeqNMF{l,n}, ids_SeqNMF{l,n}] = helper.similarity_WH_EMD(W, H_train, What_SeqNMF, Hhat_train_SeqNMF);
        [emds_W_FlexMF{l,n}, emds_H_FlexMF{l,n}, ids_FlexMF{l,n}] = helper.similarity_WH_EMD(W, H_train, What_FlexMF, Hhat_train_FlexMF);
        toc
        
        % Test significance
        disp('Test Significance')
        [pvals_SeqNMF{l,n},is_significant_SeqNMF{l,n}] = test_significance(X_test, What_SeqNMF);
        [What_FlexMF, Hhat_test_FlexMF, cost_test, errors_test, ~, ~, M_test, R_test] = FlexMF(X_test, 'K', K, 'L', L, 'W_fixed', 1, 'W_init', What_FlexMF,...
            'EMD',1, 'lambda', lambda, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'maxiter', 1, 'showPlot', 0, 'verbal', 0);
        [pvals_FlexMF{l,n},is_significant_FlexMF{l,n},~] = test_significance_EMD(X_test, What_FlexMF, M_test, 'plot', 0);
        
        Ws{l,n} = W;
        Hs_train{l,n} = H_train;
        Whats_SeqNMF{l,n} = What_SeqNMF;
        Hhats_train_SeqNMF{l,n} = Hhat_train_SeqNMF;
        Whats_FlexMF{l,n} = What_FlexMF;
        Hhats_train_FlexMF{l,n} = Hhat_train_FlexMF;
        Xs_train{l,n} = X_train;
        Xs_test{l,n} = X_test;
    end
end

save(fullfile('Simulation_Results', sprintf('Compare_FlexMF_SeqNMF_%s.mat', data_type)), ...
"emds_W_FlexMF", "emds_H_FlexMF", "emds_H_SeqNMF", "emds_W_SeqNMF", "ids_SeqNMF", "ids_FlexMF",...,
"pvals_SeqNMF", "pvals_FlexMF", "is_significant_SeqNMF", "is_significant_FlexMF", "Ws", "Hs_train",...
"Whats_SeqNMF", "Hhats_train_SeqNMF", "Whats_FlexMF", "Hhats_train_FlexMF", "Xs_train", "Xs_test")

% Shut down the parallel pool
delete(gcp('nocreate'));