clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))

data_types = {'warp', 'jitter', 'participation', 'noise', 'warpnoise', 'jitternoise'};

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

%% Generate some synthetic data
id = str2double(getenv("SLURM_ARRAY_TASK_ID"));
data_type = data_types{id}

K = 3;
Khat = 5;
T = 4000; % length of data to generate
Nneurons = 10*ones(K,1); % number of neurons in each sequence
Dt = 3.*ones(K,1); % gap between each member of the sequence
neg = 0;
noise_levels = 0:.001:.01;
participation_levels = 1:-0.1:.1; 
jitter_levels = 0:9;
warp_levels = 0:9;
gap = 100;

switch data_type
    case {'warp', 'warpnoise'}
        nLevels = numel(warp_levels);
    case {'jitter', 'jitternoise'}
        nLevels = numel(jitter_levels);
    case 'participation'
        nLevels = numel(participation_levels);
    case 'noise'
        nLevels = numel(noise_levels);
end

nSim = 10;
rng(1)
seeds = randperm(1000, nSim);
nTotal = nLevels*nSim;

emds_W_SeqNMF_flat = cell(nTotal,1);
emds_H_SeqNMF_flat = cell(nTotal,1);
emds_W_FlexMF_flat = cell(nTotal,1);
emds_H_FlexMF_flat = cell(nTotal,1);
ids_SeqNMF_flat = cell(nTotal,1);
ids_FlexMF_flat = cell(nTotal,1);
pvals_SeqNMF_flat = cell(nTotal,1);
pvals_FlexMF_flat = cell(nTotal,1);
is_significant_SeqNMF_flat = cell(nTotal,1);
is_significant_FlexMF_flat = cell(nTotal,1);
Ws_flat = cell(nTotal,1);
Hs_train_flat = cell(nTotal,1);
Whats_SeqNMF_flat = cell(nTotal,1);
Hhats_train_SeqNMF_flat = cell(nTotal,1);
Whats_FlexMF_flat = cell(nTotal,1);
Hhats_train_FlexMF_flat = cell(nTotal,1);
Xs_train_flat = cell(nTotal,1);
Xs_test_flat = cell(nTotal,1);
times_train = zeros(nTotal,1);
times_match = zeros(nTotal,1);
times_test = zeros(nTotal,1);

parfor g = 1:nTotal
    [l, n] = ind2sub([nLevels, nSim], g);
    display(['n=' num2str(n) ', level=' num2str(l)])
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
        X_train = X_train/frob_norm*Khat;
        W = W/frob_norm*Khat;
        frob_norm = norm(X_test(:));
        X_test = X_test/frob_norm*Khat;

        % Run SeqNMF
        lambda_SeqNMF = .05;
        [What_SeqNMF, Hhat_train_SeqNMF]= seqNMF(X_train,'K',Khat,'L',L,...
                'lambda', lambda_SeqNMF, 'maxiter', 50, 'showPlot', 0); 

        % Run FlexMF
        lambda_FlexMF = 1;
        lambda_M = .1;
        lambda_R = 1;
        tic
        [What_FlexMF, Hhat_train_FlexMF, ~, ~, loadings, power, M_train, R_train] = FlexMF(X_train, 'K', Khat, 'L', L, ...
        'EMD',1, 'lambda', lambda_FlexMF, 'lambda_R', lambda_R, 'lambda_M', lambda_M, ...
        'W_init', What_SeqNMF, 'H_init', Hhat_train_SeqNMF, ...
        'maxiter', 50, 'showPlot', 0, 'verbal', 0);
        times_train(g) = toc;

        % Compare algorithms
        disp('Evaluate EMDs of results')
        tic
        [emds_W_SeqNMF_flat{g}, emds_H_SeqNMF_flat{g}, ids_SeqNMF_flat{g}] = helper.similarity_WH_EMD(W, H_train, What_SeqNMF, Hhat_train_SeqNMF);
        [emds_W_FlexMF_flat{g}, emds_H_FlexMF_flat{g}, ids_FlexMF_flat{g}] = helper.similarity_WH_EMD(W, H_train, What_FlexMF, Hhat_train_FlexMF);
        times_match(g) = toc;
        
        % Test significance
        disp('Test Significance')
        [pvals_SeqNMF_flat{g},is_significant_SeqNMF_flat{g}] = test_significance(X_test, What_SeqNMF);
        tic
        [What_FlexMF, Hhat_test_FlexMF, cost_test, errors_test, ~, ~, M_test, R_test] = FlexMF(X_test, 'K', Khat, 'L', L, 'W_fixed', 1, 'W_init', What_FlexMF,...
            'EMD',1, 'lambda', lambda_FlexMF, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'maxiter', 50, 'showPlot', 0, 'verbal', 0);
        times_test(g) = toc;
        [pvals_FlexMF_flat{g},is_significant_FlexMF_flat{g},~] = test_significance_EMD(X_test, What_FlexMF, M_test, 'plot', 0);
        
        Ws_flat{g} = W;
        Hs_train_flat{g} = H_train;
        Whats_SeqNMF_flat{g} = What_SeqNMF;
        Hhats_train_SeqNMF_flat{g} = Hhat_train_SeqNMF;
        Whats_FlexMF_flat{g} = What_FlexMF;
        Hhats_train_FlexMF_flat{g} = Hhat_train_FlexMF;
        Xs_train_flat{g} = X_train;
        Xs_test_flat{g} = X_test;
end

emds_W_SeqNMF = reshape(emds_W_SeqNMF_flat, [nLevels, nSim]);
emds_H_SeqNMF = reshape(emds_H_SeqNMF_flat, [nLevels, nSim]);
emds_W_FlexMF = reshape(emds_W_FlexMF_flat, [nLevels, nSim]);
emds_H_FlexMF = reshape(emds_H_FlexMF_flat, [nLevels, nSim]);
ids_SeqNMF = reshape(ids_SeqNMF_flat, [nLevels, nSim]);
ids_FlexMF = reshape(ids_FlexMF_flat, [nLevels, nSim]);
pvals_SeqNMF = reshape(pvals_SeqNMF_flat, [nLevels, nSim]);
pvals_FlexMF = reshape(pvals_FlexMF_flat, [nLevels, nSim]);
is_significant_SeqNMF = reshape(is_significant_SeqNMF_flat, [nLevels, nSim]);
is_significant_FlexMF = reshape(is_significant_FlexMF_flat, [nLevels, nSim]);
Ws = reshape(Ws_flat, [nLevels, nSim]);
Hs_train = reshape(Hs_train_flat, [nLevels, nSim]);
Whats_SeqNMF = reshape(Whats_SeqNMF_flat, [nLevels, nSim]);
Hhats_train_SeqNMF = reshape(Hhats_train_SeqNMF_flat, [nLevels, nSim]);
Whats_FlexMF = reshape(Whats_FlexMF_flat, [nLevels, nSim]);
Hhats_train_FlexMF = reshape(Hhats_train_FlexMF_flat, [nLevels, nSim]);
Xs_train = reshape(Xs_train_flat, [nLevels, nSim]);
Xs_test = reshape(Xs_test_flat, [nLevels, nSim]);
times_train = reshape(times_train, [nLevels, nSim]);
times_match = reshape(times_match, [nLevels, nSim]);
times_test = reshape(times_test, [nLevels, nSim]);

save(fullfile('Simulation_Results', sprintf('Compare_FlexMF_SeqNMF_%s.mat', data_type)), ...
"emds_W_FlexMF", "emds_H_FlexMF", "emds_H_SeqNMF", "emds_W_SeqNMF", "ids_SeqNMF", "ids_FlexMF",...,
"pvals_SeqNMF", "pvals_FlexMF", "is_significant_SeqNMF", "is_significant_FlexMF", "Ws", "Hs_train",...
"Whats_SeqNMF", "Hhats_train_SeqNMF", "Whats_FlexMF", "Hhats_train_FlexMF", "Xs_train", "Xs_test",...
"times_train", "times_match", "times_test")

% Shut down the parallel pool
delete(gcp('nocreate'));