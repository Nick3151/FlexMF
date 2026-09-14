% Simulate data with trials, compare the performance of SeqNMF with
% multiplication update rule and Split Bregman Iteration

clear all
close all
clc
root = fileparts(pwd);
tfocs_dir  = fullfile(root, 'TFOCS');
utils_dir  = fullfile(root, 'Utils');
flexmf_dir = fullfile(root, 'FlexMF');
addpath(tfocs_dir)
addpath(utils_dir)
addpath(genpath(flexmf_dir));

% Load the cluster profile
delete(gcp('nocreate'));
rf = parcluster('rockfish');

% Set SLURM resource parameters
rf.AdditionalProperties.Partition = 'parallel';  % Specify partition
rf.AdditionalProperties.WallTime = '72:00:00';   % Wall time must match SLURM script
rf.AdditionalProperties.MemPerCPU = '3G';
rf.NumThreads = 1;                               % CPUs (threads) per worker
rf.AdditionalProperties.AdditionalSubmitArgs = '--nodes=1 --ntasks-per-node=48';

% Display properties (optional)
disp(rf.AdditionalProperties);

% Start parallel pool
% SINGLE NODE: 48 workers x 1 thread = 48 cores (one full 'parallel' node).
% Multi-node communicating pools do NOT start on this cluster (the MPI ring
% fails to form across nodes; >48 workers was confirmed to fail while 32 on a
% single node works). 48-way parallelism over the flattened loops is still a
% large speedup. If 48 ever fails to start, fall back to 32:
%   rf.AdditionalProperties.AdditionalSubmitArgs = '--nodes=1 --ntasks-per-node=32';
%   parpool(rf, 32);
disp('Starting a parallel pool...');
parpool(rf, 48);

% The addpath calls above only affect the client. The workers run in fresh
% MATLAB sessions, so the user code must be added to THEIR path too. Rockfish
% has a shared filesystem, so the workers can use the same absolute paths.
pctRunOnAll(['addpath(''' tfocs_dir ''')']);
pctRunOnAll(['addpath(''' utils_dir ''')']);
pctRunOnAll(['addpath(genpath(''' flexmf_dir '''))']);

n = str2double(getenv("SLURM_ARRAY_TASK_ID"));

%% Impact of motif shape
if n==1
    disp('Impact of motif shape on results')
    Trials = 200;
    L = 50; % length of each trial
    K = 10;
    Nmotifs = 2*(1:K);
    Nneurons = 5*ones(K, 1); % the number of neurons in each motif
    Dt = 3.*ones(K,1); % gap between each member of the motif
    
    nSim = 50;
    seeds = randperm(1000, nSim);
    pvals_SeqNMF = cell(nSim,3);
    is_significants_SeqNMF = cell(nSim,3);
    loadings_SeqNMF = cell(nSim,3);
    W_hats_SeqNMF = cell(nSim,3);
    H_hats_SeqNMF = cell(nSim,3);
    pvals_FlexMF = pvals_SeqNMF;
    is_significants_FlexMF = is_significants_SeqNMF;
    loadings_FlexMF = loadings_SeqNMF;
    W_hats_FlexMF = W_hats_SeqNMF;
    H_hats_FlexMF = H_hats_SeqNMF;
    Ws = W_hats_SeqNMF;
    Hs = H_hats_SeqNMF;
    TrainingDatas = W_hats_SeqNMF;
    
    % Flattened: one parfor over all nSim x nCond conditions. Column c selects
    % the condition (1=Calcium Transients, 2=Burst, 3=Single spikes). idx maps
    % to (i,c) of the nSim-by-3 cells via column-major ind2sub, so results match
    % the original three separate parfor blocks.
    nCond = 3;
    len_bursts = [10, 10, 1];        % continuous firing time per condition
    dynamics   = [1, 0, 0];          % consider calcium dynamics or not
    lambdas    = [.005, .005, .01];
    alpha_Hs   = [1e-2, 1e-2, 1e-3];
    labels     = {'Calcium Transients:', 'Burst:', 'Single spikes:'};
    parfor idx=1:nSim*nCond
        [i, c] = ind2sub([nSim, nCond], idx);
        tic
        disp(labels{c})
        len_burst = len_bursts(c); % Continuous firing time
        dynamic = dynamics(c); % Consider calcium dynamics or not
        [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, ...
            'len_burst', len_burst, 'dynamic', dynamic, 'seed', seeds(i));
        groups = zeros(Trials,1);
        for k=1:K
            groups(motif_ind{k}) = k;
        end
        Ws{idx} = W;
        
        % Dimension N*L*Trials
        cv = cvpartition(groups, "KFold",2);
        ind_train = find(training(cv,1));
        X_train = X(:,:,ind_train);
        ind_test = find(test(cv,1));
        X_test = X(:,:,ind_test);
        
        % Dimension N*T
        N = size(W,1);
        TrainingData = zeros(N,cv.TrainSize(1)*L);
        for t=1:cv.TrainSize(1)
            TrainingData(:,(t-1)*L+1:t*L) = squeeze(X_train(:,:,t));
        end
        TestData = zeros(N,cv.TestSize(1)*L);
        for t=1:cv.TestSize(1)
            TestData(:,(t-1)*L+1:t*L) = squeeze(X_test(:,:,t));
        end
        
        % Normalize training data
        frob_norm = norm(TrainingData(:));
        TrainingData = TrainingData/frob_norm*K;

        H_train = H(:,ind_train);
        T = cv.TrainSize(1)*L;
        H_train_full = zeros(K,T);
        H_train_full(:,1:L:T) = H_train;
        Hs{idx} = H_train_full;
    
        lambda = lambdas(c);
    
        % Run SeqNMF with multiplication rule
        [W_hat, H_hat, ~,~,loadings_SeqNMF{idx},power]= seqNMF(TrainingData,'K',K,'L',L,...
                'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
        p = .05;
        [pvals_SeqNMF{idx},is_significants_SeqNMF{idx}] = test_significance_new(TestData, W_hat,[],p);
        W_hats_SeqNMF{idx} = W_hat;
        H_hats_SeqNMF{idx} = H_hat;
        display(['SeqNMF run ' num2str(i) '/' num2str(nSim)])
    
        % Run SeqNMF with Bregman Iteration
        alpha_W = 1e-6;
        alpha_H = alpha_Hs(c);
        [W_hat, H_hat, ~,~,loadings_FlexMF{idx},power]= FlexMF(TrainingData,'K',K,'L',L,...
                'lambda', lambda, 'alpha_W', alpha_W, 'alpha_H', alpha_H, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0, 'verbal', 0); 
        p = .05;
        [pvals_FlexMF{idx},is_significants_FlexMF{idx}] = test_significance_new(TestData, W_hat,[],p);
        W_hats_FlexMF{idx} = W_hat;
        H_hats_FlexMF{idx} = H_hat;
        TrainingDatas{idx} = TrainingData;
        display(['FlexMF run ' num2str(i) '/' num2str(nSim)])
        toc
    end
    
    save('simulate_results_shape.mat', "H_hats_FlexMF", 'W_hats_FlexMF', 'is_significants_FlexMF', 'pvals_FlexMF', 'loadings_FlexMF',...
        'loadings_SeqNMF', 'W_hats_SeqNMF', 'H_hats_SeqNMF', "pvals_SeqNMF", "is_significants_SeqNMF", 'Ws', 'Hs', 'TrainingDatas')
end

%% Impact of noise
if n==2
    disp('Impact of noise level on results')
    Trials = 200;
    L = 50; % length of each trial
    K = 10;
    Nmotifs = 2*(1:K);
    Nneurons = 5*ones(K, 1); % the number of neurons in each motif
    Dt = 3.*ones(K,1); % gap between each member of the motif
    noise_levels = 0:0.01:0.1; % probability of added noise in each bin
    
    nSim = 50;
    seeds = randperm(1000, nSim);
    pvals_SeqNMF = cell(nSim,length(noise_levels));
    is_significants_SeqNMF = cell(nSim,length(noise_levels));
    loadings_SeqNMF = cell(nSim,length(noise_levels));
    W_hats_SeqNMF = cell(nSim,length(noise_levels));
    H_hats_SeqNMF = cell(nSim,length(noise_levels));
    pvals_FlexMF = pvals_SeqNMF;
    is_significants_FlexMF = is_significants_SeqNMF;
    loadings_FlexMF = loadings_SeqNMF;
    W_hats_FlexMF = W_hats_SeqNMF;
    H_hats_FlexMF = H_hats_SeqNMF;
    Ws = W_hats_SeqNMF;
    Hs = H_hats_SeqNMF;
    TrainingDatas = W_hats_SeqNMF;
    
    % Flattened (level j, sim i) -> single parfor over all nLevels*nSim tasks.
    nLevels = length(noise_levels);
    parfor idx=1:nSim*nLevels
        [i, j] = ind2sub([nSim, nLevels], idx);
        noise = noise_levels(j);
        tic
            [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, ...
            'noise', noise, 'seed', seeds(i));
            groups = zeros(Trials,1);
            for k=1:K
                groups(motif_ind{k}) = k;
            end
            Ws{idx} = W;
            
            % Dimension N*L*Trials
            cv = cvpartition(groups, "KFold",2);
            ind_train = find(training(cv,1));
            X_train = X(:,:,ind_train);
            ind_test = find(test(cv,1));
            X_test = X(:,:,ind_test);
            
            % Dimension N*T
            N = size(W,1);
            TrainingData = zeros(N,cv.TrainSize(1)*L);
            for t=1:cv.TrainSize(1)
                TrainingData(:,(t-1)*L+1:t*L) = squeeze(X_train(:,:,t));
            end
            TestData = zeros(N,cv.TestSize(1)*L);
            for t=1:cv.TestSize(1)
                TestData(:,(t-1)*L+1:t*L) = squeeze(X_test(:,:,t));
            end
            
            % Normalize training data
            frob_norm = norm(TrainingData(:));
            TrainingData = TrainingData/frob_norm*K;

            H_train = H(:,ind_train);
            T = cv.TrainSize(1)*L;
            H_train_full = zeros(K,T);
            H_train_full(:,1:L:T) = H_train;
            Hs{idx} = H_train_full;
        
            lambda = .003;
        
            % Run SeqNMF with multiplication rule
            [W_hat, H_hat, ~,~,loadings_SeqNMF{idx},power]= seqNMF(TrainingData,'K',K,'L',L,...
                    'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
            p = .05;
            [pvals_SeqNMF{idx},is_significants_SeqNMF{idx}] = test_significance_new(TestData, W_hat,[],p);
            W_hats_SeqNMF{idx} = W_hat;
            H_hats_SeqNMF{idx} = H_hat;
            display(['SeqNMF run ' num2str(i) '/' num2str(nSim)])
        
            % Run SeqNMF with Bregman Iteration
            alpha_W = 1e-5;
            alpha_H = 1e-2;
            [W_hat, H_hat, ~,~,loadings_FlexMF{idx},power]= FlexMF(TrainingData,'K',K,'L',L,...
                    'lambda', lambda, 'alpha_W', alpha_W, 'alpha_H', alpha_H, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0, 'verbal', 0); 
            p = .05;
            [pvals_FlexMF{idx},is_significants_FlexMF{idx}] = test_significance_new(TestData, W_hat,[],p);
            W_hats_FlexMF{idx} = W_hat;
            H_hats_FlexMF{idx} = H_hat;
            TrainingDatas{idx} = TrainingData;
            display(['FlexMF run ' num2str(i) '/' num2str(nSim)])
            toc
    end
    
    save('simulate_results_noise.mat', "H_hats_FlexMF", 'W_hats_FlexMF', 'is_significants_FlexMF', 'pvals_FlexMF', 'loadings_FlexMF',...
        'loadings_SeqNMF', 'W_hats_SeqNMF', 'H_hats_SeqNMF', "pvals_SeqNMF", "is_significants_SeqNMF", 'Ws', 'Hs', 'TrainingDatas')
end

%% Impact of participation rate
if n==3
    disp('Impact of participation rate on results')
    Trials = 200;
    L = 50; % length of each trial
    K = 10;
    Nmotifs = 2*(1:K);
    Nneurons = 5*ones(K, 1); % the number of neurons in each motif
    Dt = 3.*ones(K,1); % gap between each member of the motif
    
    participation_rates = 1:-0.1:0.1; % Participation probability
    
    nSim = 50;
    seeds = randperm(1000, nSim);
    pvals_SeqNMF = cell(nSim,length(participation_rates));
    is_significants_SeqNMF = cell(nSim,length(participation_rates));
    loadings_SeqNMF = cell(nSim,length(participation_rates));
    W_hats_SeqNMF = cell(nSim,length(participation_rates));
    H_hats_SeqNMF = cell(nSim,length(participation_rates));
    pvals_FlexMF = pvals_SeqNMF;
    is_significants_FlexMF = is_significants_SeqNMF;
    loadings_FlexMF = loadings_SeqNMF;
    W_hats_FlexMF = W_hats_SeqNMF;
    H_hats_FlexMF = H_hats_SeqNMF;
    Ws = W_hats_SeqNMF;
    Hs = H_hats_SeqNMF;
    TrainingDatas = W_hats_SeqNMF;
    
    % Flattened (level j, sim i) -> single parfor over all nLevels*nSim tasks.
    nLevels = length(participation_rates);
    parfor idx=1:nSim*nLevels
        [i, j] = ind2sub([nSim, nLevels], idx);
        participation = participation_rates(j).*ones(K,1);
        tic
            [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, ...
            'participation', participation, 'seed', seeds(i));
            groups = zeros(Trials,1);
            for k=1:K
                groups(motif_ind{k}) = k;
            end
            Ws{idx} = W;
            
            % Dimension N*L*Trials
            cv = cvpartition(groups, "KFold",2);
            ind_train = find(training(cv,1));
            X_train = X(:,:,ind_train);
            ind_test = find(test(cv,1));
            X_test = X(:,:,ind_test);
            
            % Dimension N*T
            N = size(W,1);
            TrainingData = zeros(N,cv.TrainSize(1)*L);
            for t=1:cv.TrainSize(1)
                TrainingData(:,(t-1)*L+1:t*L) = squeeze(X_train(:,:,t));
            end
            TestData = zeros(N,cv.TestSize(1)*L);
            for t=1:cv.TestSize(1)
                TestData(:,(t-1)*L+1:t*L) = squeeze(X_test(:,:,t));
            end
            
            % Normalize training data
            frob_norm = norm(TrainingData(:));
            TrainingData = TrainingData/frob_norm*K;

            H_train = H(:,ind_train);
            T = cv.TrainSize(1)*L;
            H_train_full = zeros(K,T);
            H_train_full(:,1:L:T) = H_train;
            Hs{idx} = H_train_full;
        
            lambda = .01;
        
            % Run SeqNMF with multiplication rule
            [W_hat, H_hat, ~,~,loadings_SeqNMF{idx},power]= seqNMF(TrainingData,'K',K,'L',L,...
                    'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
            p = .05;
            [pvals_SeqNMF{idx},is_significants_SeqNMF{idx}] = test_significance_new(TestData, W_hat,[],p);
            W_hats_SeqNMF{idx} = W_hat;
            H_hats_SeqNMF{idx} = H_hat;
            display(['SeqNMF run ' num2str(i) '/' num2str(nSim)])
        
            % Run SeqNMF with Bregman Iteration
            alpha_W = 1e-6;
            alpha_H = 1e-3;
            [W_hat, H_hat, ~,~,loadings_FlexMF{idx},power]= FlexMF(TrainingData,'K',K,'L',L,...
                    'lambda', lambda, 'alpha_W', alpha_W, 'alpha_H', alpha_H, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0, 'verbal', 0); 
            p = .05;
            [pvals_FlexMF{idx},is_significants_FlexMF{idx}] = test_significance_new(TestData, W_hat,[],p);
            W_hats_FlexMF{idx} = W_hat;
            H_hats_FlexMF{idx} = H_hat;
            TrainingDatas{idx} = TrainingData;
            display(['FlexMF run ' num2str(i) '/' num2str(nSim)])
            toc
    end
    
    save('simulate_results_participate.mat', "H_hats_FlexMF", 'W_hats_FlexMF', 'is_significants_FlexMF', 'pvals_FlexMF', 'loadings_FlexMF',...
        'loadings_SeqNMF', 'W_hats_SeqNMF', 'H_hats_SeqNMF', "pvals_SeqNMF", "is_significants_SeqNMF", 'Ws', 'Hs', 'TrainingDatas')
end

%% Impact of jittering
if n==4
    disp('Impact of jittering on results')
    Trials = 200;
    L = 100; % length of each trial
    K = 10;
    Nmotifs = 2*(1:K);
    Nneurons = 5*ones(K, 1); % the number of neurons in each motif
    Dt = 3.*ones(K,1); % gap between each member of the motif
    jitter_stds = 0:9; % Jitter time std
    
    nSim = 50;
    seeds = randperm(1000, nSim);
    pvals_SeqNMF = cell(nSim,length(jitter_stds));
    is_significants_SeqNMF = cell(nSim,length(jitter_stds));
    loadings_SeqNMF = cell(nSim,length(jitter_stds));
    W_hats_SeqNMF = cell(nSim,length(jitter_stds));
    H_hats_SeqNMF = cell(nSim,length(jitter_stds));
    pvals_FlexMF = pvals_SeqNMF;
    is_significants_FlexMF = is_significants_SeqNMF;
    loadings_FlexMF = loadings_SeqNMF;
    W_hats_FlexMF = W_hats_SeqNMF;
    H_hats_FlexMF = H_hats_SeqNMF;
    Ws = W_hats_SeqNMF;
    Hs = H_hats_SeqNMF;
    TrainingDatas = W_hats_SeqNMF;
    
    % Flattened (level j, sim i) loop -> single parfor over all nLevels*nSim
    % independent tasks, so the whole worker pool stays busy with no serial
    % per-level barriers. Linear index idx maps to (i,j) of the nSim-by-nLevels
    % cell arrays via column-major ind2sub, so results are identical.
    nLevels = length(jitter_stds);
    parfor idx=1:nSim*nLevels
        [i, j] = ind2sub([nSim, nLevels], idx);
        jitter = jitter_stds(j).*ones(K,1);
        tic
        [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, ...
        'jitter', jitter, 'seed', seeds(i));
        groups = zeros(Trials,1);
        for k=1:K
            groups(motif_ind{k}) = k;
        end
        Ws{idx} = W;

        % Dimension N*L*Trials
        cv = cvpartition(groups, "KFold",2);
        ind_train = find(training(cv,1));
        X_train = X(:,:,ind_train);
        ind_test = find(test(cv,1));
        X_test = X(:,:,ind_test);

        % Dimension N*T
        N = size(W,1);
        TrainingData = zeros(N,cv.TrainSize(1)*L);
        for t=1:cv.TrainSize(1)
            TrainingData(:,(t-1)*L+1:t*L) = squeeze(X_train(:,:,t));
        end
        TestData = zeros(N,cv.TestSize(1)*L);
        for t=1:cv.TestSize(1)
            TestData(:,(t-1)*L+1:t*L) = squeeze(X_test(:,:,t));
        end

        % Normalize training data
        frob_norm = norm(TrainingData(:));
        TrainingData = TrainingData/frob_norm*K;

        H_train = H(:,ind_train);
        T = cv.TrainSize(1)*L;
        H_train_full = zeros(K,T);
        H_train_full(:,1:L:T) = H_train;
        Hs{idx} = H_train_full;

        lambda = .01;

        % Run SeqNMF with multiplication rule
        [W_hat, H_hat, ~,~,loadings_SeqNMF{idx},power]= seqNMF(TrainingData,'K',K,'L',L,...
                'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
        p = .05;
        [pvals_SeqNMF{idx},is_significants_SeqNMF{idx}] = test_significance_new(TestData, W_hat,[],p);
        W_hats_SeqNMF{idx} = W_hat;
        H_hats_SeqNMF{idx} = H_hat;
        display(['SeqNMF run ' num2str(idx) '/' num2str(nSim*nLevels)])

        % Run SeqNMF with Bregman Iteration
        alpha_W = 1e-6;
        alpha_H = 1e-3;
        [W_hat, H_hat, ~,~,loadings_FlexMF{idx},power]= FlexMF(TrainingData,'K',K,'L',L,...
                'lambda', lambda, 'alpha_W', alpha_W, 'alpha_H', alpha_H, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0, 'verbal', 0); 
        p = .05;
        [pvals_FlexMF{idx},is_significants_FlexMF{idx}] = test_significance_new(TestData, W_hat,[],p);
        W_hats_FlexMF{idx} = W_hat;
        H_hats_FlexMF{idx} = H_hat;
        TrainingDatas{idx} = TrainingData;
        display(['FlexMF run ' num2str(idx) '/' num2str(nSim*nLevels)])
        toc
    end
    
    save('simulate_results_jitter.mat', "H_hats_FlexMF", 'W_hats_FlexMF', 'is_significants_FlexMF', 'pvals_FlexMF', 'loadings_FlexMF',...
        'loadings_SeqNMF', 'W_hats_SeqNMF', 'H_hats_SeqNMF', "pvals_SeqNMF", "is_significants_SeqNMF", 'Ws', 'Hs', 'TrainingDatas')
end

%% Impact of warping
if n==5
    disp('Impact of warping on results')
    Trials = 200;
    L = 100; % length of each trial
    K = 10;
    Nmotifs = 2*(1:K);
    Nneurons = 5*ones(K, 1); % the number of neurons in each motif
    Magnitudes = ones(K, 1); % the activation magnitudes of each motif
    Dt = 3.*ones(K,1); % gap between each member of the motif
    warp_levels = 0:9; % the maximum warping time
    
    nSim = 50;
    seeds = randperm(1000, nSim);
    pvals_SeqNMF = cell(nSim,length(warp_levels));
    is_significants_SeqNMF = cell(nSim,length(warp_levels));
    loadings_SeqNMF = cell(nSim,length(warp_levels));
    W_hats_SeqNMF = cell(nSim,length(warp_levels));
    H_hats_SeqNMF = cell(nSim,length(warp_levels));
    pvals_FlexMF = pvals_SeqNMF;
    is_significants_FlexMF = is_significants_SeqNMF;
    loadings_FlexMF = loadings_SeqNMF;
    W_hats_FlexMF = W_hats_SeqNMF;
    H_hats_FlexMF = H_hats_SeqNMF;
    Ws = W_hats_SeqNMF;
    Hs = H_hats_SeqNMF;
    TrainingDatas = W_hats_SeqNMF;
    
    % Flattened (level j, sim i) -> single parfor over all nLevels*nSim tasks.
    nLevels = length(warp_levels);
    parfor idx=1:nSim*nLevels
        [i, j] = ind2sub([nSim, nLevels], idx);
        warp = warp_levels(j);
        tic
            [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, ...
            'warp', warp, 'seed', seeds(i));
            groups = zeros(Trials,1);
            for k=1:K
                groups(motif_ind{k}) = k;
            end
            Ws{idx} = W;
            
            % Dimension N*L*Trials
            cv = cvpartition(groups, "KFold",2);
            ind_train = find(training(cv,1));
            X_train = X(:,:,ind_train);
            ind_test = find(test(cv,1));
            X_test = X(:,:,ind_test);
            
            % Dimension N*T
            N = size(W,1);
            TrainingData = zeros(N,cv.TrainSize(1)*L);
            for t=1:cv.TrainSize(1)
                TrainingData(:,(t-1)*L+1:t*L) = squeeze(X_train(:,:,t));
            end
            TestData = zeros(N,cv.TestSize(1)*L);
            for t=1:cv.TestSize(1)
                TestData(:,(t-1)*L+1:t*L) = squeeze(X_test(:,:,t));
            end
            
            % Normalize training data
            frob_norm = norm(TrainingData(:));
            TrainingData = TrainingData/frob_norm*K;

            H_train = H(:,ind_train);
            T = cv.TrainSize(1)*L;
            H_train_full = zeros(K,T);
            H_train_full(:,1:L:T) = H_train;
            Hs{idx} = H_train_full;
        
            lambda = .01;
        
            % Run SeqNMF with multiplication rule
            [W_hat, H_hat, ~,~,loadings_SeqNMF{idx},power]= seqNMF(TrainingData,'K',K,'L',L,...
                    'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
            p = .05;
            [pvals_SeqNMF{idx},is_significants_SeqNMF{idx}] = test_significance_new(TestData, W_hat,[],p);
            W_hats_SeqNMF{idx} = W_hat;
            H_hats_SeqNMF{idx} = H_hat;
            display(['SeqNMF run ' num2str(i) '/' num2str(nSim)])
        
            % Run SeqNMF with Bregman Iteration
            alpha_W = 1e-6;
            alpha_H = 1e-3;
            [W_hat, H_hat, ~,~,loadings_FlexMF{idx},power]= FlexMF(TrainingData,'K',K,'L',L,...
                    'lambda', lambda, 'alpha_W', alpha_W, 'alpha_H', alpha_H, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0, 'verbal', 0); 
            p = .05;
            [pvals_FlexMF{idx},is_significants_FlexMF{idx}] = test_significance_new(TestData, W_hat,[],p);
            W_hats_FlexMF{idx} = W_hat;
            H_hats_FlexMF{idx} = H_hat;
            TrainingDatas{idx} = TrainingData;
            display(['FlexMF run ' num2str(i) '/' num2str(nSim)])
            toc
    end
    
    save('simulate_results_warp.mat', "H_hats_FlexMF", 'W_hats_FlexMF', 'is_significants_FlexMF', 'pvals_FlexMF', 'loadings_FlexMF',...
        'loadings_SeqNMF', 'W_hats_SeqNMF', 'H_hats_SeqNMF', "pvals_SeqNMF", "is_significants_SeqNMF", 'Ws', 'Hs', 'TrainingDatas')
end

%% Impact of shared neuron
if n==6
    disp('Impact of shared neurons on results')
    Trials = 200;
    L = 50; % length of each trial
    K = 10;
    Nmotifs = 2*(1:K);
    Nneurons = 5*ones(K, 1); % the number of neurons in each motif
    Dt = 3.*ones(K,1); % gap between each member of the motif
    
    overlaps_n = 0:0.1:0.9; % Neuron overlap rate
    
    nSim = 50;
    seeds = randperm(1000, nSim);
    pvals_SeqNMF = cell(nSim,length(overlaps_n));
    is_significants_SeqNMF = cell(nSim,length(overlaps_n));
    loadings_SeqNMF = cell(nSim,length(overlaps_n));
    W_hats_SeqNMF = cell(nSim,length(overlaps_n));
    H_hats_SeqNMF = cell(nSim,length(overlaps_n));
    pvals_FlexMF = pvals_SeqNMF;
    is_significants_FlexMF = is_significants_SeqNMF;
    loadings_FlexMF = loadings_SeqNMF;
    W_hats_FlexMF = W_hats_SeqNMF;
    H_hats_FlexMF = H_hats_SeqNMF;
    Ws = W_hats_SeqNMF;
    Hs = H_hats_SeqNMF;
    TrainingDatas = W_hats_SeqNMF;
    
    % Flattened (level j, sim i) -> single parfor over all nLevels*nSim tasks.
    nLevels = length(overlaps_n);
    parfor idx=1:nSim*nLevels
        [i, j] = ind2sub([nSim, nLevels], idx);
        overlap_n = overlaps_n(j);
        tic
            [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, ...
                'overlap_n', overlap_n, 'seed', seeds(i));
            groups = zeros(Trials,1);
            for k=1:K
                groups(motif_ind{k}) = k;
            end
            Ws{idx} = W;
            
            % Dimension N*L*Trials
            cv = cvpartition(groups, "KFold",2);
            ind_train = find(training(cv,1));
            X_train = X(:,:,ind_train);
            ind_test = find(test(cv,1));
            X_test = X(:,:,ind_test);
            
            % Dimension N*T
            N = size(W,1);
            TrainingData = zeros(N,cv.TrainSize(1)*L);
            for t=1:cv.TrainSize(1)
                TrainingData(:,(t-1)*L+1:t*L) = squeeze(X_train(:,:,t));
            end
            TestData = zeros(N,cv.TestSize(1)*L);
            for t=1:cv.TestSize(1)
                TestData(:,(t-1)*L+1:t*L) = squeeze(X_test(:,:,t));
            end
            
            % Normalize training data
            frob_norm = norm(TrainingData(:));
            TrainingData = TrainingData/frob_norm*K;

            H_train = H(:,ind_train);
            T = cv.TrainSize(1)*L;
            H_train_full = zeros(K,T);
            H_train_full(:,1:L:T) = H_train;
            Hs{idx} = H_train_full;
        
            lambda = .01;
        
            % Run SeqNMF with multiplication rule
            [W_hat, H_hat, ~,~,loadings_SeqNMF{idx},power]= seqNMF(TrainingData,'K',K,'L',L,...
                    'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
            p = .05;
            [pvals_SeqNMF{idx},is_significants_SeqNMF{idx}] = test_significance_new(TestData, W_hat,[],p);
            W_hats_SeqNMF{idx} = W_hat;
            H_hats_SeqNMF{idx} = H_hat;
            display(['SeqNMF run ' num2str(i) '/' num2str(nSim)])
        
            % Run SeqNMF with Bregman Iteration
            alpha_W = 1e-6;
            alpha_H = 1e-3;
            [W_hat, H_hat, ~,~,loadings_FlexMF{idx},power]= FlexMF(TrainingData,'K',K,'L',L,...
                    'lambda', lambda, 'alpha_W', alpha_W, 'alpha_H', alpha_H, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0, 'verbal', 0); 
            p = .05;
            [pvals_FlexMF{idx},is_significants_FlexMF{idx}] = test_significance_new(TestData, W_hat,[],p);
            W_hats_FlexMF{idx} = W_hat;
            H_hats_FlexMF{idx} = H_hat;
            TrainingDatas{idx} = TrainingData;
            display(['FlexMF run ' num2str(i) '/' num2str(nSim)])
            toc
    end
    
    save('simulate_results_overlap.mat', "H_hats_FlexMF", 'W_hats_FlexMF', 'is_significants_FlexMF', 'pvals_FlexMF', 'loadings_FlexMF',...
        'loadings_SeqNMF', 'W_hats_SeqNMF', 'H_hats_SeqNMF', "pvals_SeqNMF", "is_significants_SeqNMF", 'Ws', 'Hs', 'TrainingDatas')
end

% Shut down the parallel pool
delete(gcp('nocreate'));
