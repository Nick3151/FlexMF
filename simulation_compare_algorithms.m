% Simulate data with trials, compare the performance of SeqNMF with
% multiplication update rule and Split Bregman Iteration

clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(fullfile(root, 'Utils'))
addpath(fullfile(root, 'MATLAB-tools'))
addpath(genpath(fullfile(root, 'FlexMF')));
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
    
    nSim = 100;
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
    
    parfor i=1:nSim
        tic
        disp('Calcium Transients:')
        len_burst = 10; % Continuous firing time
        dynamic = 1; % Consider calcium dynamics or not
        [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, ...
            'len_burst', len_burst, 'dynamic', dynamic, 'seed', seeds(i));
        groups = zeros(Trials,1);
        for k=1:K
            groups(motif_ind{k}) = k;
        end
        Ws{i,1} = W;
        
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
        Hs{i,1} = H_train_full;
    
        lambda = .005;
    
        % Run SeqNMF with multiplication rule
        [W_hat, H_hat, ~,~,loadings_SeqNMF{i,1},power]= seqNMF(TrainingData,'K',K,'L',L,...
                'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
        p = .05;
        [pvals_SeqNMF{i,1},is_significants_SeqNMF{i,1}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
        W_hats_SeqNMF{i,1} = W_hat;
        H_hats_SeqNMF{i,1} = H_hat;
        display(['SeqNMF run ' num2str(i) '/' num2str(nSim)])
    
        % Run SeqNMF with Bregman Iteration
        alpha_W = 1e-6;
        alpha_H = 1e-2;
        [W_hat, H_hat, ~,~,loadings_FlexMF{i,1},power]= FlexMF(TrainingData,'K',K,'L',L,...
                'lambda', lambda, 'alpha_W', alpha_W, 'alpha_H', alpha_H, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0, 'verbal', 0); 
        p = .05;
        [pvals_FlexMF{i,1},is_significants_FlexMF{i,1}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
        W_hats_FlexMF{i,1} = W_hat;
        H_hats_FlexMF{i,1} = H_hat;
        TrainingDatas{i,1} = TrainingData;
        display(['FlexMF run ' num2str(i) '/' num2str(nSim)])
        toc
    end
    
    parfor i=1:nSim
        tic
        disp('Burst:')
        len_burst = 10; % Continuous firing time
        dynamic = 0; % Consider calcium dynamics or not
        [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, ...
            'len_burst', len_burst, 'dynamic', dynamic, 'seed', seeds(i));
        groups = zeros(Trials,1);
        for k=1:K
            groups(motif_ind{k}) = k;
        end
        Ws{i,2} = W
        
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
        Hs{i,2} = H_train_full;
    
        lambda = .005;
    
        % Run SeqNMF with multiplication rule
        [W_hat, H_hat, ~,~,loadings_SeqNMF{i,2},power]= seqNMF(TrainingData,'K',K,'L',L,...
                'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
        p = .05;
        [pvals_SeqNMF{i,2},is_significants_SeqNMF{i,2}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
        W_hats_SeqNMF{i,2} = W_hat;
        H_hats_SeqNMF{i,2} = H_hat;
        display(['SeqNMF run ' num2str(i) '/' num2str(nSim)])
    
        % Run SeqNMF with Bregman Iteration
        alpha_W = 1e-6;
        alpha_H = 1e-2;
        [W_hat, H_hat, ~,~,loadings_FlexMF{i,2},power]= FlexMF(TrainingData,'K',K,'L',L,...
                'lambda', lambda, 'alpha_W', alpha_W, 'alpha_H', alpha_H, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0, 'verbal', 0); 
        p = .05;
        [pvals_FlexMF{i,2},is_significants_FlexMF{i,2}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
        W_hats_FlexMF{i,2} = W_hat;
        H_hats_FlexMF{i,2} = H_hat;
        TrainingDatas{i,2} = TrainingData;
        display(['FlexMF run ' num2str(i) '/' num2str(nSim)])
        toc
    end
    
    parfor i=1:nSim
        tic
        disp('Single spikes:')
        len_burst = 1; % Continuous firing time
        dynamic = 0; % Consider calcium dynamics or not
        [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, ...
            'len_burst', len_burst, 'dynamic', dynamic, 'seed', seeds(i));
        groups = zeros(Trials,1);
        for k=1:K
            groups(motif_ind{k}) = k;
        end
        Ws{i,3} = W;
        
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
        Hs{i,3} = H_train_full;
    
        lambda = .01;
    
        % Run SeqNMF with multiplication rule
        [W_hat, H_hat, ~,~,loadings_SeqNMF{i,3},power]= seqNMF(TrainingData,'K',K,'L',L,...
                'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
        p = .05;
        [pvals_SeqNMF{i,3},is_significants_SeqNMF{i,3}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
        W_hats_SeqNMF{i,3} = W_hat;
        H_hats_SeqNMF{i,3} = H_hat;
        display(['SeqNMF run ' num2str(i) '/' num2str(nSim)])
    
        % Run SeqNMF with Bregman Iteration
        alpha_W = 1e-6;
        alpha_H = 1e-3;
        [W_hat, H_hat, ~,~,loadings_FlexMF{i,3},power]= FlexMF(TrainingData,'K',K,'L',L,...
                'lambda', lambda, 'alpha_W', alpha_W, 'alpha_H', alpha_H, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0, 'verbal', 0); 
        p = .05;
        [pvals_FlexMF{i,3},is_significants_FlexMF{i,3}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
        W_hats_FlexMF{i,3} = W_hat;
        H_hats_FlexMF{i,3} = H_hat;
        TrainingDatas{i,3} = TrainingData;
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
    noise_levels = [0.001, 0.01, 0.1]; % probability of added noise in each bin
    
    nSim = 100;
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
    
    for j=1:length(noise_levels)
        noise = noise_levels(j);
        disp(['Noise level=', num2str(noise_levels(j))])
        parfor i=1:nSim
            tic
            [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, ...
            'noise', noise, 'seed', seeds(i));
            groups = zeros(Trials,1);
            for k=1:K
                groups(motif_ind{k}) = k;
            end
            Ws{i,j} = W;
            
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
            Hs{i,j} = H_train_full;
        
            lambda = .003;
        
            % Run SeqNMF with multiplication rule
            [W_hat, H_hat, ~,~,loadings_SeqNMF{i,j},power]= seqNMF(TrainingData,'K',K,'L',L,...
                    'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
            p = .05;
            [pvals_SeqNMF{i,j},is_significants_SeqNMF{i,j}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
            W_hats_SeqNMF{i,j} = W_hat;
            H_hats_SeqNMF{i,j} = H_hat;
            display(['SeqNMF run ' num2str(i) '/' num2str(nSim)])
        
            % Run SeqNMF with Bregman Iteration
            alpha_W = 1e-5;
            alpha_H = 1e-2;
            [W_hat, H_hat, ~,~,loadings_FlexMF{i,j},power]= FlexMF(TrainingData,'K',K,'L',L,...
                    'lambda', lambda, 'alpha_W', alpha_W, 'alpha_H', alpha_H, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0, 'verbal', 0); 
            p = .05;
            [pvals_FlexMF{i,j},is_significants_FlexMF{i,j}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
            W_hats_FlexMF{i,j} = W_hat;
            H_hats_FlexMF{i,j} = H_hat;
            TrainingDatas{i,j} = TrainingData;
            display(['FlexMF run ' num2str(i) '/' num2str(nSim)])
            toc
        end
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
    
    participation_rates = 0.5:0.1:0.9; % Participation probability
    
    nSim = 100;
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
    
    for j=1:length(participation_rates)
        participation = participation_rates(j).*ones(K,1);
        disp(['Participation rate=', num2str(participation_rates(j))])
        parfor i=1:nSim
            tic
            [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, ...
            'participation', participation, 'seed', seeds(i));
            groups = zeros(Trials,1);
            for k=1:K
                groups(motif_ind{k}) = k;
            end
            Ws{i,j} = W;
            
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
            Hs{i,j} = H_train_full;
        
            lambda = .01;
        
            % Run SeqNMF with multiplication rule
            [W_hat, H_hat, ~,~,loadings_SeqNMF{i,j},power]= seqNMF(TrainingData,'K',K,'L',L,...
                    'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
            p = .05;
            [pvals_SeqNMF{i,j},is_significants_SeqNMF{i,j}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
            W_hats_SeqNMF{i,j} = W_hat;
            H_hats_SeqNMF{i,j} = H_hat;
            display(['SeqNMF run ' num2str(i) '/' num2str(nSim)])
        
            % Run SeqNMF with Bregman Iteration
            alpha_W = 1e-6;
            alpha_H = 1e-3;
            [W_hat, H_hat, ~,~,loadings_FlexMF{i,j},power]= FlexMF(TrainingData,'K',K,'L',L,...
                    'lambda', lambda, 'alpha_W', alpha_W, 'alpha_H', alpha_H, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0, 'verbal', 0); 
            p = .05;
            [pvals_FlexMF{i,j},is_significants_FlexMF{i,j}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
            W_hats_FlexMF{i,j} = W_hat;
            H_hats_FlexMF{i,j} = H_hat;
            TrainingDatas{i,j} = TrainingData;
            display(['FlexMF run ' num2str(i) '/' num2str(nSim)])
            toc
        end
    end
    
    save('simulate_results_participate.mat', "H_hats_FlexMF", 'W_hats_FlexMF', 'is_significants_FlexMF', 'pvals_FlexMF', 'loadings_FlexMF',...
        'loadings_SeqNMF', 'W_hats_SeqNMF', 'H_hats_SeqNMF', "pvals_SeqNMF", "is_significants_SeqNMF", 'Ws', 'Hs', 'TrainingDatas')
end

%% Impact of jittering
if n==4
    disp('Impact of jittering on results')
    Trials = 200;
    L = 50; % length of each trial
    K = 10;
    Nmotifs = 2*(1:K);
    Nneurons = 5*ones(K, 1); % the number of neurons in each motif
    Dt = 3.*ones(K,1); % gap between each member of the motif
    jitter_stds = 1:3; % Jitter time std
    
    nSim = 100;
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
    
    for j=1:length(jitter_stds)
        jitter = jitter_stds(j).*ones(K,1);
        disp(['Jittering std=', num2str(jitter_stds(j))])
        parfor i=1:nSim
            tic
            [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, ...
            'jitter', jitter, 'seed', seeds(i));
            groups = zeros(Trials,1);
            for k=1:K
                groups(motif_ind{k}) = k;
            end
            Ws{i,j} = W;
            
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
            Hs{i,j} = H_train_full;
        
            lambda = .01;
        
            % Run SeqNMF with multiplication rule
            [W_hat, H_hat, ~,~,loadings_SeqNMF{i,j},power]= seqNMF(TrainingData,'K',K,'L',L,...
                    'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
            p = .05;
            [pvals_SeqNMF{i,j},is_significants_SeqNMF{i,j}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
            W_hats_SeqNMF{i,j} = W_hat;
            H_hats_SeqNMF{i,j} = H_hat;
            display(['SeqNMF run ' num2str(i) '/' num2str(nSim)])
        
            % Run SeqNMF with Bregman Iteration
            alpha_W = 1e-6;
            alpha_H = 1e-3;
            [W_hat, H_hat, ~,~,loadings_FlexMF{i,j},power]= FlexMF(TrainingData,'K',K,'L',L,...
                    'lambda', lambda, 'alpha_W', alpha_W, 'alpha_H', alpha_H, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0, 'verbal', 0); 
            p = .05;
            [pvals_FlexMF{i,j},is_significants_FlexMF{i,j}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
            W_hats_FlexMF{i,j} = W_hat;
            H_hats_FlexMF{i,j} = H_hat;
            TrainingDatas{i,j} = TrainingData;
            display(['FlexMF run ' num2str(i) '/' num2str(nSim)])
            toc
        end
    end
    
    save('simulate_results_jitter.mat', "H_hats_FlexMF", 'W_hats_FlexMF', 'is_significants_FlexMF', 'pvals_FlexMF', 'loadings_FlexMF',...
        'loadings_SeqNMF', 'W_hats_SeqNMF', 'H_hats_SeqNMF', "pvals_SeqNMF", "is_significants_SeqNMF", 'Ws', 'Hs', 'TrainingDatas')
end

%% Impact of warping
if n==5
    disp('Impact of warping on results')
    Trials = 200;
    L = 50; % length of each trial
    K = 10;
    Nmotifs = 2*(1:K);
    Nneurons = 5*ones(K, 1); % the number of neurons in each motif
    Magnitudes = ones(K, 1); % the activation magnitudes of each motif
    Dt = 3.*ones(K,1); % gap between each member of the motif
    warp_levels = [1,2,3]; % the maximum warping time
    
    nSim = 100;
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
    
    for j=1:length(warp_levels)
        warp = warp_levels(j);
        disp(['Warping level=', num2str(warp_levels(j))])
        parfor i=1:nSim
            tic
            [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, ...
            'warp', warp, 'seed', seeds(i));
            groups = zeros(Trials,1);
            for k=1:K
                groups(motif_ind{k}) = k;
            end
            Ws{i,j} = W;
            
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
            Hs{i,j} = H_train_full;
        
            lambda = .01;
        
            % Run SeqNMF with multiplication rule
            [W_hat, H_hat, ~,~,loadings_SeqNMF{i,j},power]= seqNMF(TrainingData,'K',K,'L',L,...
                    'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
            p = .05;
            [pvals_SeqNMF{i,j},is_significants_SeqNMF{i,j}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
            W_hats_SeqNMF{i,j} = W_hat;
            H_hats_SeqNMF{i,j} = H_hat;
            display(['SeqNMF run ' num2str(i) '/' num2str(nSim)])
        
            % Run SeqNMF with Bregman Iteration
            alpha_W = 1e-6;
            alpha_H = 1e-3;
            [W_hat, H_hat, ~,~,loadings_FlexMF{i,j},power]= FlexMF(TrainingData,'K',K,'L',L,...
                    'lambda', lambda, 'alpha_W', alpha_W, 'alpha_H', alpha_H, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0, 'verbal', 0); 
            p = .05;
            [pvals_FlexMF{i,j},is_significants_FlexMF{i,j}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
            W_hats_FlexMF{i,j} = W_hat;
            H_hats_FlexMF{i,j} = H_hat;
            TrainingDatas{i,j} = TrainingData;
            display(['FlexMF run ' num2str(i) '/' num2str(nSim)])
            toc
        end
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
    
    overlaps_n = 0:0.2:0.8; % Participation probability
    
    nSim = 100;
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
    
    for j=1:length(overlaps_n)
        overlap_n = overlaps_n(j);
        disp(['neurons overlap rate=', num2str(overlaps_n(j))])
        parfor i=1:nSim
            tic
            [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, ...
                'overlap_n', overlap_n, 'seed', seeds(i));
            groups = zeros(Trials,1);
            for k=1:K
                groups(motif_ind{k}) = k;
            end
            Ws{i,j} = W;
            
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
            Hs{i,j} = H_train_full;
        
            lambda = .01;
        
            % Run SeqNMF with multiplication rule
            [W_hat, H_hat, ~,~,loadings_SeqNMF{i,j},power]= seqNMF(TrainingData,'K',K,'L',L,...
                    'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
            p = .05;
            [pvals_SeqNMF{i,j},is_significants_SeqNMF{i,j}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
            W_hats_SeqNMF{i,j} = W_hat;
            H_hats_SeqNMF{i,j} = H_hat;
            display(['SeqNMF run ' num2str(i) '/' num2str(nSim)])
        
            % Run SeqNMF with Bregman Iteration
            alpha_W = 1e-6;
            alpha_H = 1e-3;
            [W_hat, H_hat, ~,~,loadings_FlexMF{i,j},power]= FlexMF(TrainingData,'K',K,'L',L,...
                    'lambda', lambda, 'alpha_W', alpha_W, 'alpha_H', alpha_H, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0, 'verbal', 0); 
            p = .05;
            [pvals_FlexMF{i,j},is_significants_FlexMF{i,j}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
            W_hats_FlexMF{i,j} = W_hat;
            H_hats_FlexMF{i,j} = H_hat;
            TrainingDatas{i,j} = TrainingData;
            display(['FlexMF run ' num2str(i) '/' num2str(nSim)])
            toc
        end
    end
    
    save('simulate_results_overlap.mat', "H_hats_FlexMF", 'W_hats_FlexMF', 'is_significants_FlexMF', 'pvals_FlexMF', 'loadings_FlexMF',...
        'loadings_SeqNMF', 'W_hats_SeqNMF', 'H_hats_SeqNMF', "pvals_SeqNMF", "is_significants_SeqNMF", 'Ws', 'Hs', 'TrainingDatas')
end