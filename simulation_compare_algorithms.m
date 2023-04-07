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

%% Simulate data for one condition
Trials = 200;
% Trials = 10;
L = 50; % length of each trial

Nmotifs = 2*(1:10);
K = 10;
Nneurons = 5*ones(K, 1); % the number of neurons in each motif
Magnitudes = ones(K, 1); % the activation magnitudes of each motif
Dt = 3.*ones(K,1); % gap between each member of the motif
noise = 0.01; % probability of added noise in each bin
jitter = 0*ones(K,1); % Jitter time std
participation = .7.*ones(K,1); % Participation probability = 100%
warp = 0; % the maximum warping time
len_burst = 20; % Burst time
dynamic = 1; % Consider calcium dynamics or not
overlap = 0;
neg = 0; % Proportion of negative indices in W

nSim = 100;
seeds = randperm(1000, nSim);
pvals_SeqNMF = cell(nSim,1);
is_significant_SeqNMF = cell(nSim,1);
loadings_SeqNMF = cell(nSim,1);
W_hats_SeqNMF = cell(nSim,1);
H_hats_SeqNMF = cell(nSim,1);
pvals_FlexMF = pvals_SeqNMF;
is_significant_FlexMF = is_significant_SeqNMF;
loadings_FlexMF = loadings_SeqNMF;
W_hats_FlexMF = W_hats_SeqNMF;
H_hats_FlexMF = H_hats_SeqNMF;
Ws = W_hats_SeqNMF;
Hs = H_hats_SeqNMF;

parpool(48)
parfor i=1:nSim
    tic
    [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Magnitudes, Dt, noise, jitter, participation, warp, len_burst, dynamic, overlap, neg, seeds(i));
    groups = zeros(Trials,1);
    for k=1:K
        groups(motif_ind{k}) = k;
    end
    Ws{i} = W;
    Hs{i} = H;
    
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
    nuc_norm = norm(svd(TrainingData),1);
    TrainingData = TrainingData/nuc_norm*N;

    lambda = .005;
    alpha = 5e-5;

    % Run SeqNMF with multiplication rule
    [W_hat, H_hat, ~,loadings_SeqNMF{i},power]= seqNMF(TrainingData,'K',K,'L',L,...
            'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
    p = .05;
    [pvals_SeqNMF{i},is_significant_SeqNMF{i}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
    W_hats_SeqNMF{i} = W_hat;
    H_hats_SeqNMF{i} = H_hat;
    display(['SeqNMF run ' num2str(i) '/' num2str(nSim)])

    % Run SeqNMF with Bregman Iteration
    [W_hat, H_hat, ~,loadings_FlexMF{i},power]= FlexMF(TrainingData,'K',K,'L',L,...
            'lambda', lambda, 'alpha', alpha, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0, 'verbal', 0); 
    p = .05;
    [pvals_FlexMF{i},is_significant_FlexMF{i}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
    W_hats_FlexMF{i} = W_hat;
    H_hats_FlexMF{i} = H_hat;
    display(['FlexMF run ' num2str(i) '/' num2str(nSim)])
    toc
end

save('simulate_results.mat', "H_hats_FlexMF", 'W_hats_FlexMF', 'is_significant_FlexMF', 'pvals_FlexMF', 'loadings_FlexMF',...
    'loadings_SeqNMF', 'W_hats_SeqNMF', 'H_hats_SeqNMF', "pvals_SeqNMF", "is_significant_SeqNMF", 'Ws', 'Hs')

%% Impact of motif shape
disp('Impact of motif shape on results')
Trials = 200;
L = 50; % length of each trial
Nmotifs = 2*(1:10);
K = 10;
Nneurons = 5*ones(K, 1); % the number of neurons in each motif
Magnitudes = ones(K, 1); % the activation magnitudes of each motif
Dt = 3.*ones(K,1); % gap between each member of the motif
noise = 0.01; % probability of added noise in each bin
jitter = 0*ones(K,1); % Jitter time std
participation = .7.*ones(K,1); % Participation probability = 100%
warp = 0; % the maximum warping time
len_burst = 20; % Continuous firing time
dynamic = 1; % Consider calcium dynamics or not
overlap = 0;
neg = 0; % Proportion of negative indices in W

nSim = 100;
seeds = randperm(1000, nSim);
pvals_SeqNMF = cell(nSim,3);
is_significant_SeqNMF = cell(nSim,3);
loadings_SeqNMF = cell(nSim,3);
W_hats_SeqNMF = cell(nSim,3);
H_hats_SeqNMF = cell(nSim,3);
pvals_FlexMF = pvals_SeqNMF;
is_significant_FlexMF = is_significant_SeqNMF;
loadings_FlexMF = loadings_SeqNMF;
W_hats_FlexMF = W_hats_SeqNMF;
H_hats_FlexMF = H_hats_SeqNMF;
Ws = W_hats_SeqNMF;
Hs = H_hats_SeqNMF;

parfor i=1:nSim
    tic
    disp('Calcium Transients:')
    len_burst = 20; % Continuous firing time
    dynamic = 1; % Consider calcium dynamics or not
    [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Magnitudes, Dt, noise, jitter, participation, warp, len_burst, dynamic, overlap, neg, seeds(i));
    groups = zeros(Trials,1);
    for k=1:K
        groups(motif_ind{k}) = k;
    end
    Ws{i} = W;
    Hs{i} = H;
    
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
    nuc_norm = norm(svd(TrainingData),1);
    TrainingData = TrainingData/nuc_norm*N;

    lambda = .005;
    alpha = 5e-5;

    % Run SeqNMF with multiplication rule
    [W_hat, H_hat, ~,loadings_SeqNMF{i,1},power]= seqNMF(TrainingData,'K',K,'L',L,...
            'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
    p = .05;
    [pvals_SeqNMF{i,1},is_significant_SeqNMF{i,1}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
    W_hats_SeqNMF{i,1} = W_hat;
    H_hats_SeqNMF{i,1} = H_hat;
    display(['SeqNMF run ' num2str(i) '/' num2str(nSim)])

    % Run SeqNMF with Bregman Iteration
    [W_hat, H_hat, ~,loadings_FlexMF{i,1},power]= FlexMF(TrainingData,'K',K,'L',L,...
            'lambda', lambda, 'alpha', alpha, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0, 'verbal', 0); 
    p = .05;
    [pvals_FlexMF{i,1},is_significant_FlexMF{i,1}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
    W_hats_FlexMF{i,1} = W_hat;
    H_hats_FlexMF{i,1} = H_hat;
    display(['FlexMF run ' num2str(i) '/' num2str(nSim)])
    toc
end

parfor i=1:nSim
    tic
    disp('Burst:')
    len_burst = 20; % Continuous firing time
    dynamic = 0; % Consider calcium dynamics or not
    [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Magnitudes, Dt, noise, jitter, participation, warp, len_burst, dynamic, overlap, neg, seeds(i));
    groups = zeros(Trials,1);
    for k=1:K
        groups(motif_ind{k}) = k;
    end
    Ws{i} = W;
    Hs{i} = H;
    
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
    nuc_norm = norm(svd(TrainingData),1);
    TrainingData = TrainingData/nuc_norm*N;

    lambda = .005;
    alpha = 5e-5;

    % Run SeqNMF with multiplication rule
    [W_hat, H_hat, ~,loadings_SeqNMF{i,2},power]= seqNMF(TrainingData,'K',K,'L',L,...
            'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
    p = .05;
    [pvals_SeqNMF{i,2},is_significant_SeqNMF{i,2}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
    W_hats_SeqNMF{i,2} = W_hat;
    H_hats_SeqNMF{i,2} = H_hat;
    display(['SeqNMF run ' num2str(i) '/' num2str(nSim)])

    % Run SeqNMF with Bregman Iteration
    [W_hat, H_hat, ~,loadings_FlexMF{i,2},power]= FlexMF(TrainingData,'K',K,'L',L,...
            'lambda', lambda, 'alpha', alpha, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0, 'verbal', 0); 
    p = .05;
    [pvals_FlexMF{i,2},is_significant_FlexMF{i,2}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
    W_hats_FlexMF{i,2} = W_hat;
    H_hats_FlexMF{i,2} = H_hat;
    display(['FlexMF run ' num2str(i) '/' num2str(nSim)])
    toc
end

parfor i=1:nSim
    tic
    disp('Single spikes:')
    len_burst = 1; % Continuous firing time
    dynamic = 0; % Consider calcium dynamics or not
    [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Magnitudes, Dt, noise, jitter, participation, warp, len_burst, dynamic, overlap, neg, seeds(i));
    groups = zeros(Trials,1);
    for k=1:K
        groups(motif_ind{k}) = k;
    end
    Ws{i} = W;
    Hs{i} = H;
    
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
    nuc_norm = norm(svd(TrainingData),1);
    TrainingData = TrainingData/nuc_norm*N;

    lambda = .005;
    alpha = 5e-5;

    % Run SeqNMF with multiplication rule
    [W_hat, H_hat, ~,loadings_SeqNMF{i,3},power]= seqNMF(TrainingData,'K',K,'L',L,...
            'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
    p = .05;
    [pvals_SeqNMF{i,3},is_significant_SeqNMF{i,3}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
    W_hats_SeqNMF{i,3} = W_hat;
    H_hats_SeqNMF{i,3} = H_hat;
    display(['SeqNMF run ' num2str(i) '/' num2str(nSim)])

    % Run SeqNMF with Bregman Iteration
    [W_hat, H_hat, ~,loadings_FlexMF{i,3},power]= FlexMF(TrainingData,'K',K,'L',L,...
            'lambda', lambda, 'alpha', alpha, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0, 'verbal', 0); 
    p = .05;
    [pvals_FlexMF{i,3},is_significant_FlexMF{i,3}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
    W_hats_FlexMF{i,3} = W_hat;
    H_hats_FlexMF{i,3} = H_hat;
    display(['FlexMF run ' num2str(i) '/' num2str(nSim)])
    toc
end

save('simulate_results_shape.mat', "H_hats_FlexMF", 'W_hats_FlexMF', 'is_significant_FlexMF', 'pvals_FlexMF', 'loadings_FlexMF',...
    'loadings_SeqNMF', 'W_hats_SeqNMF', 'H_hats_SeqNMF', "pvals_SeqNMF", "is_significant_SeqNMF", 'Ws', 'Hs')

%% Impact of participation rate
disp('Impact of participation rate on results')
Trials = 200;
L = 50; % length of each trial
Nmotifs = 2*(1:10);
K = 10;
Nneurons = 5*ones(K, 1); % the number of neurons in each motif
Magnitudes = ones(K, 1); % the activation magnitudes of each motif
Dt = 3.*ones(K,1); % gap between each member of the motif
noise = 0.001; % probability of added noise in each bin
jitter = 0*ones(K,1); % Jitter time std
participation_rates = 0.5:0.1:0.9; % Participation probability
warp = 0; % the maximum warping time
len_burst = 10; % Continuous firing time
dynamic = 1; % Consider calcium dynamics or not
overlap = 0;
neg = 0; % Proportion of negative indices in W

nSim = 100;
seeds = randperm(1000, nSim);
pvals_SeqNMF = cell(nSim,length(participation_rates));
is_significant_SeqNMF = cell(nSim,length(participation_rates));
loadings_SeqNMF = cell(nSim,length(participation_rates));
W_hats_SeqNMF = cell(nSim,length(participation_rates));
H_hats_SeqNMF = cell(nSim,length(participation_rates));
pvals_FlexMF = pvals_SeqNMF;
is_significant_FlexMF = is_significant_SeqNMF;
loadings_FlexMF = loadings_SeqNMF;
W_hats_FlexMF = W_hats_SeqNMF;
H_hats_FlexMF = H_hats_SeqNMF;
Ws = W_hats_SeqNMF;
Hs = H_hats_SeqNMF;

for j=1:length(participation_rates)
    participation = participation_rates(j).*ones(K,1);
    disp(['Participation rate=', num2str(participation_rates(j))])
    parfor i=1:nSim
        tic
        len_burst = 20; % Continuous firing time
        dynamic = 1; % Consider calcium dynamics or not
        [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Magnitudes, Dt, noise, jitter, participation, warp, len_burst, dynamic, overlap, neg, seeds(i));
        groups = zeros(Trials,1);
        for k=1:K
            groups(motif_ind{k}) = k;
        end
        Ws{i} = W;
        Hs{i} = H;
        
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
        nuc_norm = norm(svd(TrainingData),1);
        TrainingData = TrainingData/nuc_norm*N;
    
        lambda = .005;
        alpha = 5e-5;
    
        % Run SeqNMF with multiplication rule
        [W_hat, H_hat, ~,loadings_SeqNMF{i,j},power]= seqNMF(TrainingData,'K',K,'L',L,...
                'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
        p = .05;
        [pvals_SeqNMF{i,j},is_significant_SeqNMF{i,j}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
        W_hats_SeqNMF{i,j} = W_hat;
        H_hats_SeqNMF{i,j} = H_hat;
        display(['SeqNMF run ' num2str(i) '/' num2str(nSim)])
    
        % Run SeqNMF with Bregman Iteration
        [W_hat, H_hat, ~,loadings_FlexMF{i,j},power]= FlexMF(TrainingData,'K',K,'L',L,...
                'lambda', lambda, 'alpha', alpha, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0, 'verbal', 0); 
        p = .05;
        [pvals_FlexMF{i,j},is_significant_FlexMF{i,j}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
        W_hats_FlexMF{i,j} = W_hat;
        H_hats_FlexMF{i,j} = H_hat;
        display(['FlexMF run ' num2str(i) '/' num2str(nSim)])
        toc
    end
end

save('simulate_results_participate.mat', "H_hats_FlexMF", 'W_hats_FlexMF', 'is_significant_FlexMF', 'pvals_FlexMF', 'loadings_FlexMF',...
    'loadings_SeqNMF', 'W_hats_SeqNMF', 'H_hats_SeqNMF', "pvals_SeqNMF", "is_significant_SeqNMF", 'Ws', 'Hs')