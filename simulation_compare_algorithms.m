% Simulate data with trials, compare the performance of SeqNMF with
% multiplication update rule and Split Bregman Iteration

clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(fullfile(root, 'Utils'))
addpath(fullfile(root, 'MATLAB-tools'))
rmpath(genpath(fullfile(root, 'seqNMF-master')));
addpath(genpath(fullfile(root, 'FlexMF')));

%% Simulate data
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
len_spike = 20; % Continuous firing time
dynamic = 1; % Consider calcium dynamics or not
overlap = 0;
neg = 0; % Proportion of negative indices in W

nSim = 100;
seeds = randperm(1000, nSim);
pvals_SeqNMF = ones(nSim,K);
is_significant_SeqNMF = zeros(nSim,K);
loadings_SeqNMF = zeros(nSim,K);
W_hats_SeqNMF = cell(nSim,1);
H_hats_SeqNMF = cell(nSim,1);
pvals_FlexMF = pvals_SeqNMF;
is_significant_FlexMF = is_significant_SeqNMF;
loadings_FlexMF = loadings_SeqNMF;
W_hats_FlexMF = W_hats_SeqNMF;
H_hats_FlexMF = H_hats_SeqNMF;
Ws = W_hats_SeqNMF;
Hs = H_hats_SeqNMF;

for i=1:nSim
    [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Magnitudes, Dt, noise, jitter, participation, warp, len_spike, dynamic, overlap, neg, seeds(i));
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
    [W_hat, H_hat, ~,loadings_SeqNMF,power]= seqNMF(TrainingData,'K',K,'L',L,...
            'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
    p = .05;
    [pvals_SeqNMF(i,:),is_significant_SeqNMF(i,:)] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
    W_hats_SeqNMF{i} = W_hat;
    H_hats_SeqNMF{i} = H_hat;
    display(['SeqNMF run ' num2str(i) '/' num2str(nSim)])

    % Run SeqNMF with Bregman Iteration
    [W_hat, H_hat, ~,loadings_FlexMF(i,:),power]= FlexMF(TrainingData,'K',K,'L',L,...
            'lambda', lambda, 'alpha', alpha, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0); 
    p = .05;
    [pvals_FlexMF(i,:),is_significant_FlexMF(i,:)] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
    W_hats_FlexMF{i} = W_hat;
    H_hats_FlexMF{i} = H_hat;
    display(['FlexMF run ' num2str(i) '/' num2str(nSim)])
end

save('simulate_results.mat', "H_hats_FlexMF", 'W_hats_FlexMF', 'is_significant_FlexMF', 'pvals_FlexMF', 'loadings_FlexMF',...
    'loadings_SeqNMF', 'W_hats_SeqNMF', 'H_hats_SeqNMF', "pvals_SeqNMF", "is_significant_SeqNMF", 'Ws', 'Hs')

%% Compare results
load('simulate_results.mat')

i = 1;
W = Ws{i};
W_hat_SeqNMF = W_hats_SeqNMF{i};
W_hat_FlexMF = W_hats_FlexMF{i};
[coeff_SeqNMF, ids_SeqNMF] = helper.similarity_W(W, W_hat_SeqNMF);
[coeff_FlexMF, ids_FlexMF] = helper.similarity_W(W, W_hat_FlexMF);