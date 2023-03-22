% Compare Matrix Multiplication Update rule vs Split Bregman Iteration
clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(fullfile(root, 'Utils'))
addpath(fullfile(root, 'MATLAB-tools'))

%%
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

nsim = 100;
seeds = randperm(1000, nsim);
pvals = zeros(nsim,K);
is_significant = zeros(nsim,K);

[X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Magnitudes, Dt, noise, jitter, participation, warp, len_spike, dynamic, overlap, neg, seeds(1));
groups = zeros(Trials,1);
for k=1:K
    groups(motif_ind{k}) = k;
end

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
SimpleXplot_patch([TrainingData, TestData], [cv.TrainSize(1), cv.TestSize(1)], L); 

set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])