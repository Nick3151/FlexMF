clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
%% Generate some synthetic data, example
Trials = 30; % total number of trials
L = 50; % length of each trial
K = 3;
Nmotifs = 6*ones(K, 1); % the number of occurences of motifs
Nneurons = 5*ones(K, 1); % the number of neurons in each motif
Magnitudes = ones(K, 1); % the activation magnitudes of each motif
Dt = 3.*ones(K,1); % gap between each member of the motif
noise = 0.001; % probability of added noise in each bin
jitter = 2*ones(K,1); % Jitter time std
participation = .9.*ones(K,1); % Participation probability = 100%
warp = 0; % the maximum warping time
bin = 0; % Binary data or not
neg = 0; % Proportion of negative indices in W
seed = 0;
[X, W, H, V_hat] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Magnitudes, Dt, noise, jitter, participation, warp, bin, neg, seed);

N = size(W,1);
TestData = zeros(N,Trials*L);
for t=1:Trials
    TestData(:,(t-1)*L+1:t*L) = squeeze(X(:,:,t));
end

figure; SimpleWHPlot_trials(W, H, [], 1, X, 1); title('generated data','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
[pvals,is_significant] = test_significance_trials(TestData, Trials, L, W);

%% In how much percent of trials should the motif present in order to be significant
% Trials = 100;
Trials = 10;
L = 50; % length of each trial
presence_rate = 0.5:-0.05:0.05;
% Nmotifs = round(Trials*(presence_rate));
Nmotifs = 1:10;
K = length(Nmotifs);
Nneurons = 5*ones(K, 1); % the number of neurons in each motif
Magnitudes = ones(K, 1); % the activation magnitudes of each motif
Dt = 3.*ones(K,1); % gap between each member of the motif
noise = 0.1; % probability of added noise in each bin
jitter = 0*ones(K,1); % Jitter time std
participation = .9.*ones(K,1); % Participation probability = 100%
warp = 0; % the maximum warping time
bin = 1; % Binary data or not
neg = 0; % Proportion of negative indices in W
seed = 0;
[X, W, H, V_hat] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Magnitudes, Dt, noise, jitter, participation, warp, bin, neg, seed);

N = size(W,1);
TestData = zeros(N,Trials*L);
for t=1:Trials
    TestData(:,(t-1)*L+1:t*L) = squeeze(X(:,:,t));
end

figure; SimpleWHPlot_trials(W, H, [], 1, X, 1); title('generated data','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
[pvals,is_significant] = test_significance_trials(TestData, Trials, L, W);

figure;
plot(Nmotifs, pvals)
xlabel('# Motifs')
ylabel('p values')

%% The effect of magnitudes
Trials = 30;
L = 50; % length of each trial
K = 10;
Nmotifs = 6*ones(1,K);
Nneurons = 5*ones(K, 1); % the number of neurons in each motif
Magnitudes = 1:K; % the activation magnitudes of each motif
Dt = 3.*ones(K,1); % gap between each member of the motif
noise = 0.1; % probability of added noise in each bin
jitter = 0*ones(K,1); % Jitter time std
participation = .9.*ones(K,1); % Participation probability = 100%
warp = 0; % the maximum warping time
bin = 1; % Binary data or not
neg = 0; % Proportion of negative indices in W
seed = 0;
[X, W, H, V_hat] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Magnitudes, Dt, noise, jitter, participation, warp, bin, neg, seed);

N = size(W,1);
TestData = zeros(N,Trials*L);
for t=1:Trials
    TestData(:,(t-1)*L+1:t*L) = squeeze(X(:,:,t));
end

figure; SimpleWHPlot_trials(W, H, [], 1, X, 1); title('generated data','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
[pvals,is_significant] = test_significance_trials(TestData, Trials, L, W);

figure;
plot(Magnitudes, pvals)
xlabel('Magnitude of motif')
ylabel('p values')

%% The effect of noise
Trials = 30;
L = 50; % length of each trial
K = 5;
Nmotifs = 6*ones(1,K);
Nneurons = 5*ones(K, 1); % the number of neurons in each motif
Magnitudes = ones(1,K);% the activation magnitudes of each motif
Dt = 3.*ones(K,1); % gap between each member of the motif
jitter = 0*ones(K,1); % Jitter time std
participation = .9.*ones(K,1); % Participation probability = 100%
warp = 0; % the maximum warping time
bin = 1; % Binary data or not
neg = 0; % Proportion of negative indices in W
seed = 0;

noise = 0:0.1:0.5;
for i=1:length(noise)
    [X, W, H, V_hat] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Magnitudes, Dt, noise(i), jitter, participation, warp, bin, neg, seed);
    
    N = size(W,1);
    TestData = zeros(N,Trials*L);
    for t=1:Trials
        TestData(:,(t-1)*L+1:t*L) = squeeze(X(:,:,t));
    end
    [pvals(i,:),is_significant(i,:)] = test_significance_trials(TestData, Trials, L, W);
end

figure; SimpleWHPlot_trials(W, H, [], 1, X, 1); title('generated data','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

figure;
plot(noise, pvals(:,1))