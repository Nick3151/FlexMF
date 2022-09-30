clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
%% Generate some synthetic data
trials = 30; % total number of trials
frames = 150; % length of each trial
number_of_seqences = 3;
Nsequences = 6*ones(number_of_seqences, 1); % the number of occurences of sequences
Nneurons = 10*ones(number_of_seqences, 1); % the number of neurons in each sequence
Dt = 3.*ones(number_of_seqences,1); % gap between each member of the sequence
noise = 0.001; % probability of added noise in each bin
jitter = 2*ones(number_of_seqences,1); % Jitter time std
participation = .9.*ones(number_of_seqences,1); % Participation probability = 100%
warp = 2; % the maximum warping time
bin = 0; % Binary data or not
neg = .2; % Proportion of negative indices in W
seed = 2;
[X, W, H, V_hat] = generate_data_trials(trials, frames, Nsequences, Nneurons, Dt, noise, jitter, participation, warp, bin, neg, seed);
figure; SimpleWHPlot(W, H, trials, frames, 1, X); title('generated data','Fontsize',16)

set(gcf,'position',[200,200,1600,900])
nuc_norm = norm(svd(X),1);
X = X/nuc_norm*size(X,1);