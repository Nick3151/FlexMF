clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
%% Generate some synthetic data
T = 30; % total number of trials
L = 50; % length of each trial
number_of_seqences = 3;
Nmotifs = 6*ones(number_of_seqences, 1); % the number of occurences of sequences
Nneurons = 5*ones(number_of_seqences, 1); % the number of neurons in each sequence
Dt = 3.*ones(number_of_seqences,1); % gap between each member of the sequence
noise = 0.001; % probability of added noise in each bin
jitter = 2*ones(number_of_seqences,1); % Jitter time std
participation = .9.*ones(number_of_seqences,1); % Participation probability = 100%
warp = 2; % the maximum warping time
bin = 0; % Binary data or not
neg = .2; % Proportion of negative indices in W
seed = 0;
[X, W, H, V_hat] = generate_data_trials(T, L, Nmotifs, Nneurons, Dt, noise, jitter, participation, warp, bin, neg, seed);
figure; SimpleWHPlot_trials(W, H, [], 1, X, 1); title('generated data','Fontsize',16)

set(gcf,'position',[200,200,1600,900])
