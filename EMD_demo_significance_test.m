%% Demo script: FlexMF with EMD significance test
clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))
addpath(genpath(fullfile(root, 'FlexMF')));

%% Generate some sequences with temporal warping
number_of_seqences = 3;
T = 2000; % length of data to generate
Nneurons = 5*ones(number_of_seqences,1); % number of neurons in each sequence
Dt = 3.*ones(number_of_seqences,1); % gap between each member of the sequence
noise = 0; % probability of added noise in each bin
warp = 2; % stretch should be less than Dt
gap = 100;
neg = 0;
seed = 1;
[Xwarp, Wwarp, Hwarp, ~] = generate_data(T,Nneurons,Dt, 'noise',noise, 'warp', warp, 'seed', seed, 'len_burst', 1, 'dynamic', 0);

% Split into training and test set
Xtrain = Xwarp(:,1:round(T/2));
Xtest = Xwarp(:,1+round(T/2):end);

figure; SimpleWHPlot(Wwarp,Hwarp,'Data',Xwarp,'plotAll', 1, 'onsets', round(T/2)); 

title('generated data warping','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.5])

%% Find sequence with FlexMF in training data
L = size(Wwarp,3);
K = 3;
frob_norm = norm(Xtrain(:));
Xtrain = Xtrain/frob_norm*K;
Wwarp = Wwarp/frob_norm*K;

lambda = 1e-2;
lambda_M = 1e-2;
lambda_R = 1e2;
figure;
[What, Hhat_train, cost_train, errors_train, loadings, power, M_train, R_train] = FlexMF(Xtrain, 'K', K, 'L', L, ...
    'EMD',1, 'lambda', lambda, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'maxiter', 50);

figure; SimpleWHPlot(What, Hhat_train, 'plotAll', 1); title('FlexMF recon')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

%% Plot M, R
figure; plot_MR(M_train, R_train)

%% Fix What, rerun FlexMF on test data
frob_norm = norm(Xtest(:));
Xtest = Xtest/frob_norm*K;

figure;
[What, Hhat_test, cost_test, errors_test, ~, ~, M_test, R_test] = FlexMF(Xtest, 'K', K, 'L', L, 'W_fixed', 1, 'W_init', What,...
    'EMD',1, 'lambda', lambda, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'maxiter', 1);

figure; SimpleWHPlot(What, Hhat_test, 'plotAll', 1); title('FlexMF test recon')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; plot_MR(M_test, R_test)

%% test significance
[pvals,is_significant,is_single] = test_significance_EMD(Xtest, What, M_test, 'plot', 1);