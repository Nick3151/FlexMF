%% Test FlexMF parameters (sacling of lambda, lambda_R, lambda_M)
clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))
rmpath(genpath(fullfile(root, 'seqNMF-master')));
addpath(genpath(fullfile(root, 'FlexMF')));

%% Generate some synthetic data with temporal jittering or time warping
number_of_seqences = 3;
T = 800; % length of data to generate
Nneurons = 10*ones(number_of_seqences,1); % number of neurons in each sequence
Dt = 3.*ones(number_of_seqences,1); % gap between each member of the sequence
noise = 0; % probability of added noise in each bin
jitter = 5*ones(number_of_seqences,1); % Jitter std
participation = 1.*ones(number_of_seqences,1); % Participation parameter = 100%
warp = 2; % stretch should be less than Dt
gap = 100;
neg = 0;
seed = 1;
[X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',noise, 'seed', seed, 'len_burst', 1, 'dynamic', 0);
[Xwarp, Wwarp, Hwarp, ~] = generate_data(T,Nneurons,Dt, 'noise',noise, 'warp', warp, 'seed', seed, 'len_burst', 1, 'dynamic', 0);
[Xjit, Wjit, Hjit, ~] = generate_data(T,Nneurons,Dt, 'noise',noise, 'jitter', jitter, 'seed', seed, 'len_burst', 1, 'dynamic', 0);
L = size(W,3);

%% Normalize data
K = 3;
% frob_norm = norm(Xwarp(:));
% Xwarp = Xwarp/frob_norm*K;
% Wwarp = Wwarp/frob_norm*K;

%% Run SeqNMF
lambda = .05;
lambdaL1H = 0;
lambdaL1W = 0;
lambdaOrthoH = 0;

figure;
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
[What_SeqNMF, Hhat_SeqNMF, ~, errors_SeqNMF,loadings,power]= seqNMF(Xwarp,'K',K,'L',L,...
            'lambda', lambda, 'maxiter', 50, 'showPlot', 1); 

%% Run FlexMF with different parameters
nlambdas = 5;
lambdas = logspace(-4,0,nlambdas);
lambda_Ms = lambdas;
lambda_R = 1;

nSim = 10;
Whats = cell(nlambdas,nSim);
Hhats = cell(nlambdas,nSim);
emds_W = cell(nlambdas, nSim);
emds_H = cell(nlambdas, nSim);
num_detected = cell(nlambdas, nSim);
ids_match = cell(nlambdas, nSim);

for i=1:length(lambdas)
    display(['Testing lambda ' num2str(i) '/' num2str(nlambdas)])
    lambda = lambdas(i);
    lambda_M = lambda_Ms(i);
    parfor n=1:nSim
        display(n)
        tic
        [What, Hhat, cost, errors, loadings, power, M, R] = FlexMF(Xwarp, 'K', K, 'L', L, ...
            'EMD',1, 'lambda', lambda, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'maxiter', 50, 'tolerance', 1e-4, ...
            'W_init', What_SeqNMF, 'H_init', Hhat_SeqNMF, 'showPlot', 0, 'verbal', 0);
        toc
        Whats{i,n} = What;
        Hhats{i,n} = Hhat;

        disp('Evaluate EMDs of results')
        tic
        [emds_W{i,n}, emds_H{i,n}, ids] = helper.similarity_WH_EMD(Wwarp, Hwarp, What, Hhat);
        num_detected{i,n} = length(ids);
        ids_match{i,n} = ids;
        toc
    end
end
save('Simulation_Results/EMD_params_test1.mat', "Hhats", "Whats", "lambda_R", "lambda_Ms", "lambdas", "Xwarp")

%% Look at factors
% clear all
load('Simulation_Results/EMD_params_test1.mat')
i = 3;
n = 6;
plotAll = 1;
figure; SimpleWHPlot_patch(Whats{i,n}, Hhats{i,n}, 'plotAll', plotAll); title('FlexMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])