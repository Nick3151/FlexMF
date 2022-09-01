clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
%% Generate some synthetic data
number_of_seqences = 3;
T = 3000; % length of data to generate
Nneurons = 10*ones(number_of_seqences,1); % number of neurons in each sequence
Dt = 3.*ones(number_of_seqences,1); % gap between each member of the sequence
NeuronNoise = 0.001; % probability of added noise in each bin
SeqNoiseTime = zeros(number_of_seqences,1); % Jitter parameter = 0%
SeqNoiseNeuron = 1.*ones(number_of_seqences,1); % Participation parameter = 100%
neg = 0.2; % Proportion of negative indices in W
[X, W, H, V_hat] = generate_data(T,Nneurons,Dt,NeuronNoise,SeqNoiseTime,SeqNoiseNeuron,0,0,neg,1);
nuc_norm = norm(svd(X),1);
X = X/nuc_norm*size(X,1);

%% Choose lambda with BayesOpt
K = 5;
L = 50;
lambdaL1H = 0;
lambda_var = optimizableVariable('lambda', [1e-4, 1e-1], 'Transform', 'log');
lambdaW_var = optimizableVariable('lambdaL1W', [1e-3, 1], 'Transform', 'log');
alpha_var = optimizableVariable('alpha', [1e-6, 1e-2], 'Transform', 'log');
fun = @(x)compute_error_balance_score(X,K,L,x.lambda,lambdaL1H,x.lambdaL1W,x.alpha);
results = bayesopt(fun, [lambda_var,lambdaW_var,alpha_var],'AcquisitionFunctionName','expected-improvement-plus')