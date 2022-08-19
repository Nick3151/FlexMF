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

%% Procedure for choosing lambda
nLambdas = 9; % increase if you're patient
nAlphas = 9;
K = 5; 
L = 50;
lambdas = sort([logspace(-1,-5,nLambdas)], 'ascend'); 
alphas = sort([logspace(-1,-5,nAlphas)], 'ascend'); 
loadings = zeros(nLambdas, nAlphas, K);
regularization = zeros(nLambdas, nAlphas);
costs = zeros(nLambdas, nAlphas); 
times = zeros(nLambdas, nAlphas); 
scores = zeros(nLambdas, nAlphas); 
W_hats = cell(nLambdas, nAlphas);
H_hats = cell(nLambdas, nAlphas);

for li = 1:length(lambdas)
    for ai = 1:length(alphas)
        tic
        [W_hat, H_hat, ~,loadings(li,ai,:),power]= FlexMF(X,'K',K,'L',L, 'maxiter', 50,...
            'lambdaL1W', .1, 'lambda', lambdas(li), 'alpha', alphas(ai), 'showPlot', 0); 

        [costs(li,ai),regularization(li,ai),~] = helper.get_FlexMF_cost(X,W_hat,H_hat);
        display(['Testing lambda ' num2str(li) '/' num2str(length(lambdas))])
        display(['Testing alpha ' num2str(ai) '/' num2str(length(alphas))])
        time = toc
        times(li,ai) = time;
        W_hats{li,ai} = W_hat;
        H_hats{li,ai} = H_hat;
        scores(li,ai) = helper.similarity(W, H, W_hat, H_hat);
    end
end
save('choose_lambda.mat', 'costs', 'regularization', 'times', 'loadings',...
    'lambdas', 'alphas', 'W_hats', 'H_hats', 'scores')