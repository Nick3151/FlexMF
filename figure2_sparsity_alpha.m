clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))

%% Generate some synthetic data
nLambdas = 5;
lambdas = sort([logspace(-1,-3,nLambdas)], 'ascend'); 
n = str2double(getenv("SLURM_ARRAY_TASK_ID"));
lambda = lambdas(n);

nSim = 100;
rng(1)
seeds = randperm(1000, nSim);

number_of_seqences = 5;
T = 5000; % length of data to generate
Nneurons = 5*ones(number_of_seqences,1); % number of neurons in each sequence
Dt = 3.*ones(number_of_seqences,1); % gap between each member of the sequence
neg = 0;
gap = 200;

%% Run simulation on different combinations of alpha_W, alpha_H
nHs = 9;
nWs = 9;
K = 5;
alpha_Ws = sort([logspace(-2,-6,nWs)], 'ascend');
alpha_Hs = sort([logspace(-2,-6,nHs)], 'ascend');
loadings = cell(nSim, nWs, nHs);
reg_crosses = cell(nSim, nWs, nHs);
recon_errors = cell(nSim, nWs, nHs);
times = cell(nSim, nWs, nHs);
scores_W = cell(nSim, nWs, nHs);
scores_H = cell(nSim, nWs, nHs);
num_detected = cell(nSim, nWs, nHs);
num_success = cell(nSim, nWs, nHs);
sparsity_Ws = cell(nSim, nWs, nHs);
sparsity_Hs = cell(nSim, nWs, nHs);

W_hats = cell(nSim, nWs, nHs);
H_hats = cell(nSim, nWs, nHs);
Ws = cell(nSim, 1);
Hs = cell(nSim, 1);


for Wi = 1:length(alpha_Ws)
    for Hi = 1:length(alpha_Hs)
        display(['Testing alpha_W ' num2str(Wi) '/' num2str(length(alpha_Ws))])
        display(['Testing alpha_H ' num2str(Hi) '/' num2str(length(alpha_Hs))])
        alpha_W = alpha_Ws(Wi);
        alpha_H = alpha_Hs(Hi);
        for n=1:nSim
            [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'seed', seeds(n), 'gap', gap);
            nuc_norm = norm(svd(X),1);
            X = X/nuc_norm*size(X,1);
            Ws{n} = W;
            Hs{n} = H;
            L = size(W,3);
            tic
            [W_hat, H_hat, ~, ~, loadings{n,Wi,Hi}, ~]= FlexMF(X,'K',K,'L',L, 'maxiter', 50,...
                'lambda', lambda, 'alpha_W', alpha_W, 'alpha_H', alpha_H, 'neg_prop', neg, 'showPlot', 0, 'verbal', 0); 
    
            [recon_errors{n,Wi,Hi},reg_crosses{n,Wi,Hi},~,~] = helper.get_FlexMF_cost(X,W_hat,H_hat);
            
            time = toc
            times{n,Wi,Hi} = time;
            W_hats{n,Wi,Hi} = W_hat;
            H_hats{n,Wi,Hi} = H_hat;
            sparsity_Ws{n,Wi,Hi} = mean(W_hat==0, 'all');
            sparsity_Hs{n,Wi,Hi} = mean(H_hat==0, 'all');

            [coeffs_W, coeffs_H, ids] = helper.similarity_WH(W, H, W_hat, H_hat);
            scores_W{n,Wi,Hi} = coeffs_W;
            scores_H{n,Wi,Hi} = coeffs_H;
            num_detected{n,Wi,Hi} = length(ids);
            success_thresh = .8;
            ids_success = ids(coeffs_H>success_thresh & coeffs_W>success_thresh);
            num_success{n,Wi,Hi} = length(ids_success);
        end
    end
end

save(sprintf('sparsity_alpha_lambda=%0.3e.mat', lambda), 'recon_errors', 'reg_crosses', 'times', 'loadings',...
    'alpha_Ws', 'alpha_Hs', 'lambda', 'W_hats', 'H_hats', 'scores_W', 'scores_H',...
    "num_detected", 'num_success', 'Ws', 'Hs', 'sparsity_Ws', 'sparsity_Hs')