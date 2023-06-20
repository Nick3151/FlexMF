clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))

%% Generate some synthetic data
conditions = ["clean"; "noise"; "participation"; "jitter"; "warp"; "overlap"];
n = str2double(getenv("SLURM_ARRAY_TASK_ID"));
condition = conditions(n);

nSim = 100;
rng(1)
seeds = randperm(1000, nSim);

number_of_seqences = 5;
T = 5000; % length of data to generate
Nneurons = 5*ones(number_of_seqences,1); % number of neurons in each sequence
Dt = 3.*ones(number_of_seqences,1); % gap between each member of the sequence
noise = .01; % probability of added noise in each bin
overlap_n = .5;
neg = 0;
participation = .7.*ones(number_of_seqences,1); % Participation parameter = 100%
warp = 2; % stretch should be less than Dt
jitter = 2.*ones(number_of_seqences,1);
gap = 200;

%% Procedure for choosing lambda
nLambdas = 9; % increase if you're patient
nAlphas = 9;
K = 5; 
lambdas = sort([logspace(0,-4,nLambdas)], 'ascend'); 
alphas = sort([logspace(-2,-6,nAlphas)], 'ascend'); 
loadings = cell(nSim, nLambdas, nAlphas);
reg_crosses = cell(nSim, nLambdas, nAlphas);
reg_Ws = cell(nSim, nLambdas, nAlphas);
reg_Hs = cell(nSim, nLambdas, nAlphas);
recon_errors = cell(nSim, nLambdas, nAlphas);
times = cell(nSim, nLambdas, nAlphas);
scores = cell(nSim, nLambdas, nAlphas);
num_detected = cell(nSim, nLambdas, nAlphas);
W_hats = cell(nSim, nLambdas, nAlphas);
H_hats = cell(nSim, nLambdas, nAlphas);
Ws = cell(nSim, 1);
Hs = cell(nSim, 1);

for li = 1:length(lambdas)
    for ai = 1:length(alphas)
        display(['Testing lambda ' num2str(li) '/' num2str(length(lambdas))])
        display(['Testing alpha ' num2str(ai) '/' num2str(length(alphas))])
        lambda = lambdas(li);
        alpha = alphas(ai);
        parfor n=1:nSim
            if strcmpi(condition, 'clean')
                [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'seed', seeds(n), 'gap', gap);
            elseif strcmpi(condition, 'noise')
                [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'seed', seeds(n), 'noise', noise, 'gap', gap);
            elseif strcmpi(condition, 'participation')
                [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'seed', seeds(n), 'participation', participation, 'gap', gap);
            elseif strcmpi(condition, 'jitter')
                [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'seed', seeds(n), 'jitter', jitter, 'gap', gap);
            elseif strcmpi(condition, 'warp')
                [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'seed', seeds(n), 'warp', warp, 'gap', gap);
            elseif strcmpi(condition, 'overlap')
                [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'seed', seeds(n), 'overlap_n', overlap_n, 'gap', gap);
            end
            nuc_norm = norm(svd(X),1);
            X = X/nuc_norm*size(X,1);
            Ws{n} = W;
            Hs{n} = H;
            L = size(W,3);
            tic
            [W_hat, H_hat, ~, ~, loadings{n,li,ai}, ~]= FlexMF(X,'K',K,'L',L, 'maxiter', 50,...
                'lambda', lambda, 'alpha', alpha, 'neg_prop', neg, 'showPlot', 0, 'verbal', 0); 
    
            [recon_errors{n,li,ai},reg_crosses{n,li,ai},reg_Ws{n,li,ai},reg_Hs{n,li,ai}] = helper.get_FlexMF_cost(X,W_hat,H_hat);
            
            time = toc
            times{n,li,ai} = time;
            W_hats{n,li,ai} = W_hat;
            H_hats{n,li,ai} = H_hat;
            [coeff, ids] = helper.similarity_W(W, W_hat);
            scores{n,li,ai} = coeff;
            num_detected{n,li,ai} = length(ids);
        end
    end
end

save(strcat('choose_lambda_', condition, '.mat'), 'recon_errors', 'reg_crosses', 'reg_Ws', 'reg_Hs', 'times', 'loadings',...
    'lambdas', 'alphas', 'W_hats', 'H_hats', 'scores', "num_detected", "Ws", "Hs")