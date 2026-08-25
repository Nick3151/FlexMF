clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))

data_type = 'warp';

profile on

%% Generate some synthetic data, with warping
arrayn = 4;
nSim = 1;
rng(1)
seeds = randperm(1000, nSim);

K = 3;
T = 2000; % length of data to generate
Nneurons = 5*ones(K,1); % number of neurons in each sequence
Dt = 3.*ones(K,1); % gap between each member of the sequence
neg = 0;
noise = .01;
participation = .7.*ones(K,1); 
jitter = 5*ones(K,1);
warp = 2;
gap = 100;

%% Run simulation on different combinations of lambda_M, lambda_R
nLambdas = 9;
nMs = 9;
nRs = 9;

lambdas = logspace(-4,0,nLambdas);
lambda_Ms = logspace(-4,0,nMs);
lambda_Rs = logspace(0,4,nRs);

times = cell(nSim, 1);
emds_W = cell(nSim, 1);
emds_H = cell(nSim, 1);
num_detected = cell(nSim, 1);
W_hats = cell(nSim, 1);
H_hats = cell(nSim, 1);
Ws = cell(nSim, 1);
Hs = cell(nSim, 1);

Li = floor(arrayn/(nMs*nRs))+1;
Mi = floor((arrayn - (Li-1)*nMs*nRs)/nRs)+1;
Ri = arrayn - (Li-1)*nMs*nRs - (Mi-1)*nRs;
display(['Testing lambda ' num2str(Li) '/' num2str(nLambdas)])
display(['Testing lambda_M ' num2str(Mi) '/' num2str(nMs)])
display(['Testing lambda_R ' num2str(Ri) '/' num2str(nRs)])
lambda = lambdas(Li);
lambda_M = lambda_Ms(Mi);
lambda_R = lambda_Rs(Ri);
for n=1:nSim
    switch data_type
        case 'warp'
            [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',0, 'warp', warp, 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
        case 'jitter'
            [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',0, 'jitter', jitter, 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
        case 'participation'
            [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',0, 'participation', participation, 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
        case 'noise'
            [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise', noise, 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
    end
    frob_norm = norm(X(:));
    X = X/frob_norm*K;
    W = W/frob_norm*K;
    Ws{n} = W;
    Hs{n} = H;
    L = size(W,3);
    tic
    [W_hat, H_hat]= FlexMF(X,'K',K,'L',L, 'EMD',1, 'maxiter', 50,...
        'lambda', lambda, 'lambda_M', lambda_M, 'lambda_R', lambda_R, 'showPlot', 0, 'verbal', 0); 
    
    time = toc
    times{n} = time;
    W_hats{n} = W_hat;
    H_hats{n} = H_hat;
    
    disp('Evaluate EMDs of results')
    tic
    [emds_W{n}, emds_H{n}, ids] = helper.similarity_WH_EMD(W, H, W_hat, H_hat);
    num_detected{n} = length(ids);
    toc
end

profile viewer
profsave
p = profile('info')
