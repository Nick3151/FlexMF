clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))

data_type = 'warp';

% Start parallel pool
disp('Starting a parallel pool...');
parpool('local', 48);  % 2 nodes × 12 tasks

%% Generate some synthetic data, with warping
nLambdas = 9;
lambdas = sort(logspace(0,-4,nLambdas), 'ascend'); 
arrayn = str2double(getenv("SLURM_ARRAY_TASK_ID"))+1
%n=1;
lambda = lambdas(arrayn);

nSim = 48;
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
nMs = 9;
nRs = 9;
K = 3;
lambda_Ms = logspace(-4,0,nMs);
lambda_Rs = logspace(0,4,nRs);

tmp_Ms = repmat(lambda_Ms.',1,nRs);
tmp_Rs = repmat(lambda_Rs,nMs,1);

tmp=cell(nMs*nRs,2);
tmp=[tmp_Ms(:) tmp_Rs(:)];

times = cell(nSim, nRs, nMs);
emds_W = cell(nSim, nRs, nMs);
emds_H = cell(nSim, nRs, nMs);
% num_detected = cell(nSim, nRs, nMs);
W_hats = cell(nSim, nRs, nMs);
H_hats = cell(nSim, nRs, nMs);
Ws = cell(nSim, 1);
Hs = cell(nSim, 1);


%for Mi = 1:length(lambda_Ms)
%    for Ri = 1:length(lambda_Rs)

arrayn
Mi = floor(arrayn/nMs)+1
Ri = arrayn - nRs*(Mi-1)

        display(['Testing lambda_M ' num2str(Mi) '/' num2str(length(lambda_Ms))])
        display(['Testing lambda_R ' num2str(Ri) '/' num2str(length(lambda_Rs))])
        lambda_M = lambda_Ms(Mi);
        lambda_R = lambda_Rs(Ri);

        parfor n=1:nSim

	    n 
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
            times{n,Mi,Ri} = time;
            W_hats{n,Mi,Ri} = W_hat;
            H_hats{n,Mi,Ri} = H_hat;
            
%            disp('Evaluate EMDs of results')
%            tic
%            [emds_W{n,Wi,Hi}, emds_H{n,Wi,Hi}, ids] = helper.similarity_WH_EMD(W, H, W_hat, H_hat);
%            num_detected{n,Wi,Hi} = length(ids);
%            toc
%        end
%    end
end

save(sprintf('testarray_%s_lambda=%0.3e.mat', data_type, lambda), 'times',...
    'lambda_Ms', 'lambda_Rs', 'lambda', 'W_hats', 'H_hats', 'emds_W', 'emds_H',...
    'Ws', 'Hs')

% Shut down the parallel pool
delete(gcp('nocreate'))

exit
