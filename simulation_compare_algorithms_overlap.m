clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(fullfile(root, 'Utils'))
addpath(fullfile(root, 'MATLAB-tools'))
addpath(genpath(fullfile(root, 'FlexMF')));

%% Impact of shared neuron
disp('Impact of shared neurons on results')
Trials = 200;
L = 50; % length of each trial
K = 10;
Nmotifs = 2*(1:K);
Nneurons = 5*ones(K, 1); % the number of neurons in each motif
Dt = 3.*ones(K,1); % gap between each member of the motif

overlaps_n = 0:0.2:0.8; % Participation probability

nSim = 100;
seeds = randperm(1000, nSim);
pvals_SeqNMF = cell(nSim,length(overlaps_n));
is_significants_SeqNMF = cell(nSim,length(overlaps_n));
loadings_SeqNMF = cell(nSim,length(overlaps_n));
W_hats_SeqNMF = cell(nSim,length(overlaps_n));
H_hats_SeqNMF = cell(nSim,length(overlaps_n));
pvals_FlexMF = pvals_SeqNMF;
is_significants_FlexMF = is_significants_SeqNMF;
loadings_FlexMF = loadings_SeqNMF;
W_hats_FlexMF = W_hats_SeqNMF;
H_hats_FlexMF = H_hats_SeqNMF;
Ws = W_hats_SeqNMF;
Hs = H_hats_SeqNMF;
TrainingDatas = W_hats_SeqNMF;

for j=1:length(overlaps_n)
    overlap_n = overlaps_n(j);
    disp(['neurons overlap rate=', num2str(overlaps_n(j))])
    parfor i=1:nSim
        tic
        [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, ...
            'overlap_n', overlap_n, 'seed', seeds(i));
        groups = zeros(Trials,1);
        for k=1:K
            groups(motif_ind{k}) = k;
        end
        Ws{i,j} = W;
        Hs{i,j} = H;
        
        % Dimension N*L*Trials
        cv = cvpartition(groups, "KFold",2);
        ind_train = find(training(cv,1));
        X_train = X(:,:,ind_train);
        ind_test = find(test(cv,1));
        X_test = X(:,:,ind_test);
        
        % Dimension N*T
        N = size(W,1);
        TrainingData = zeros(N,cv.TrainSize(1)*L);
        for t=1:cv.TrainSize(1)
            TrainingData(:,(t-1)*L+1:t*L) = squeeze(X_train(:,:,t));
        end
        TestData = zeros(N,cv.TestSize(1)*L);
        for t=1:cv.TestSize(1)
            TestData(:,(t-1)*L+1:t*L) = squeeze(X_test(:,:,t));
        end
        
        % Normalize training data
        nuc_norm = norm(svd(TrainingData),1);
        TrainingData = TrainingData/nuc_norm*N;
    
        lambda = .005;
        alpha = 5e-5;
    
        % Run SeqNMF with multiplication rule
        [W_hat, H_hat, ~,loadings_SeqNMF{i,j},power]= seqNMF(TrainingData,'K',K,'L',L,...
                'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
        p = .05;
        [pvals_SeqNMF{i,j},is_significants_SeqNMF{i,j}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
        W_hats_SeqNMF{i,j} = W_hat;
        H_hats_SeqNMF{i,j} = H_hat;
        display(['SeqNMF run ' num2str(i) '/' num2str(nSim)])
    
        % Run SeqNMF with Bregman Iteration
        [W_hat, H_hat, ~,loadings_FlexMF{i,j},power]= FlexMF(TrainingData,'K',K,'L',L,...
                'lambda', lambda, 'alpha', alpha, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0, 'verbal', 0); 
        p = .05;
        [pvals_FlexMF{i,j},is_significants_FlexMF{i,j}] = test_significance_trials(TestData, cv.TestSize(1), L, W_hat,[],p);
        W_hats_FlexMF{i,j} = W_hat;
        H_hats_FlexMF{i,j} = H_hat;
        TrainingDatas{i,j} = TrainingData;
        display(['FlexMF run ' num2str(i) '/' num2str(nSim)])
        toc
    end
end

save('simulate_results_overlap.mat', "H_hats_FlexMF", 'W_hats_FlexMF', 'is_significants_FlexMF', 'pvals_FlexMF', 'loadings_FlexMF',...
    'loadings_SeqNMF', 'W_hats_SeqNMF', 'H_hats_SeqNMF', "pvals_SeqNMF", "is_significants_SeqNMF", 'Ws', 'Hs', 'TrainingDatas')

%% Analysis
% load('simulate_results_overlap.mat')
% i = 1;
% TrainingData = TrainingDatas{i,3};
% W_hat = W_hats_SeqNMF{i,3};
% H_hat = H_hats_SeqNMF{i,3};
% lambda=.005;
% [recon_error_SeqNMF, reg_cross, reg_W, reg_H] = helper.get_FlexMF_cost(TrainingData,W_hat,H_hat);
% reg_cross_SeqNMF = reg_cross*lambda;
% 
% is_significant = is_significants_SeqNMF{i,3};
% 
% plotAll = 1;
% figure; SimpleWHPlot_patch(W_hat, H_hat, 100, 50, is_significant, [], plotAll); title('SeqNMF reconstruction')
% set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
% figure; SimpleWHPlot_patch(W_hat, H_hat, 100, 50, is_significant, TrainingData, plotAll); title('SeqNMF factors, with raw data')
% set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])