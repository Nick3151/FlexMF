% Compare Matrix Multiplication Update rule vs Split Bregman Iteration
clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(fullfile(root, 'Utils'))
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))
rmpath(genpath(fullfile(root, 'seqNMF-master')));
addpath(genpath(fullfile(root, 'FlexMF')));

%%
Trials = 200;
% Trials = 10;
L = 50; % length of each trial

K = 10;
Nmotifs = 2*(K:-1:1);
Nneurons = 5*ones(K, 1); % the number of neurons in each motif
Dt = 3.*ones(K,1); % gap between each member of the motif
noise = 0; % probability of added noise in each bin
participation = .7.*ones(K,1); % Participation probability = 100%
warp = 0; % the maximum warping time
len_spike = 1; % Continuous firing time
dynamic = 0; % Consider calcium dynamics or not
overlap_t = 0;
overlap_n = .6;
neg = 0; % Proportion of negative indices in W

nsim = 100;
seeds = randperm(1000, nsim);
pvals = zeros(nsim,K);
is_significant = zeros(nsim,K);

[X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, ...
    'overlap_n', overlap_n);
% [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, ...
%     'participation', participation);
groups = zeros(Trials,1);
for k=1:K
    groups(motif_ind{k}) = k;
end

% Dimension N*L*Trials
cv = cvpartition(groups, "KFold",2);
ind_train = training(cv,1);
X_train = X(:,:,ind_train);
ind_test = test(cv,1);
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
f1 = figure;
SimpleXplot_patch([TrainingData, TestData], [cv.TrainSize(1), cv.TestSize(1)], L); 
set(f1,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

f2 = figure;
H_train = H(:,ind_train);
SimpleWHPlot_trials(W, H_train, [], X_train, 1); title('generated data','Fontsize',16)
set(f2,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

save2pdf('Simulated_data_overlap.pdf', f2)
% save2pdf('Simulated_data_participation.pdf', f2)

% Normalize training data
nuc_norm = norm(svd(TrainingData),1);
TrainingData = TrainingData/nuc_norm*N;

%% Procedure for choosing lambda
nLambdas = 20; % increase if you're patient
lambdaL1H = 0;
lambdaL1W = 0;
lambdaOrthoH = 0;
lambdas = sort(logspace(-1,-5,nLambdas), 'ascend'); 
loadings = [];
regularization = []; 
cost = []; 
for li = 1:length(lambdas)
    [N,T] = size(TrainingData);
    [What, Hhat, ~,~,loadings(li,:),power]= seqNMF(TrainingData,'K',K,'L',L,...
        'lambdaL1W', lambdaL1W, 'lambda', lambdas(li), 'lambdaOrthoH', lambdaOrthoH, 'maxiter', 100, 'showPlot', 0); 
    [cost(li),regularization(li),~] = helper.get_seqNMF_cost(TrainingData,What,Hhat);
    display(['Testing lambda ' num2str(li) '/' num2str(length(lambdas))])
end
%% plot costs as a function of lambda
windowSize = 3; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
Rs = filtfilt(b,a,regularization); 
minRs = prctile(regularization,10); maxRs= prctile(regularization,90);
Rs = (Rs-minRs)/(maxRs-minRs); 
R = (regularization-minRs)/(maxRs-minRs); 
Cs = filtfilt(b,a,cost); 
minCs =  prctile(cost,10); maxCs =  prctile(cost,90); 
Cs = (Cs -minCs)/(maxCs-minCs); 
C = (cost -minCs)/(maxCs-minCs); 

figure; hold on
plot(lambdas,Rs, 'b')
plot(lambdas,Cs,'r')
scatter(lambdas, R, 'b', 'markerfacecolor', 'flat');
scatter(lambdas, C, 'r', 'markerfacecolor', 'flat');
xlabel('Lambda'); ylabel('Cost (au)')
set(legend('Correlation cost', 'Reconstruction cost'), 'Box', 'on')
set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])

%% Run SeqNMF
lambda = .005;
% lambda = .01;
lambdaL1H = 0;
lambdaL1W = 0;
lambdaOrthoH = 0;

figure;
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
[What, Hhat, ~, errors_SeqNMF,loadings,power]= seqNMF(TrainingData,'K',K,'L',L,...
            'lambdaL1W', lambdaL1W, 'lambda', lambda, 'lambdaOrthoH', lambdaOrthoH, 'maxiter', 100, 'showPlot', 1); 

[recon_error_SeqNMF, reg_cross, reg_W, reg_H] = helper.get_FlexMF_cost(TrainingData,What,Hhat);
reg_cross_SeqNMF = reg_cross*lambda;
reg_W_SeqNMF = reg_W*lambdaL1W;
reg_H_SeqNMF = reg_H*lambdaL1H;

p = .05; % desired p value for factors
display('Testing significance of factors on held-out data')
[pvals,is_significant] = test_significance_trials(TestData, cv.TestSize(1), L, What,[],p);

% plot, sorting neurons by latency within each factor
[max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(What(:,:,:),1);
indSort = hybrid(:,3);

%% Look at factors
plotAll = 1;
figure; SimpleWHPlot_patch(What, Hhat, cv.TrainSize(1), L, is_significant, [], plotAll); title('SeqNMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot_patch(What, Hhat, cv.TrainSize(1), L, is_significant, TrainingData, plotAll); title('SeqNMF factors, with raw data')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

% save2pdf('Simulated_data_SeqNMF.pdf', gcf)
% save2pdf('Simulated_data_HOrth.pdf', gcf)

% Compute similarity to ground truth
[coeff_SeqNMF, ids_SeqNMF] = helper.similarity_W(W, What);

%% Run FlexMF
% alpha = 1e-4;
alpha = 5e-5;

% load([exp_num, '_BayesOpt.mat'])
% lambda = results.XAtMinEstimatedObjective.lambda;
% alpha = results.XAtMinEstimatedObjective.alpha;

lambdaL1H = 0;
lambdaL1W = 0;

display('Running FlexMF on 2p data')
figure;
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
tic
[What,Hhat,~,errors_FlexMF,loadings,power] = FlexMF(TrainingData,'K',K, 'L', L, 'maxiter', 50,...
    'lambda', lambda, 'alpha', alpha, 'lambdaL1W', lambdaL1W, 'lambdaL1H', lambdaL1H, 'neg_prop', 0, 'showPlot', 1);
toc

% 
figure;
hold on
plot(errors_SeqNMF(2:length(errors_FlexMF), 1), 'b')
plot(errors_SeqNMF(2:length(errors_FlexMF), 2)*lambda, 'b--')
plot(errors_FlexMF(2:end, 1), 'r')
plot(errors_FlexMF(2:end, 2)*lambda, 'r--')
xlabel('# Iteration')
legend({'MUR Reconstruction Error', 'MUR Regularizaion Error', ...
    'SBI Reconstruction Error', 'SBI Regularizaion Error'})
save2pdf('Simulation_Error_Curves.pdf')


[recon_error_FlexMF, reg_cross, reg_W, reg_H] = helper.get_FlexMF_cost(TrainingData,What,Hhat);
reg_cross_FlexMF = reg_cross*lambda;
reg_W_FlexMF = reg_W*lambdaL1W;
reg_H_FlexMF = reg_H*lambdaL1H;

sparseness = sum(Hhat>0, 2)/size(Hhat,2);

p = .05; % desired p value for factors
display('Testing significance of factors')
[pvals,is_significant] = test_significance_trials(TestData, cv.TestSize(1), L, What,[],p);

%% Look at factors
plotAll = 1;
figure; SimpleWHPlot_patch(What, Hhat, cv.TrainSize(1), L, is_significant, [], plotAll); title('FlexMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot_patch(What, Hhat, cv.TrainSize(1), L, is_significant, TrainingData, plotAll); title('FlexMF factors, with raw data')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

% save2pdf('Simulated_data_FlexMF.pdf', gcf)

% Compute similarity to ground truth
[coeff_FlexMF, ids_FlexMF] = helper.similarity_W(W, What);

%% Compare algorithm results
coeff_all = zeros(2,K);
coeff_all(1, ids_SeqNMF) = coeff_SeqNMF;
coeff_all(2, ids_FlexMF) = coeff_FlexMF;

figure; bar(1:K, coeff_all);
xlabel('# Motif Occurences')
ylabel('Correlation to ground truth')
legend({'SeqNMF', 'FlexMF'}, 'Location', 'northwest')
set(gca, 'FontSize', 14)