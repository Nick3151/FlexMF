clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))

%% Generate some synthetic data
number_of_seqences = 5;
T = 5000; % length of data to generate
Nneurons = 5*ones(number_of_seqences,1); % number of neurons in each sequence
Dt = 3.*ones(number_of_seqences,1); % gap between each member of the sequence
noise = 0.001; % probability of added noise in each bin
jitter = 5*ones(number_of_seqences,1); % Jitter std
participation = 1.*ones(number_of_seqences,1); % Participation parameter = 100%
warp = 2; % stretch should be less than Dt
gap = 200;
neg = 0;
bin = 0;
seed = 1;
[X, W, H, ~] = generate_data(T,Nneurons,Dt, 'seed', seed, 'gap', gap, 'overlap_n', .5);

plotAll = 1;
figure; SimpleWHPlot_patch(W,H,[],[],[],X,plotAll); title('generated data raw','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

nuc_norm = norm(svd(X),1);
X = X/nuc_norm*size(X,1);

%% Choose lambda with BayesOpt
K = 5;
L = 50;
lambdaL1H = 0;
lambda_var = optimizableVariable('lambda', [1e-4, 1e0], 'Transform', 'log');
lambdaW_var = optimizableVariable('lambdaL1W', [1e-3, 1], 'Transform', 'log');
alpha_var = optimizableVariable('alpha', [1e-6, 1e-2], 'Transform', 'log');

% fun = @(x)compute_error_balance_score_3d(X,K,L,x.lambda,x.lambdaL1W,x.alpha);
% results = bayesopt(fun, [lambda_var,lambdaW_var,alpha_var],'AcquisitionFunctionName','expected-improvement-plus','UseParallel',true, 'MaxObjectiveEvaluations',100)

etas = [1e-6];
for i=1:length(etas)
    fun = @(x)compute_error_balance_score_2d(X,K,L,x.lambda,x.alpha, 'eta', etas(i));
    results{i} = bayesopt(fun, [lambda_var,alpha_var],'AcquisitionFunctionName','expected-improvement-plus','UseParallel',true, 'MaxObjectiveEvaluations',100)
end

lambdas = cellfun(@(x) x.XAtMinEstimatedObjective.lambda, results);
alphas = cellfun(@(x) x.XAtMinEstimatedObjective.alpha, results);

%% Run with best parameters
lambda = results{1}.XAtMinEstimatedObjective.lambda;
alpha = results{1}.XAtMinEstimatedObjective.alpha;
[W_hat,H_hat,~,~,loadings,power] = FlexMF(X,'K',K, 'L', L, 'maxiter', 50,...
    'lambda', lambda, 'alpha', alpha, 'lambdaL1W', 0, 'lambdaL1H', 0, 'neg_prop', 0, 'showPlot', 0);
[recon_error, reg_cross, reg_W, reg_H] = helper.get_FlexMF_cost(X,W_hat,H_hat);
reg_cross = reg_cross*lambda;

plotAll = 1;
figure; SimpleWHPlot_patch(W_hat,H_hat,[],[],[],[],plotAll); title('generated data raw','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])