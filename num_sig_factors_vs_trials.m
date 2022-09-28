% How many trials needed to find significant sequence
% Running on simluated data
% Choose best parameter for each sequence gap, using BayesOpt
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
gap = 200; % maximum gap between sequences
neg = 0; % Proportion of negative indices in W
[X, W, H, V_hat] = generate_data(T,Nneurons,Dt,NeuronNoise,SeqNoiseTime,SeqNoiseNeuron,gap,0,0,neg,1);

figure; SimpleWHPlot(W,H,X); title('generated data','Fontsize',16)
set(gcf,'position',[200,200,1200,900])
nuc_norm = norm(svd(X),1);
X = X/nuc_norm*size(X,1);

%% Split into training and test set, find best params with BayesOpt
split = floor(T*.75);
Xtrain = X(:,1:split);
Xtest = X(:,(split+1):end);

K = 5;
L = 50;
neg = 0;
lambda_var = optimizableVariable('lambda', [1e-5, 1e-1], 'Transform', 'log');
alpha_var = optimizableVariable('alpha', [1e-6, 1e-2], 'Transform', 'log');
fun = @(x)compute_error_balance_score_2d(X,K,L,x.lambda,x.alpha);
results = bayesopt(fun, [lambda_var,alpha_var],'AcquisitionFunctionName','expected-improvement-plus','UseParallel',true, 'MaxObjectiveEvaluations',30, 'PlotFcn',[]);

lambda = results.bestPoint.lambda;
alpha = results.bestPoint.alpha;

%% Run FlexMF with best params 20 times, show number of significant factors
nIter = 20;
for iteri = 1:nIter
    [W_hat, H_hat, ~,loadings(iteri,:),power]= FlexMF(Xtrain,'K',K,'L',L,...
            'lambda', lambda, 'alpha', alpha, 'neg_prop', neg, 'maxiter', 50, 'showPlot', 0); 
    p = .05;
    [pvals(iteri,:),is_significant(iteri,:)] = test_significance(Xtest,W_hat,p);
    display(['FlexMF run ' num2str(iteri) '/' num2str(nIter)])
    display(['lambda=' num2str(lambda) ', alpha=' num2str(alpha)])
end

figure; hold on
h = histogram(sum(is_significant,2), 'edgecolor', 'w', 'facecolor', .7*[1 1 1]); 
h.BinCounts = h.BinCounts/sum(h.BinCounts)*100; 
xlim([0 K]); 
xlabel('# significant factors')
ylabel('% FlexMF runs')