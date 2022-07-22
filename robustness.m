clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))

% Common settings
number_of_seqences = 3;
T = 3000; % length of data to generate
Nneurons = 10*ones(number_of_seqences,1); % number of neurons in each sequence
Dt = 3.*ones(number_of_seqences,1); % gap between each member of the sequence

%% Robustness to participation probability
N_trials = 1;
prob_dropout = 0:.1:.9;
scores = zeros(N_trials, length(prob_dropout));

for i = 1:length(prob_dropout)
    for n = 1:N_trials
        NeuronNoise = 0; % probability of added noise in each bin
        SeqNoiseTime = zeros(number_of_seqences,1); % Jitter parameter = 0%
        SeqNoiseNeuron = (1-prob_dropout(i)).*ones(number_of_seqences,1); % Participation parameter = 100%
        neg = 0.2; % Proportion of negative indices in W
        [X, W, H, V_hat] = generate_data(T,Nneurons,Dt,NeuronNoise,SeqNoiseTime,SeqNoiseNeuron,0,0,neg,1);
    %     figure; SimpleWHPlot(W,H,X); title('generated data','Fontsize',16)
    %     set(gcf,'position',[200,200,1200,900])
        nuc_norm = norm(svd(X),1);
        X = X/nuc_norm*size(X,1);

        % Run FlexMF
        shg; clf
        K = 5;
        L = 50;
        lambda = 1e-2;
        alpha = 1e-2;
        lambdaL1H = 0;
        lambdaL1W = .1;
        display(['Participation probability ' num2str(1-prob_dropout(i))])
        tic
        [W_hat,H_hat,cost,loadings,power] = FlexMF(X,'K',K, 'L', L, 'maxiter', 50,...
            'lambda', lambda, 'alpha', alpha, 'lambdaL1W', lambdaL1W, 'lambdaL1H', lambdaL1H);
        toc
        set(gcf,'position',[200,200,1200,900])
        scores(n,i) = helper.similarity(W, H, W_hat, H_hat);
    end
end
figure;
plot(prob_dropout, scores)
xlabel('% Dropout')
ylabel('Similarity to ground truth')