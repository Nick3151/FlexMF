function [true_pos_rates, false_pos_rates] = simulation_spont(Trials, Length, Nmotifs, Nneurons, Dt, varargin)
% Simulate motifs with spontaneous activities
% Get motif neurons true positive and false positive rates under different
% activation threshold

%% Parse inputs
K = length(Nmotifs); % The number of motifs
if size(Nmotifs,1) ~= 1
    Nmotifs = Nmotifs';
end
validArrayK = @(x) isnumeric(x) && length(x)==K;

p  = inputParser;
addOptional(p, 'noise', .001, @isscalar)
addOptional(p, 'participation', .7.*ones(K, 1), validArrayK)
addOptional(p, 'additional_neurons', 50, @isscalar)
addOptional(p, 'spontaneous_rate', .001, @isscalar)
addOptional(p, 'nSim', 100, @isscalar)
addOptional(p, 'thresh_active_freq_all', 0:.05:1)
parse(p, varargin{:})

noise = p.Results.noise;
participation = p.Results.participation;
additional_neurons = p.Results.additional_neurons;
spontaneous_rate = p.Results.spontaneous_rate;
nSim = p.Results.nSim;
thresh_active_freq_all = p.Results.thresh_active_freq_all;

%%
warp = 0; % the maximum warping time
jitter = 0*ones(K,1);
len_burst = 1; % Continuous firing time
dynamic = 1; % Consider calcium dynamics or not
overlap_t = 0;
overlap_n = 0;
seeds = randperm(1000, nSim);

n_thresh = length(thresh_active_freq_all);
true_pos_rates = zeros(K,n_thresh,nSim);
false_pos_rates = zeros(K,n_thresh,nSim);
Ws = cell(nSim,1);
Whats = cell(nSim,1);
Hhat_trains = cell(nSim,1);
Hhat_tests = cell(nSim,1);
is_significants = cell(nSim,1);
X_trains = cell(nSim,1);
X_tests = cell(nSim,1);

parfor n=1:nSim
    [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, Length, Nmotifs, Nneurons, Dt, ...
        'participation', participation, 'additional_neurons', additional_neurons, 'spontaneous_rate', spontaneous_rate, 'seed', seeds(n));
    groups = zeros(Trials,1);
    for k=1:K
        groups(motif_ind{k}) = k;
    end
    
    % Dimension N*L*Trials
    rng(seeds(n))
    cv = cvpartition(groups, "Holdout",.5);
    ind_train = find(training(cv,1));
    X_train = X(:,:,ind_train);
    ind_test = find(test(cv,1));
    X_test = X(:,:,ind_test);
    
    % Dimension N*T
    N = size(W,1);
    TrainingData = zeros(N,cv.TrainSize*Length);
    for t=1:cv.TrainSize
        TrainingData(:,(t-1)*Length+1:t*Length) = squeeze(X_train(:,:,t));
    end
    TestData = zeros(N,cv.TestSize*Length);
    for t=1:cv.TestSize
        TestData(:,(t-1)*Length+1:t*Length) = squeeze(X_test(:,:,t));
    end
    [N,T] = size(TrainingData);

% f1 = figure;
% SimpleXplot_patch([TrainingData, TestData], [cv.TrainSize, cv.TestSize], Length); 
% set(f1,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
% 
% f2 = figure;
% H_train = H(:,ind_train);
% H_test = H(:,ind_test);
% SimpleWHPlot_trials(W, H_train, [], X_train, 1); title('generated data','Fontsize',16)
% set(f2,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
% 
% save2pdf('Simulated_data.pdf', f2)

    % Normalize training data
    frob_norm = norm(TrainingData(:));
    TrainingData = TrainingData/frob_norm*K;

%% Run FlexMF
    lambda = 0.01;
    alpha_W = 1e-6;
    alpha_H = 1e-3;
    
    lambdaL1H = 0;
    lambdaL1W = 0;
    
    display('Running FlexMF on 2p data')
%     figure;
%     set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
    tic
    [What,Hhat_train,~,errors_FlexMF,loadings,power] = FlexMF(TrainingData,'K',K, 'L', Length, 'maxiter', 50, 'tolerance', 1e-2,...
        'lambda', lambda, 'alpha_W', alpha_W, 'alpha_H', alpha_H, 'lambdaL1W', lambdaL1W, 'lambdaL1H', lambdaL1H, 'neg_prop', 0, 'showPlot', 0, 'verbal', 0);
    toc
    
    p = .05; % desired p value for factors
    display('Testing significance of factors on held-out data')
    [pvals,is_significant,is_single] = test_significance_new(TestData, What, [], p);

% %% Look at factors
% plotAll = 1;
% figure; SimpleWHPlot_patch(What, Hhat_train, 'trials', cv.TrainSize, 'frames', Length, 'is_significant', is_significant, 'plotAll', plotAll); title('FlexMF reconstruction')
% set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
% figure; SimpleWHPlot_patch(What, Hhat_train, 'trials', cv.TrainSize, 'frames', Length, 'is_significant', is_significant, 'Data', TrainingData, 'plotAll', plotAll); title('SeqNMF factors, with raw data')
% set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

%% refit W to test data
    frob_norm = norm(TestData(:));
    TestData = TestData/frob_norm*K;
    
    display('Refit W to test data')
    
    % figure; set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
    [What,Hhat_test,~,~,loadings_test,power_test] = FlexMF(TestData,'K',K, 'L', Length, 'maxiter', 1,...
        'lambda', lambda, 'alpha_W', alpha_W, 'alpha_H', alpha_H, ...
        'W_fixed', 1, 'W_init', What, 'SortFactors', 0, 'neg_prop', 0, 'showPlot', 0, 'verbal', 0);
    
    % figure; SimpleWHPlot_patch(What, Hhat_test, 'trials', cv.TestSize, 'frames', Length, 'is_significant', is_significant, 'plotAll', plotAll)
    % title('FlexMF reconstruction')
    % set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
    % figure; SimpleWHPlot_patch(What, Hhat_test, 'trials', cv.TestSize, 'frames', Length, 'Data', TestData, 'is_significant', is_significant, 'plotAll', plotAll)
    % title('FlexMF factors, with raw data')
    % set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

%% Compute useful things    
    neurons_active_true = cell(K,1);
    num_overlap = floor(overlap_n*Nneurons);
    j = 1;
    for k = 1:K
        neurons_active_true{k} = j:j+Nneurons(k)-1;
        j = j+Nneurons(k)-num_overlap(k);
    end
    neurons_spont = j:N;
%% Plot ROC for different neurons active frequency threshold
    for i=1:n_thresh
        thresh_active_freq = thresh_active_freq_all(i);
        [active_time_motifs_all, magnitudes_motifs_all, neurons_active_motifs, onsets_motifs] = ...
        analyze_motifs_activation([TrainingData, TestData], [Hhat_train,Hhat_test], What, 'thresh_active_freq', thresh_active_freq);
        [true_pos_rates_tmp, false_pos_rates_tmp, ids] = ...
            helper.match_motif_neurons_overlap(neurons_active_true, neurons_active_motifs, neurons_spont);
        for k=1:K
            id_tmp = find(ids==k);
            if ~isempty(id_tmp)
                true_pos_rates(k,i,n) = true_pos_rates_tmp(id_tmp);
                false_pos_rates(k,i,n) = false_pos_rates_tmp(id_tmp);
            end
        end
    end
    Ws{n} = W;
    Whats{n} = What;
    is_significants{n} = is_significant;
    Hhat_trains{n} = Hhat_train;
    Hhat_tests{n} = Hhat_test;
    
end

true_pos_rates = mean(true_pos_rates,3);
false_pos_rates = mean(false_pos_rates,3);

save('simulate_results_spont.mat', 'Ws', 'Whats', 'is_significants', ...
    'Hhat_trains', 'Hhat_tests', 'X_trains', 'X_tests', 'true_pos_rates', 'false_pos_rates')