% Significance test of ground truth factor vs number of occurrences in the
% data
clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(fullfile(root, 'Utils'))
addpath(fullfile(root, 'MATLAB-tools'))

%% Generate some synthetic data, example
Trials = 30; % total number of trials
L = 50; % length of each trial
K = 3;
Nmotifs = 6*ones(K, 1); % the number of occurences of motifs
Nneurons = 5*ones(K, 1); % the number of neurons in each motif
Dt = 3.*ones(K,1); % gap between each member of the motif
len_burst = 1; % Continuous firing time
dynamic = 0; % Consider calcium dynamics or not
seed = 1;
[X, W, H, X_hat, ~] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, 'len_burst', len_burst, 'dynamic', dynamic, 'seed', seed);

N = size(W,1);
TestData = zeros(N,Trials*L);
for t=1:Trials
    TestData(:,(t-1)*L+1:t*L) = squeeze(X(:,:,t));
end

figure; SimpleWHPlot_trials(W, H, [], X, 1); title('generated data','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
% [pvals,is_significant] = test_significance_trials(TestData, Trials, L, W, 1);

%% In how much percent of trials should the motif present in order to be significant
Trials = 100;
% Trials = 10;
L = 50; % length of each trial

Nmotifs = 1:10;
K = 10;
Nneurons = 5*ones(K, 1); % the number of neurons in each motif
Magnitudes = ones(K, 1); % the activation magnitudes of each motif
Dt = 3.*ones(K,1); % gap between each member of the motif
noise = 0.01; % probability of added noise in each bin
jitter = 0*ones(K,1); % Jitter time std
participation = .7.*ones(K,1); % Participation probability = 100%
warp = 0; % the maximum warping time
neg = 0; % Proportion of negative indices in W
nsim = 100;
seeds = randperm(1000, nsim);
pvals_calc = zeros(nsim,K);
is_significant_calc = zeros(nsim,K);
pvals_burst = zeros(nsim,K);
is_significant_burst = zeros(nsim,K);
pvals_spike = zeros(nsim,K);
is_significant_spike = zeros(nsim,K);
parfor i=1:nsim
    len_burst = 20; % Continuous firing time
    dynamic = 1; % Consider calcium dynamics or not
    [X, W, H, X_hat, ~] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, ...
        'noise', noise, 'participation', participation, 'len_burst', len_burst, 'dynamic', dynamic, 'seed', seeds(i));
    N = size(W,1);
    TestData = zeros(N,Trials*L);
    for t=1:Trials
        TestData(:,(t-1)*L+1:t*L) = squeeze(X(:,:,t));
    end
    [pvals_calc(i,:),is_significant_calc(i,:)] = test_significance_trials(TestData, Trials, L, W);
    % [pvals,is_significant] = test_significance(TestData, W);

    len_burst = 20; % Continuous firing time
    dynamic = 0; % Consider calcium dynamics or not
    [X, W, H, X_hat, ~] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, ...
        'noise', noise, 'participation', participation, 'len_burst', len_burst, 'dynamic', dynamic, 'seed', seeds(i));
    N = size(W,1);
    TestData = zeros(N,Trials*L);
    for t=1:Trials
        TestData(:,(t-1)*L+1:t*L) = squeeze(X(:,:,t));
    end
    [pvals_burst(i,:),is_significant_burst(i,:)] = test_significance_trials(TestData, Trials, L, W);

    len_burst = 1; % Continuous firing time
    dynamic = 0; % Consider calcium dynamics or not
    [X, W, H, X_hat, ~] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, ...
        'noise', noise, 'participation', participation, 'len_burst', len_burst, 'dynamic', dynamic, 'seed', seeds(i));
    N = size(W,1);
    TestData = zeros(N,Trials*L);
    for t=1:Trials
        TestData(:,(t-1)*L+1:t*L) = squeeze(X(:,:,t));
    end
    [pvals_spike(i,:),is_significant_spike(i,:)] = test_significance_trials(TestData, Trials, L, W);
end

colors = distinguishable_colors(10,{'w','k'});
figure; hold on
swarmchart(Nmotifs, pvals_spike', 24, colors(1,:), 'filled')
p1 = plot(Nmotifs, mean(pvals_spike), 'linewidth', 2, 'Color', colors(1,:));
swarmchart(Nmotifs, pvals_burst', 24, colors(2,:), 'filled')
p2 = plot(Nmotifs, mean(pvals_burst), 'linewidth', 2, 'Color', colors(2,:));
swarmchart(Nmotifs, pvals_calc', 24, colors(3,:), 'filled')
p3 = plot(Nmotifs, mean(pvals_calc), 'linewidth', 2, 'Color', colors(3,:));
xlabel('# Motif Occurences')
ylabel('p values')
xlim([0, max(Nmotifs)]); ylim([0,0.2])
legend([p1,p2,p3], {'Single spike', 'Sustained spikes', 'Calcium Transient'})

figure; hold on
plot(Nmotifs, mean(is_significant_spike), 'linewidth', 2, 'Color', colors(1,:))  
plot(Nmotifs, mean(is_significant_burst), 'linewidth', 2, 'Color', colors(2,:))  
plot(Nmotifs, mean(is_significant_calc), 'linewidth', 2, 'Color', colors(3,:))    
xlabel('# Motif Occurences')
ylabel('Probability of being significant')
xlim([0, max(Nmotifs)]); ylim([0,1])
legend({'Single spike', 'Burst', 'Calcium Transient'}, 'Location','southeast')

%% The effect of magnitudes
Trials = 30;
L = 50; % length of each trial
K = 10;
Nmotifs = 6*ones(1,K);
Nneurons = 5*ones(K, 1); % the number of neurons in each motif
Magnitudes = 1:K; % the activation magnitudes of each motif
Dt = 3.*ones(K,1); % gap between each member of the motif
noise = 0.01; % probability of added noise in each bin
participation = .7.*ones(K,1); % Participation probability = 100%
overlap_t = 1;
len_burst = 20; % Continuous firing time
dynamic = 1; % Consider calcium dynamics or not
nsim = 100;
seeds = randperm(1000, nsim);
pvals_calc = zeros(nsim,K);
is_significant_calc = zeros(nsim,K);

parfor i=1:nsim
    [X, W, H, X_hat, ~] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, 'Magnitudes', Magnitudes,...
        'noise', noise, 'participation', participation, 'overlap_t', overlap_t, 'len_burst', len_burst, 'dynamic', dynamic, 'seed', seeds(i));
    
    N = size(W,1);
    TestData = zeros(N,Trials*L);
    for t=1:Trials
        TestData(:,(t-1)*L+1:t*L) = squeeze(X(:,:,t));
    end
    
    [pvals_calc(i,:),is_significant_calc(i,:)] = test_significance_trials(TestData, Trials, L, W);
end

figure; SimpleWHPlot_trials(W, H, [], X, 1); title('generated data','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

figure; hold on
swarmchart(Magnitudes, pvals_calc', 24, colors(1,:), 'filled')
plot(Magnitudes, mean(pvals_calc), 'linewidth', 2, 'Color', colors(1,:));
xlabel('Magnitude of motif')
ylabel('p values')
xlim([0, max(Magnitudes)]); ylim([0,0.2])

figure; hold on
plot(Magnitudes, mean(is_significant_calc), 'linewidth', 2, 'Color', colors(1,:)) 
xlabel('Magnitude of motif')
ylabel('Probability of being significant')
xlim([0, max(Magnitudes)]); ylim([0,1])

%% Impact of noise
Trials = 30;
% Trials = 10;
L = 50; % length of each trial

Nmotifs = 1:10;
K = 10;
Nneurons = 5*ones(K, 1); % the number of neurons in each motif
Dt = 3.*ones(K,1); % gap between each member of the motif
noise_levels = [0.001, 0.01, 0.1]; % probability of added noise in each bin
overlap_t = 1;
len_burst = 20; % Continuous firing time
dynamic = 1; % Consider calcium dynamics or not
nsim = 100;
seeds = randperm(1000, nsim);
pvals = zeros(nsim,K);
is_significant = zeros(nsim,K);
f1 = figure; hold on
ax1 = gca;
f2 = figure; hold on
ax2 = gca;
for j=1:length(noise_levels)
    noise = noise_levels(j);
    parfor i=1:nsim
        [X, W, H, X_hat, ~] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt,...
        'noise', noise, 'overlap_t', overlap_t, 'len_burst', len_burst, 'dynamic', dynamic, 'seed', seeds(i));
        N = size(W,1);
        TestData = zeros(N,Trials*L);
        for t=1:Trials
            TestData(:,(t-1)*L+1:t*L) = squeeze(X(:,:,t));
        end
        [pvals(i,:),is_significant(i,:)] = test_significance_trials(TestData, Trials, L, W);
    end

    colors = distinguishable_colors(10,{'w','k'});
    swarmchart(ax1, Nmotifs, pvals', 24, colors(j,:), 'filled')
    p(j) = plot(ax1, Nmotifs, mean(pvals), 'linewidth', 2, 'Color', colors(j,:));
    plot(ax2, Nmotifs, mean(is_significant), 'linewidth', 2, 'Color', colors(j,:))
end

xlabel(ax1, '# Motif Occurences')
ylabel(ax1, 'p values')
xlim(ax1, [0, max(Nmotifs)]); ylim(ax1, [0,0.2])
legend(ax1, p, {'noise=.001', 'noise=.01', 'noise=.1'})
 
xlabel(ax2, '# Motif Occurences')
ylabel(ax2, 'Probability of being significant')
xlim(ax2, [0, max(Nmotifs)]); ylim([0,1])
legend(ax2, {'noise=.001', 'noise=.01', 'noise=.1'})

%% Impact of spike length
Trials = 30;
% Trials = 10;
L = 50; % length of each trial

Nmotifs = 1:10;
K = 10;
Nneurons = 5*ones(K, 1); % the number of neurons in each motif
Dt = 3.*ones(K,1); % gap between each member of the motif
overlap_t = 1;
lens_burst = [1,10,20,30]; % Continuous firing time
dynamic = 1; % Consider calcium dynamics or not
nsim = 100;
seeds = randperm(1000, nsim);
pvals = zeros(nsim,K);
is_significant = zeros(nsim,K);
f1 = figure; hold on
ax1 = gca;
f2 = figure; hold on
ax2 = gca;
for j=1:length(lens_burst)
    len_burst = lens_burst(j);
    parfor i=1:nsim
        [X, W, H, X_hat, ~] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt,...
        'overlap_t', overlap_t, 'len_burst', len_burst, 'dynamic', dynamic, 'seed', seeds(i));
        N = size(W,1);
        TestData = zeros(N,Trials*L);
        for t=1:Trials
            TestData(:,(t-1)*L+1:t*L) = squeeze(X(:,:,t));
        end
        [pvals(i,:),is_significant(i,:)] = test_significance_trials(TestData, Trials, L, W);
    end

    colors = distinguishable_colors(10,{'w','k'});
    swarmchart(ax1, Nmotifs, pvals', 24, colors(j,:), 'filled')
    p(j) = plot(ax1, Nmotifs, mean(pvals), 'linewidth', 2, 'Color', colors(j,:));
    plot(ax2, Nmotifs, mean(is_significant), 'linewidth', 2, 'Color', colors(j,:))
end

xlabel(ax1, '# Motif Occurences')
ylabel(ax1, 'p values')
xlim(ax1, [0, max(Nmotifs)]); ylim(ax1, [0,0.2])
legend(ax1, p, {'spike length = 1','spike length = 10','spike length = 20','spike length = 30'})
 
xlabel(ax2, '# Motif Occurences')
ylabel(ax2, 'Probability of being significant')
xlim(ax2, [0, max(Nmotifs)]); ylim([0,1])
legend(ax2, {'spike length = 1','spike length = 10','spike length = 20','spike length = 30'})

%% Impact of participation rate
Trials = 30;
% Trials = 10;
L = 50; % length of each trial

Nmotifs = 1:10;
K = 10;
Nneurons = 5*ones(K, 1); % the number of neurons in each motif
Dt = 3.*ones(K,1); % gap between each member of the motif
noise = 0.001; % probability of added noise in each bin
participation_rates = 0.5:0.1:0.9; % Participation probability
overlap_t = 1;
len_burst = 10; % Continuous firing time
dynamic = 1; % Consider calcium dynamics or not
nsim = 100;
seeds = randperm(1000, nsim);
pvals = zeros(nsim,K);
is_significant = zeros(nsim,K);
f1 = figure; hold on
ax1 = gca;
f2 = figure; hold on
ax2 = gca;
for j=1:length(participation_rates)
    participation = participation_rates(j).*ones(K,1);
    parfor i=1:nsim
        [X, W, H, X_hat, ~] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt,...
        'participation', participation, 'overlap_t', overlap_t, 'len_burst', len_burst, 'dynamic', dynamic, 'seed', seeds(i));
        N = size(W,1);
        TestData = zeros(N,Trials*L);
        for t=1:Trials
            TestData(:,(t-1)*L+1:t*L) = squeeze(X(:,:,t));
        end
        [pvals(i,:),is_significant(i,:)] = test_significance_trials(TestData, Trials, L, W);
    end

    colors = distinguishable_colors(10,{'w','k'});
    swarmchart(ax1, Nmotifs, pvals', 24, colors(j,:), 'filled')
    p(j) = plot(ax1, Nmotifs, mean(pvals), 'linewidth', 2, 'Color', colors(j,:));
    plot(ax2, Nmotifs, mean(is_significant), 'linewidth', 2, 'Color', colors(j,:))
end

xlabel(ax1, '# Motif Occurences')
ylabel(ax1, 'p values')
xlim(ax1, [0, max(Nmotifs)]); ylim(ax1, [0,0.2])
legend(ax1, p, {'participation=0.5','participation=0.6','participation=0.7','participation=0.8','participation=0.9'})
 
xlabel(ax2, '# Motif Occurences')
ylabel(ax2, 'Probability of being significant')
xlim(ax2, [0, max(Nmotifs)]); ylim([0,1])
legend(ax2, {'participation=0.5','participation=0.6','participation=0.7','participation=0.8','participation=0.9'})

%% Impact of temporal jittering
Trials = 30;
% Trials = 10;
L = 50; % length of each trial

Nmotifs = 1:10;
K = 10;
Nneurons = 5*ones(K, 1); % the number of neurons in each motif
Dt = 3.*ones(K,1); % gap between each member of
jitter_stds = 1:3; % Jitter time std
warp = 0; % the maximum warping time
overlap_t = 1;
len_burst = 10; % Continuous firing time
dynamic = 1; % Consider calcium dynamics or not
nsim = 100;
seeds = randperm(1000, nsim);
pvals = zeros(nsim,K);
is_significant = zeros(nsim,K);
f1 = figure; hold on
ax1 = gca;
f2 = figure; hold on
ax2 = gca;
for j=1:length(jitter_stds)
    jitter = jitter_stds(j).*ones(K,1);
    parfor i=1:nsim
        [X, W, H, X_hat, ~] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt,...
        'jitter', jitter, 'overlap_t', overlap_t, 'len_burst', len_burst, 'dynamic', dynamic, 'seed', seeds(i));
        N = size(W,1);
        TestData = zeros(N,Trials*L);
        for t=1:Trials
            TestData(:,(t-1)*L+1:t*L) = squeeze(X(:,:,t));
        end
        [pvals(i,:),is_significant(i,:)] = test_significance_trials(TestData, Trials, L, W);
    end

    colors = distinguishable_colors(10,{'w','k'});
    swarmchart(ax1, Nmotifs, pvals', 24, colors(j,:), 'filled')
    p(j) = plot(ax1, Nmotifs, mean(pvals), 'linewidth', 2, 'Color', colors(j,:));
    plot(ax2, Nmotifs, mean(is_significant), 'linewidth', 2, 'Color', colors(j,:))
end

xlabel(ax1, '# Motif Occurences')
ylabel(ax1, 'p values')
xlim(ax1, [0, max(Nmotifs)]); ylim(ax1, [0,0.2])
legend(ax1, p, {'jitter std = 1', 'jitter std = 2', 'jitter std = 3'})
 
xlabel(ax2, '# Motif Occurences')
ylabel(ax2, 'Probability of being significant')
xlim(ax2, [0, max(Nmotifs)]); ylim([0,1])
legend(ax2, {'jitter std = 1', 'jitter std = 2', 'jitter std = 3'})

%% Impact of temporal warping
Trials = 30;
% Trials = 10;
L = 50; % length of each trial

Nmotifs = 1:10;
K = 10;
Nneurons = 5*ones(K, 1); % the number of neurons in each motif
Dt = 3.*ones(K,1); % gap between each member of the motif
warp_levels = [1,2,3]; % the maximum warping time
overlap_t = 1;
len_burst = 10; % Continuous firing time
dynamic = 1; % Consider calcium dynamics or not
nsim = 100;
seeds = randperm(1000, nsim);
pvals = zeros(nsim,K);
is_significant = zeros(nsim,K);
f1 = figure; hold on
ax1 = gca;
f2 = figure; hold on
ax2 = gca;
for j=1:length(warp_levels)
    warp = warp_levels(j);
    parfor i=1:nsim
        [X, W, H, X_hat, ~] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt,...
        'warp', warp, 'overlap_t', overlap_t, 'len_burst', len_burst, 'dynamic', dynamic, 'seed', seeds(i));
        N = size(W,1);
        TestData = zeros(N,Trials*L);
        for t=1:Trials
            TestData(:,(t-1)*L+1:t*L) = squeeze(X(:,:,t));
        end
        [pvals(i,:),is_significant(i,:)] = test_significance_trials(TestData, Trials, L, W);
    end

    colors = distinguishable_colors(10,{'w','k'});
    swarmchart(ax1, Nmotifs, pvals', 24, colors(j,:), 'filled')
    p(j) = plot(ax1, Nmotifs, mean(pvals), 'linewidth', 2, 'Color', colors(j,:));
    plot(ax2, Nmotifs, mean(is_significant), 'linewidth', 2, 'Color', colors(j,:))
end

xlabel(ax1, '# Motif Occurences')
ylabel(ax1, 'p values')
xlim(ax1, [0, max(Nmotifs)]); ylim(ax1, [0,0.2])
legend(ax1, p, {'warp = 1', 'warp = 2', 'warp = 3'})
 
xlabel(ax2, '# Motif Occurences')
ylabel(ax2, 'Probability of being significant')
xlim(ax2, [0, max(Nmotifs)]); ylim([0,1])
legend(ax2, {'warp = 1', 'warp = 2', 'warp = 3'})
