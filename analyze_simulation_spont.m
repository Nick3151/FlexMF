clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(fullfile(root, 'Utils'))
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))
addpath(genpath(fullfile(root, 'FlexMF')));

Trials = 200;
K = 5;
Length = 50;
Nmotifs = 16:-2:8;
Nneurons = 5*ones(K, 1); % the number of neurons in each motif
Dt = 3.*ones(K,1); % gap between each member of the motif
thresh_active_freq_all = 0:.05:1;
spontaneous_rate = .003;

[true_pos_rates, false_pos_rates] = simulation_spont(Trials, Length, Nmotifs, Nneurons, Dt, ...
    'thresh_active_freq_all', thresh_active_freq_all, 'spontaneous_rate', spontaneous_rate);

%% load and plot ROC
thresh_active_freq_all = 0:.05:1;
spontaneous_rate = .003;
load(sprintf('simulate_results_spont=%0.3f.mat', spontaneous_rate))
figure; 
colors_tp = [0 0 1;
    .15 .15 1;
    .3 .3 1;
    .45 .45 1;
    .6 .6 1];
colors_fp = [1 0 0;
    1 .15 .15;
    1 .3 .3;
    1 .45 .45;
    1 .6 .6];
colororder([colors_tp;colors_fp])
plot(thresh_active_freq_all, true_pos_rates, 'LineWidth',2)
hold on
plot(thresh_active_freq_all, false_pos_rates, 'LineWidth',2)
xlim([0,1]), ylim([0,1])
legend({'TPR, N=16', 'TPR, N=14', 'TPR, N=12', 'TPR, N=10', 'TPR, N=8', ...
    'FPR, N=16', 'FPR, N=14', 'FPR, N=12', 'FPR, N=10', 'FPR, N=8'})
xlabel('Neuron Active thresh')
save2pdf(sprintf('Determine neurons active thresh to exclude spontaneous activities spont=%0.3f', spontaneous_rate))

figure;
hold on
for i=1:5
    plot(false_pos_rates(i,:), true_pos_rates(i,:), 'LineWidth',2)
end
xlim([0,.2])
ylim([0,1])
legend({'N=16', 'N=14', 'N=12', 'N=10', 'N=8'})