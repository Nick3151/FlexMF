%% Compare results: prob of being significant vs #occur 
clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(fullfile(root, 'Utils'))
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))
addpath(genpath(fullfile(root, 'FlexMF')));

outdir = 'Simulation_SBI';
if ~exist(outdir, 'dir'); mkdir(outdir); end

color_palet_blue = [[.8 .8 1]; [.5 .5 1];[.2 .2 1]; [0 0 .8]; [0 0 .5]];
color_palet_red = flip(color_palet_blue,2);

%% Shape conditions
sel = 1;
noise_levels = {'Calcium transients', 'Bursts', 'Single spikes'};
noise_type = 'Motif shape';
[f1, f2, f3] = simulation_compare_algorithm_plot('simulate_results_shape.mat', noise_levels, noise_type, sel, 'lambda', .005);
save2pdf(fullfile(outdir, 'Simulation_results_shape'), f1);
save2pdf(fullfile(outdir, 'Simulation_results_sparsity_shape'), f2);
save2pdf(fullfile(outdir, 'Simulation_results_error_shape'), f3);
% save2pdf(fullfile(outdir, 'Simulation_results_min_motif_occur_shape'), f4);

%% Noise conditions
sel = 5;
noise_levels = 100*(0:0.01:0.1); % probability of added noise in each bin
noise_type = '% Additive Noise';
[f1, f2, f3] = simulation_compare_algorithm_plot("simulate_results_noise.mat", noise_levels, noise_type, sel, 'lambda', .003);
save2pdf(fullfile(outdir, 'Simulation_results_noise'), f1);
save2pdf(fullfile(outdir, 'Simulation_results_sparsity_noise'), f2);
save2pdf(fullfile(outdir, 'Simulation_results_error_noise'), f3);
% save2pdf(fullfile(outdir, 'Simulation_results_min_motif_occur_noise'), f4);

%% Participation conditions
sel = 5;
participation_rates = 1:-0.1:.1; % Participation probability
dropout_rates = 100*(1-participation_rates);
noise_type = '% Dropout';
[f1, f2, f3] = simulation_compare_algorithm_plot("simulate_results_participate.mat", dropout_rates, noise_type, sel);
% [f1, f2, f3] = simulation_compare_algorithm_plot("simulate_results_participate.mat", dropout_rates, noise_type, sel, 'metric', 'F1');

save2pdf(fullfile(outdir, 'Simulation_results_participate'), f1);
save2pdf(fullfile(outdir, 'Simulation_results_sparsity_participate'), f2);
save2pdf(fullfile(outdir, 'Simulation_results_error_participate'), f3);
% save2pdf(fullfile(outdir, 'Simulation_results_min_motif_occur_participate'), f4);

%% Jittering conditions
sel = 5;
jitter_stds = 0:9; % Jitter time std
noise_type = 'Jittering Std';
% [f1, f2, f3] = simulation_compare_algorithm_plot("simulate_results_jitter.mat", jitter_stds, noise_type, sel);
[f1, f2, f3] = simulation_compare_algorithm_plot("simulate_results_jitter.mat", jitter_stds, noise_type, sel, 'metric', 'F1');

% save2pdf(fullfile(outdir, 'Simulation_results_jittering'), f1);
save2pdf(fullfile(outdir, 'Simulation_F1_results_jittering'), f1);
save2pdf(fullfile(outdir, 'Simulation_results_sparsity_jittering'), f2);
save2pdf(fullfile(outdir, 'Simulation_results_error_jittering'), f3);
% save2pdf(fullfile(outdir, 'Simulation_results_min_motif_occur_jittering'), f4);

%% Warping conditions
sel = 5;
warp_levels = 0:9; % the maximum warping time
noise_type = 'Warping Level';
% [f1, f2, f3] = simulation_compare_algorithm_plot("simulate_results_warp.mat", warp_levels, noise_type, sel);
[f1, f2, f3] = simulation_compare_algorithm_plot("simulate_results_warp.mat", warp_levels, noise_type, sel, 'metric', 'F1');

% save2pdf(fullfile(outdir, 'Simulation_results_warping'), f1);
save2pdf(fullfile(outdir, 'Simulation_F1_results_warping'), f1);
save2pdf(fullfile(outdir, 'Simulation_results_sparsity_warping'), f2);
save2pdf(fullfile(outdir, 'Simulation_results_error_warping'), f3);
% save2pdf(fullfile(outdir, 'Simulation_results_min_motif_occur_warping'), f4);

%% Overlapping conditions
sel = 5;
overlaps_n = 0:0.1:0.9; % Neuron overlap rate
noise_type = 'Neurons Overlap Ratio';
[f1, f2, f3] = simulation_compare_algorithm_plot("simulate_results_overlap.mat", overlaps_n, noise_type, sel);

save2pdf(fullfile(outdir, 'Simulation_results_overlap'), f1);
save2pdf(fullfile(outdir, 'Simulation_results_sparsity_overlap'), f2);
save2pdf(fullfile(outdir, 'Simulation_results_error_overlap'), f3);
% save2pdf(fullfile(outdir, 'Simulation_results_min_motif_occur_overlap'), f4);
