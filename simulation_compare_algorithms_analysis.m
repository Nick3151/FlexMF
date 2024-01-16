%% Compare results: prob of being significant vs #occur 
clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(fullfile(root, 'Utils'))
addpath(fullfile(root, 'MATLAB-tools'))
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))
rmpath(genpath(fullfile(root, 'seqNMF-master')));
addpath(genpath(fullfile(root, 'FlexMF')));

color_palet_blue = [[.8 .8 1]; [.5 .5 1];[.2 .2 1]; [0 0 .8]; [0 0 .5]];
color_palet_red = flip(color_palet_blue,2);

%% Shape conditions
sel = 1;
[f1, f2, f3] = simulation_compare_algorithm_plot('simulate_results_shape.mat', 3, {'Calcium transients', 'Bursts', 'Single spikes'}, sel, 'lambda', .005);
save2pdf('Simulation_results_shape', f1);

%% Noise conditions
sel = 2;
[f1, f2, f3] = simulation_compare_algorithm_plot("simulate_results_noise.mat", 3, {'Noise=.001', 'Noise=.01', 'Noise=.1'}, sel, 'lambda', .003)
save2pdf('Simulation_results_noise', f1);

%% Participation conditions
sel = 3;
[f1, f2, f3] = simulation_compare_algorithm_plot("simulate_results_participate.mat", 5, ...
    {'Participate=.5', 'Participate=.6', 'Participate=.7', 'Participate=.8', 'Participate=.9'}, sel);
save2pdf('Simulation_results_participate', f1);
save2pdf('Simulation_results_sparsity', f2);
save2pdf('Simulation_results_error', f3);

%% Jittering conditions
sel = 2;
simulation_compare_algorithm_plot("simulate_results_jitter.mat", 3, ...
    {'Jitter Std=1', 'Jitter Std=2', 'Jitter Std=3'}, sel)

%% Warping conditions
sel = 2;
simulation_compare_algorithm_plot("simulate_results_warp.mat", 3, ...
    {'Warp=1', 'Warp=2', 'Warp=3'}, sel)

%% Overlapping conditions
sel = 3;
[f1, f2, f3] = simulation_compare_algorithm_plot("simulate_results_overlap.mat", 5, ...
    {'Overlap=0%', 'Overlap=20%', 'Overlap=40%', 'Overlap=60%', 'Overlap=80%'}, sel);
save2pdf('Simulation_results_overlap', f1);
save2pdf('Simulation_results_sparsity', f2);
save2pdf('Simulation_results_error', f3);
