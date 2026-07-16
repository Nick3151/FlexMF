%% Visualize SBI simulation examples
clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'Utils'))
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))

outdir = 'Simulation_SBI';
if ~exist(outdir, 'dir'); mkdir(outdir); end

%% Noise conditions
sel = 1;    % Selected noise level
noise_type = 'shape';
[f1,f2,f3] = simulation_plot_motifs(noise_type, sel, outdir);
close all

sel = 5;    % Selected noise level
noise_type = 'noise';
[f1,f2,f3] = simulation_plot_motifs(noise_type, sel, outdir);
close all

sel = 5;    % Selected noise level
noise_type = 'participate';
[f1,f2,f3] = simulation_plot_motifs(noise_type, sel, outdir);
close all

sel = 5;    % Selected noise level
noise_type = 'jitter';
[f1,f2,f3] = simulation_plot_motifs(noise_type, sel, outdir);
close all

sel = 5;    % Selected noise level
noise_type = 'warp';
[f1,f2,f3] = simulation_plot_motifs(noise_type, sel, outdir);
close all

sel = 5;    % Selected noise level
noise_type = 'overlap';
[f1,f2,f3] = simulation_plot_motifs(noise_type, sel, outdir);


function [f1,f2,f3] = simulation_plot_motifs(noise_type, sel, outdir)
load(sprintf('simulate_results_%s.mat', noise_type))
i = 1;  % Simulation seed

W = Ws{i,sel};
H = Hs{i,sel};
W_hat_SeqNMF = W_hats_SeqNMF{i,sel};
W_hat_SBI = W_hats_FlexMF{i,sel};
H_hat_SeqNMF = H_hats_SeqNMF{i,sel};
H_hat_SBI = H_hats_FlexMF{i,sel};
is_significant_SeqNMF = is_significants_SeqNMF{i,sel};
is_significant_FlexMF = is_significants_FlexMF{i,sel};
[N,K,L] = size(W);
TrainingData = TrainingDatas{i,sel};
[~,T] = size(TrainingData);
training_trials = T/L;

% pair with ground truth motifs
[coeff_W_SeqNMF, coeff_H_SeqNMF, ids_SeqNMF] = helper.similarity_WH(W, H, W_hat_SeqNMF, H_hat_SeqNMF);
[coeff_W_SBI, coeff_H_SBI, ids_SBI] = helper.similarity_WH(W, H, W_hat_SBI, H_hat_SBI);
ids_SeqNMF_match = ids_SeqNMF;
ids_SeqNMF_match(~ids_SeqNMF) = setdiff(1:K,ids_SeqNMF);
ids_SBI_match = ids_SBI;
ids_SBI_match(~ids_SBI) = setdiff(1:K,ids_SBI);
W_hat_SeqNMF_match = zeros(N,K,L);
H_hat_SeqNMF_match = zeros(K,T);
is_significant_SeqNMF_match = zeros(K,1);
W_hat_SeqNMF_match(:,ids_SeqNMF_match,:) = W_hat_SeqNMF;
H_hat_SeqNMF_match(ids_SeqNMF_match,:) = H_hat_SeqNMF;
for k=1:K
    if ids_SeqNMF(k)
        is_significant_SeqNMF_match(ids_SeqNMF(k),:) = is_significant_SeqNMF(k);
    end
end
W_hat_SBI_match = zeros(N,K,L);
H_hat_SBI_match = zeros(K,T);
is_significant_FlexMF_match = zeros(K,1);
W_hat_SBI_match(:,ids_SBI_match,:) = W_hat_SBI;
H_hat_SBI_match(ids_SBI_match,:) = H_hat_SBI;
for k=1:K
    if ids_SBI(k)
        is_significant_FlexMF_match(ids_SBI(k),:) = is_significant_FlexMF(k);
    end
end

% if strcmp(noise_type, 'jitter')||strcmp(noise_type, 'warp')
%     L = 100;
% else
%     L = 50;
% end

plotAll = 1;
f1 = figure;
SimpleWHPlot_patch(W,H,'trials', training_trials, 'frames', L, 'Data', TrainingData, 'plotAll', plotAll)
title(sprintf('Simulated data %s', noise_type))
set(f1,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
f2 = figure;
SimpleWHPlot_patch(W_hat_SeqNMF_match,H_hat_SeqNMF_match,'trials', training_trials, 'frames', L, 'is_significant', is_significant_SeqNMF_match, 'plotAll', plotAll)
title(sprintf('Simulated data %s SeqNMF', noise_type))
set(f2,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
f3 = figure;
SimpleWHPlot_patch(W_hat_SBI_match,H_hat_SBI_match,'trials', training_trials, 'frames', L, 'is_significant', is_significant_FlexMF_match, 'plotAll', plotAll)
title(sprintf('Simulated data %s SBI', noise_type))
set(f3,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

exportgraphics(f1, fullfile(outdir, sprintf('Simulated_data_%s.pdf', noise_type)), "ContentType", "vector")
exportgraphics(f2, fullfile(outdir, sprintf('Simulated_data_%s_SeqNMF.pdf', noise_type)), "ContentType", "vector")
exportgraphics(f3, fullfile(outdir, sprintf('Simulated_data_%s_SBI.pdf', noise_type)), "ContentType", "vector")
end
