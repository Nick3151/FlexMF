clear all
close all
clc
root = fileparts(pwd);
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))
addpath(fullfile(root, 'Utils'))

data_type = 'noise';
load(fullfile('Simulation_Results', sprintf('Compare_FlexMF_SeqNMF_%s.mat', data_type)))

noise_levels = 0:.005:.02;
participation_levels = 1:-0.1:.6; 
jitter_levels = 0:4;
warp_levels = 0:4;

switch data_type
    case 'warp'
        noise_levels_plot = warp_levels;
    case 'jitter'
        noise_levels_plot = jitter_levels;
    case 'participation'
        noise_levels_plot = participation_levels;
    case 'noise'
        noise_levels_plot = noise_levels;
    case 'warpnoise'
        noise_levels_plot = warp_levels;
    case 'jitternoise'
        noise_levels_plot = jitter_levels;
end      


mean_emds_W_SeqNMF = cellfun(@mean, emds_W_SeqNMF);
mean_emds_H_SeqNMF = cellfun(@mean, emds_H_SeqNMF);
mean_emds_W_FlexMF = cellfun(@mean, emds_W_FlexMF);
mean_emds_H_FlexMF = cellfun(@mean, emds_H_FlexMF);

[nNoise, nSim] = size(mean_emds_W_SeqNMF);

% Build matching x-positions with small side-by-side offsets
offset = 0.18;          % horizontal separation between methods
x_seq  = repmat(1:nNoise,1,nSim) - offset;
x_flex = repmat(1:nNoise,1,nSim) + offset;

% Colors (change if you prefer grayscale)
c_seq  = [0.25 0.55 0.95];
c_flex = [0.92 0.35 0.35];

% ---- Plot ----
figure; 
set(gcf, 'Position', [100, 100, 500, 1000])
ax1 = subplot('Position', [0.1 0.65 0.85 0.3]);
hold on; box off; 
% Swarm points (jitter only within each x-position)
s1 = swarmchart(x_seq,  mean_emds_W_SeqNMF(:),  20, c_seq,  'filled', ...
    'XJitter','density','XJitterWidth',0.2,'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.7);
s2 = swarmchart(x_flex, mean_emds_W_FlexMF(:), 20, c_flex, 'filled', ...
    'XJitter','density','XJitterWidth',0.2,'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.7);

% Axes and labels
xlim([1-0.6, nNoise+0.6]);
set(ax1, 'XTick', []);
ylabel('Mean EMD W');
% ylabel('Mean EMD H');

% ---- Optional: summary overlays ----
showMedian = false;   % set true/false
showMeanSEM = true; % set true/false

bw = 0.1;  % half-width for horizontal summary lines

for i = 1:nNoise
    % data at this level
    yi_seq  = mean_emds_W_SeqNMF(i, :);
    yi_flex = mean_emds_W_FlexMF(i, :);

    % centers
    xc_seq  = i - offset;
    xc_flex = i + offset;

    if showMedian
        mseq  = median(yi_seq,  'omitnan');
        mflex = median(yi_flex, 'omitnan');
        plot([xc_seq-bw,  xc_seq+bw],  [mseq  mseq],  'k-', 'LineWidth',1.3);
        plot([xc_flex-bw, xc_flex+bw], [mflex mflex], 'k-', 'LineWidth',1.3);
    end

    if showMeanSEM
        mu_seq  = mean(yi_seq,  'omitnan');
        mu_flex = mean(yi_flex, 'omitnan');
        se_seq  = std(yi_seq, 0, 'omitnan') / sqrt(sum(isfinite(yi_seq)));
        se_flex = std(yi_flex,0, 'omitnan') / sqrt(sum(isfinite(yi_flex)));

        % mean line
        plot([xc_seq-bw,  xc_seq+bw],  [mu_seq mu_seq],   '-', 'Color',[0.1 0.1 0.1], 'LineWidth',1.1);
        plot([xc_flex-bw, xc_flex+bw], [mu_flex mu_flex], '-', 'Color',[0.1 0.1 0.1], 'LineWidth',1.1);
        % SEM bars
        line([xc_seq xc_seq],   [mu_seq-se_seq,  mu_seq+se_seq],   'Color',[0.1 0.1 0.1], 'LineWidth',1.1);
        line([xc_flex xc_flex], [mu_flex-se_flex,mu_flex+se_flex], 'Color',[0.1 0.1 0.1], 'LineWidth',1.1);
    end
end
hold off

ax2 = subplot('Position', [0.1 0.3 0.85 0.3]);
hold on; box off; 
% Swarm points (jitter only within each x-position)
s1 = swarmchart(x_seq,  mean_emds_H_SeqNMF(:),  20, c_seq,  'filled', ...
    'XJitter','density','XJitterWidth',0.2,'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.7);
s2 = swarmchart(x_flex, mean_emds_H_FlexMF(:), 20, c_flex, 'filled', ...
    'XJitter','density','XJitterWidth',0.2,'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.7);

% Axes and labels
xlim([1-0.6, nNoise+0.6]);
set(ax2, 'XTick', []);
ylabel('Mean EMD H');

% ---- Optional: summary overlays ----
showMedian = false;   % set true/false
showMeanSEM = true; % set true/false

bw = 0.1;  % half-width for horizontal summary lines

for i = 1:nNoise
    % data at this level
    yi_seq  = mean_emds_H_SeqNMF(i, :);
    yi_flex = mean_emds_H_FlexMF(i, :);

    % centers
    xc_seq  = i - offset;
    xc_flex = i + offset;

    if showMedian
        mseq  = median(yi_seq,  'omitnan');
        mflex = median(yi_flex, 'omitnan');
        plot([xc_seq-bw,  xc_seq+bw],  [mseq  mseq],  'k-', 'LineWidth',1.3);
        plot([xc_flex-bw, xc_flex+bw], [mflex mflex], 'k-', 'LineWidth',1.3);
    end

    if showMeanSEM
        mu_seq  = mean(yi_seq,  'omitnan');
        mu_flex = mean(yi_flex, 'omitnan');
        se_seq  = std(yi_seq, 0, 'omitnan') / sqrt(sum(isfinite(yi_seq)));
        se_flex = std(yi_flex,0, 'omitnan') / sqrt(sum(isfinite(yi_flex)));

        % mean line
        plot([xc_seq-bw,  xc_seq+bw],  [mu_seq mu_seq],   '-', 'Color',[0.1 0.1 0.1], 'LineWidth',1.1);
        plot([xc_flex-bw, xc_flex+bw], [mu_flex mu_flex], '-', 'Color',[0.1 0.1 0.1], 'LineWidth',1.1);
        % SEM bars
        line([xc_seq xc_seq],   [mu_seq-se_seq,  mu_seq+se_seq],   'Color',[0.1 0.1 0.1], 'LineWidth',1.1);
        line([xc_flex xc_flex], [mu_flex-se_flex,mu_flex+se_flex], 'Color',[0.1 0.1 0.1], 'LineWidth',1.1);
    end
end

legend([s1 s2], {'SeqNMF','FlexMF'}, 'Location','northwest', 'Box','off');

hold off;

% Number of significant motifs
num_significant_seq = cellfun(@nnz, is_significant_SeqNMF);
num_significant_flex = cellfun(@nnz, is_significant_FlexMF);

ax3 = subplot('Position', [0.1 0.1 0.85 0.15]);
hold on
s3 = swarmchart(x_seq, num_significant_seq(:), 20, c_seq,  'filled', ...
    'XJitter','density','XJitterWidth',0.2,'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.7);
s4 = swarmchart(x_flex, num_significant_flex(:), 20, c_flex,  'filled', ...
    'XJitter','density','XJitterWidth',0.2,'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.7);
% Axes and labels
xlim([1-0.6, nNoise+0.6]);
set(ax3, 'XTick', 1:nNoise, 'XTickLabel', string(noise_levels_plot));
xlabel([data_type, ' Level']);
ylabel('Number Significant');

save2pdf(['Compare_FlexMF_SeqNMF_ ', data_type])

%% Look at factors
i = 4;
n = 1;
plotAll = 1;
figure; SimpleWHPlot(Whats_FlexMF{i,n}, Hhats_train_FlexMF{i,n}, 'plotAll', plotAll); title('FlexMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot(Whats_FlexMF{i,n}, Hhats_train_FlexMF{i,n}, 'Data', Xs_train{i,n}, 'plotAll', plotAll); title('FlexMF factors, with raw data')
% set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
save2pdf(sprintf('FlexMF_EMD_%s = %0.2f.pdf', data_type, noise_levels_plot(i)))

figure; SimpleWHPlot(Whats_SeqNMF{i,n}, Hhats_train_SeqNMF{i,n}, 'plotAll', plotAll); title('SeqNMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot(Whats_SeqNMF{i,n}, Hhats_train_SeqNMF{i,n}, 'Data', Xs_train{i,n}, 'plotAll', plotAll); title('SeqNMF factors, with raw data')
% set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
save2pdf(sprintf('SeqNMF_%s = %0.2f.pdf', data_type, noise_levels_plot(i)))

%%
[emds_W_SeqNMF{i,n}, emds_H_SeqNMF{i,n}, ids_SeqNMF{i,n}] = helper.similarity_WH_EMD(Ws{i,n}, Hs_train{i,n}, Whats_SeqNMF{i,n}, Hhats_train_SeqNMF{i,n});
[emds_W_FlexMF{i,n}, emds_H_FlexMF{i,n}, ids_FlexMF{i,n}] = helper.similarity_WH_EMD(Ws{i,n}, Hs_train{i,n}, Whats_FlexMF{i,n}, Hhats_train_FlexMF{i,n});