clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(fullfile(root, 'Utils'))
addpath(fullfile(root, 'MATLAB-tools'))
rmpath(genpath(fullfile(root, 'seqNMF-master')));
addpath(genpath(fullfile(root, 'FlexMF')));

%% Compare results: prob of being significant vs #occur
load('simulate_results_shape.mat')
nSim = 100;
K = 10;
nCond = 3;
num_detected_SeqNMF = zeros(K,nCond);
num_significant_SeqNMF = zeros(K,nCond);
num_detected_FlexMF = zeros(K,nCond);
num_significant_FlexMF = zeros(K,nCond);
coeffs_SeqNMF = nan(nSim, K, nCond);
coeffs_FlexMF = nan(nSim, K, nCond);
coeffs_SeqNMF_significant = nan(nSim, K, nCond);
coeffs_FlexMF_significant = nan(nSim, K, nCond);

for i=1:nSim
    for j=1:nCond
        is_significant_SeqNMF = logical(is_significants_SeqNMF{i,j});
        is_significant_FlexMF = logical(is_significants_FlexMF{1,j});
        W = Ws{i,j};
        W_hat_SeqNMF = W_hats_SeqNMF{i,j};
        W_hat_FlexMF = W_hats_FlexMF{i,j};
        % pair with ground truth factor
        [coeff_SeqNMF, ids_SeqNMF] = helper.similarity_W(W, W_hat_SeqNMF);
        [coeff_FlexMF, ids_FlexMF] = helper.similarity_W(W, W_hat_FlexMF);
        num_detected_SeqNMF(ids_SeqNMF,j) = num_detected_SeqNMF(ids_SeqNMF,j)+1;
        num_detected_FlexMF(ids_FlexMF,j) = num_detected_FlexMF(ids_FlexMF,j)+1;
    
        ids_SeqNMF_significant = ids_SeqNMF(is_significant_SeqNMF(1:length(ids_SeqNMF)));
        ids_FlexMF_significant = ids_FlexMF(is_significant_FlexMF(1:length(ids_FlexMF)));
    
        num_significant_SeqNMF(ids_SeqNMF_significant) = num_significant_SeqNMF(ids_SeqNMF_significant)+1;
        num_significant_FlexMF(ids_FlexMF_significant) = num_significant_FlexMF(ids_FlexMF_significant)+1;
    
        coeffs_SeqNMF(i,ids_SeqNMF,j) = coeff_SeqNMF;
        coeffs_FlexMF(i,ids_FlexMF,j) = coeff_FlexMF;
    
        coeffs_SeqNMF_significant(i,ids_SeqNMF_significant,j) = coeff_SeqNMF(is_significant_SeqNMF(1:length(ids_SeqNMF)));
        coeffs_FlexMF_significant(i,ids_FlexMF_significant,j) = coeff_FlexMF(is_significant_FlexMF(1:length(ids_FlexMF)));
    end
end

%% Make plots
color_palet = [[1 .5 .5]; [1 .2 .2]; [.5 .5 1]; [.2 .2 1]]; 
figure; hold on
plot(1:K, num_detected_SeqNMF/nSim, 'LineWidth',2,'Color',color_palet(1,:))
plot(1:K, num_significant_SeqNMF/nSim, 'LineWidth',2,'Color',color_palet(2,:))
plot(1:K, num_detected_FlexMF/nSim, 'LineWidth',2,'Color',color_palet(3,:))
plot(1:K, num_significant_FlexMF/nSim, 'LineWidth',2,'Color',color_palet(4,:))
xlabel('# Motif Occurences')
ylabel('Proportion')
ylim([0,1])
legend({'SeqNMF Detected', 'SeqNMF Significant', 'FlexMF Detected', 'FlexMF Significant'}, 'Location', 'southeast')
set(gca, 'FontSize', 14)

%%
figure;
boxplot(coeffs_SeqNMF, 'Colors',color_palet(1,:), 'Plotstyle', 'compact', 'Positions', 2*(1:K)-.25);
hold on
boxplot(coeffs_FlexMF, 'Colors',color_palet(3,:), 'Plotstyle', 'compact', 'Positions', 2*(1:K)+.25);
xlabel('# Motif Occurences')
ylabel('Motif Similarity Coefficients')
ylim([0,1])
set(gca, 'FontSize', 12, 'XTickLabel', num2str([1:K]'), 'XTick', 2*(1:K))
boxes_SeqNMF = findobj(gca,'Tag','Box','Color',color_palet(1,:));
boxes_FlexMF = findobj(gca,'Tag','Box','Color',color_palet(3,:));
legend([boxes_SeqNMF(1), boxes_FlexMF(1)], {'SeqNMF Detected', 'FlexMF Detected'}, 'Location', 'bestoutside')

%%
figure;
boxplot(coeffs_SeqNMF_significant, 'Colors',color_palet(2,:), 'Plotstyle', 'compact', 'Positions', 2*(1:K)-.25);
hold on
boxplot(coeffs_FlexMF_significant, 'Colors',color_palet(4,:), 'Plotstyle', 'compact', 'Positions', 2*(1:K)+.25);
xlabel('# Motif Occurences')
ylabel('Motif Similarity Coefficients')
ylim([0,1])
set(gca, 'FontSize', 12, 'XTickLabel', num2str([1:K]'), 'XTick', 2*(1:K))
boxes_SeqNMF = findobj(gca,'Tag','Box','Color',color_palet(2,:));
boxes_FlexMF = findobj(gca,'Tag','Box','Color',color_palet(4,:));
legend([boxes_SeqNMF(1), boxes_FlexMF(1)], {'SeqNMF Significant', 'FlexMF Significant'}, 'Location', 'bestoutside')