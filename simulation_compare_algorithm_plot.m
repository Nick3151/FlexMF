function [f1,f2,f3] = simulation_compare_algorithm_plot(file_name, nCond, legends, sel, varargin)
p = inputParser;
addOptional(p, 'nSim', 100);
addOptional(p, 'K', 10);
addOptional(p, 'c1', [0 0.4549 0.2902]);
addOptional(p, 'c2', [0.9020 0 0.7843]);
addOptional(p, "linestyles", ["-", "--", ":", "-.", "-"]);
addOptional(p, "markers", [".", ".", ".", "." "o"]);
addOptional(p, 'thresh_x', 3.5)
addOptional(p, 'thresh_y', .8)
parse(p, varargin{:});
nSim = p.Results.nSim;
K = p.Results.K;
c1 = p.Results.c1;
c2 = p.Results.c2;
linestyles = p.Results.linestyles;
markers = p.Results.markers;
thresh_x = p.Results.thresh_x;
thresh_y = p.Results.thresh_y;

load(file_name)

num_success_SeqNMF = zeros(K,nCond);
num_success_FlexMF = zeros(K,nCond);
coeffs_W_SeqNMF = nan(nSim, K, nCond);
coeffs_W_FlexMF = nan(nSim, K, nCond);
coeffs_H_SeqNMF = nan(nSim, K, nCond);
coeffs_H_FlexMF = nan(nSim, K, nCond);
sparsity_H_SeqNMF = zeros(nSim, nCond);
sparsity_W_SeqNMF = zeros(nSim, nCond);
sparsity_H_FlexMF = zeros(nSim, nCond);
sparsity_W_FlexMF = zeros(nSim, nCond);
sparsity_reg_SeqNMF = zeros(nSim, nCond);
sparsity_reg_FlexMF = zeros(nSim, nCond);
recon_errors_SeqNMF = zeros(nSim, nCond);
reg_costs_SeqNMF = zeros(nSim, nCond);
recon_errors_FlexMF = zeros(nSim, nCond);
reg_costs_FlexMF = zeros(nSim, nCond);
L = size(W_hats_SeqNMF{1,1}, 3);

for i=1:nSim
    for j=1:nCond
        is_significant_SeqNMF = logical(is_significants_SeqNMF{i,j});
        is_significant_FlexMF = logical(is_significants_FlexMF{1,j});
        W = Ws{i,j};
        H = Hs{i,j};
        W_hat_SeqNMF = W_hats_SeqNMF{i,j};
        W_hat_FlexMF = W_hats_FlexMF{i,j};
        H_hat_SeqNMF = H_hats_SeqNMF{i,j};
        H_hat_FlexMF = H_hats_FlexMF{i,j};
        % pair with ground truth factor
        [coeff_W_SeqNMF, coeff_H_SeqNMF, ids_SeqNMF] = helper.similarity_WH(W, H, W_hat_SeqNMF, H_hat_SeqNMF);
        [coeff_W_FlexMF, coeff_H_FlexMF, ids_FlexMF] = helper.similarity_WH(W, H, W_hat_FlexMF, H_hat_FlexMF);
        ids_success_SeqNMF = ids_SeqNMF(is_significant_SeqNMF(1:length(ids_SeqNMF)) & (coeff_W_SeqNMF>.8) & (coeff_H_SeqNMF>.8));
        ids_success_FlexMF = ids_FlexMF(is_significant_FlexMF(1:length(ids_FlexMF)) & (coeff_W_FlexMF>.8) & (coeff_H_FlexMF>.8));
    
        num_success_SeqNMF(ids_success_SeqNMF,j) = num_success_SeqNMF(ids_success_SeqNMF,j)+1;
        num_success_FlexMF(ids_success_FlexMF,j) = num_success_FlexMF(ids_success_FlexMF,j)+1;
    
        coeffs_W_SeqNMF(i,ids_SeqNMF,j) = coeff_W_SeqNMF;
        coeffs_W_FlexMF(i,ids_FlexMF,j) = coeff_W_FlexMF;
        coeffs_H_SeqNMF(i,ids_SeqNMF,j) = coeff_H_SeqNMF;
        coeffs_H_FlexMF(i,ids_FlexMF,j) = coeff_H_FlexMF;

        sparsity_W_SeqNMF(i,j) = mean(W_hat_SeqNMF==0, "all");
        sparsity_H_SeqNMF(i,j) = mean(H_hat_SeqNMF==0, "all");
        sparsity_W_FlexMF(i,j) = mean(W_hat_FlexMF==0, "all");
        sparsity_H_FlexMF(i,j) = mean(H_hat_FlexMF==0, "all");
        
        % Sparsity of regularization term
        TrainingData = TrainingDatas{i,j};
        Q = ones(K);
        Q(1:K+1:end) = 0;   % off diagonal mask
        smoothkernel = ones(1,(2*L)-1); 
        WTX_SeqNMF = helper.transconv(W_hat_SeqNMF, TrainingData);
        WTXS_SeqNMF = conv2(abs(WTX_SeqNMF), smoothkernel, 'same');
        reg_SeqNMF = WTXS_SeqNMF*H_hat_SeqNMF';
        sparsity_reg_SeqNMF(i,j) = mean(reg_SeqNMF.*Q==0, "all");
        WTX_FlexMF = helper.transconv(W_hat_FlexMF, TrainingData);
        WTXS_FlexMF = conv2(abs(WTX_FlexMF), smoothkernel, 'same');
        reg_FlexMF = WTXS_FlexMF*H_hat_SeqNMF';
        sparsity_reg_FlexMF(i,j) = mean(reg_FlexMF.*Q==0, "all");

        lambda=.005;
        [recon_errors_SeqNMF(i,j), reg_cross, ~, ~] = helper.get_FlexMF_cost(TrainingData,W_hat_SeqNMF,H_hat_SeqNMF);
        reg_costs_SeqNMF(i,j) = reg_cross*lambda;
        [recon_errors_FlexMF(i,j), reg_cross, ~, ~] = helper.get_FlexMF_cost(TrainingData,W_hat_FlexMF,H_hat_FlexMF);
        reg_costs_FlexMF(i,j) = reg_cross*lambda;
    end
end

% Make plots
legends_MR = cellfun(@(x) [x ' MUR'], legends,'UniformOutput',false);
legends_SB = cellfun(@(x) [x ' SBI'], legends,'UniformOutput',false);   

f1 = figure; 
ax1 = subplot('Position', [0.1 0.6 0.6 0.3]);
hold on
for j=1:nCond
    p1(j) = plot(2*(1:K), num_success_SeqNMF(:,j)./nSim, 'LineStyle', linestyles(j), 'Marker', markers(j), 'LineWidth',2,'Color',c1);
    p2(j) = plot(2*(1:K), num_success_FlexMF(:,j)./nSim, 'LineStyle', linestyles(j), 'Marker', markers(j), 'LineWidth',2,'Color',c2);
%     lh(j).Color(4) = alphas(j);
end

% ylabel('P(significant)')
patch(2*[thresh_x, K+.5, K+.5, thresh_x], [0, 0, 1, 1], [.5 .5 .5], 'EdgeColor', 'none', 'FaceAlpha', .3)
legend([p1 p2], [legends_MR legends_SB], 'Orientation','vertical', 'NumColumns',1, 'Position',[0.7 0.1 0.25 0.8], 'Box','off')
set(gca, 'FontSize', 12, 'XTick', [])

j = sel;
ax2 = subplot('Position', [0.1 0.35 0.6 0.2]);
hold on
boxplot(coeffs_W_SeqNMF(:,:,j), 'Plotstyle', 'compact', 'Colors', c1, 'Positions', 2*(1:K)-.25, 'Jitter', 0);
boxplot(coeffs_W_FlexMF(:,:,j), 'Plotstyle', 'compact', 'Colors', c2, 'Positions', 2*(1:K)+.25, 'Jitter', 0);
% title('Correlation Coeff Multiplicative Update Rule')
patch(2*[thresh_x, K+.5, K+.5, thresh_x], [0, 0, 1, 1], [.5 .5 .5], 'EdgeColor', 'none', 'FaceAlpha', .3)
yline(thresh_y, 'LineWidth', 2, 'LineStyle','--');
set(gca, 'FontSize', 12, 'XTick', 2*(1:K), 'YTick', [0, .5, thresh_y, 1], 'box','off', 'Position', [0.1 0.35 0.6 0.2])

ax3 = subplot('Position', [0.1 0.1 0.6 0.2]);
hold on
boxplot(coeffs_H_SeqNMF(:,:,j), 'Plotstyle', 'compact', 'Colors', c1, 'Positions', 2*(1:K)-.25, 'Jitter', 0);
boxplot(coeffs_H_FlexMF(:,:,j), 'Plotstyle', 'compact', 'Colors', c2, 'Positions', 2*(1:K)+.25, 'Jitter', 0);
% title('Correlation Coeff Split Bregman Iteration')
set(gca, 'FontSize', 12, 'XTick', 2*(1:K), 'XTickLabel', num2str((1:K)'), 'YTick', [0, .5, thresh_y, 1], 'box','off', 'Position', [0.1 0.1 0.6 0.2])
patch(2*[thresh_x, K+.5, K+.5, thresh_x], [0, 0, 1, 1], [.5 .5 .5], 'EdgeColor', 'none', 'FaceAlpha', .3)
yline(thresh_y, 'LineWidth', 2, 'LineStyle','--');
% xtickangle(ax, 90)
% xlabel('# Motif Occurences')
set(gcf, 'Position', [100, 100, 800, 800])
linkaxes([ax1, ax2, ax3], 'xy')
set(gca, 'XLim', 2*[0, K+.5], 'YLim', [0,1])


% Compare sparsity levels
f2 = figure;
subplot(131)
hold on
boxplot(sparsity_H_SeqNMF, 'Plotstyle', 'compact', 'Colors', c1, 'Positions', 2*(1:nCond)-.25, 'Jitter', 0);
boxplot(sparsity_H_FlexMF, 'Plotstyle', 'compact', 'Colors', c2, 'Positions', 2*(1:nCond)+.25, 'Jitter', 0);
title('Sparsity of H')
ylim([0,1])
set(gca, 'FontSize', 12, 'XTickLabel', legends, 'XTick', 2*(1:nCond), 'box','off', 'Position', [0.05 0.15 0.2 0.75])
xtickangle(45)

subplot(132)
hold on
boxplot(sparsity_W_SeqNMF, 'Plotstyle', 'compact', 'Colors', c1, 'Positions', 2*(1:nCond)-.25, 'Jitter', 0);
boxplot(sparsity_W_FlexMF, 'Plotstyle', 'compact', 'Colors', c2, 'Positions', 2*(1:nCond)+.25, 'Jitter', 0);
title('Sparsity of W')
ylim([0,1])
set(gca, 'FontSize', 12, 'XTickLabel', legends, 'XTick', 2*(1:nCond), 'YTick', [], 'box','off', 'Position', [0.3 0.15 0.2 0.75])
xtickangle(45)

subplot(133)
hold on
boxplot(sparsity_reg_SeqNMF, 'Plotstyle', 'compact', 'Colors', c1, 'Positions', 2*(1:nCond)-.25, 'Jitter', 0);
boxplot(sparsity_reg_FlexMF, 'Plotstyle', 'compact', 'Colors', c2, 'Positions', 2*(1:nCond)+.25, 'Jitter', 0);
title('Sparsity of Regularization')
ylim([0,1])
set(gca, 'FontSize', 12, 'XTickLabel', legends, 'XTick', 2*(1:nCond), 'YTick', [], 'box','off', 'Position', [0.55 0.15 0.2 0.75])
xtickangle(45)

boxes_SeqNMF = findobj(gca,'Tag','Box','Color',c1);
boxes_FlexMF = findobj(gca,'Tag','Box','Color',c2);
legend([boxes_SeqNMF(1), boxes_FlexMF(1)], {'MUR', 'SBI'}, 'box','off', 'Position', [0.8 0.75 0.15 0.1])
set(gcf, 'Position', [100, 100, 1500, 800])

% Compare cost functions
f3 = figure; 
subplot(121)
hold on
boxplot(recon_errors_SeqNMF, 'Plotstyle', 'compact', 'Colors', c1, 'Positions', 2*(1:nCond)-.25, 'Jitter', 0);
boxplot(recon_errors_FlexMF, 'Plotstyle', 'compact', 'Colors', c2, 'Positions', 2*(1:nCond)+.25, 'Jitter', 0);
title('Reconstruction Costs')
set(gca, 'FontSize', 12, 'XTickLabel', legends, 'XTick', 2*(1:nCond), 'box','off', 'Position', [0.05 0.15 0.3 0.75])
xtickangle(45)

subplot(122)
hold on
boxplot(reg_costs_SeqNMF, 'Plotstyle', 'compact', 'Colors', c1, 'Positions', 2*(1:nCond)-.25, 'Jitter', 0);
boxplot(reg_costs_FlexMF, 'Plotstyle', 'compact', 'Colors', c2, 'Positions', 2*(1:nCond)+.25, 'Jitter', 0);
title('Regularization Costs')
set(gca, 'FontSize', 12, 'XTickLabel', legends, 'XTick', 2*(1:nCond), 'YTick', [], 'box','off', 'Position', [0.4 0.15 0.3 0.75])
xtickangle(45)
boxes_SeqNMF = findobj(gca,'Tag','Box','Color',c1);
boxes_FlexMF = findobj(gca,'Tag','Box','Color',c2);
legend([boxes_SeqNMF(1), boxes_FlexMF(1)], {'MUR', 'SBI'}, 'box','off', 'Position', [0.75 0.75 0.2 0.1])
set(gcf, 'Position', [100, 100, 1500, 800])