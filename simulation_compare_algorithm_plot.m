function simulation_compare_algorithm_plot(file_name, nCond, legends, sel, varargin)
p = inputParser;
addOptional(p, 'nSim', 100);
addOptional(p, 'K', 10);
addOptional(p, 'c1', 'b');
addOptional(p, 'c2', 'r');
addOptional(p, "linestyles", ["-", "--", ":", "-.", "-"]);
addOptional(p, "markers", [".", ".", ".", "." "o"]);
parse(p, varargin{:});
nSim = p.Results.nSim;
K = p.Results.K;
c1 = p.Results.c1;
c2 = p.Results.c2;
linestyles = p.Results.linestyles;
markers = p.Results.markers;

load(file_name)

num_detected_SeqNMF = zeros(K,nCond);
num_significant_SeqNMF = zeros(K,nCond);
num_detected_FlexMF = zeros(K,nCond);
num_significant_FlexMF = zeros(K,nCond);
coeffs_SeqNMF = nan(nSim, K, nCond);
coeffs_FlexMF = nan(nSim, K, nCond);
coeffs_SeqNMF_significant = nan(nSim, K, nCond);
coeffs_FlexMF_significant = nan(nSim, K, nCond);
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
        W_hat_SeqNMF = W_hats_SeqNMF{i,j};
        W_hat_FlexMF = W_hats_FlexMF{i,j};
        H_hat_SeqNMF = H_hats_SeqNMF{i,j};
        H_hat_FlexMF = H_hats_FlexMF{i,j};
        % pair with ground truth factor
        [coeff_SeqNMF, ids_SeqNMF] = helper.similarity_W(W, W_hat_SeqNMF);
        [coeff_FlexMF, ids_FlexMF] = helper.similarity_W(W, W_hat_FlexMF);
        num_detected_SeqNMF(ids_SeqNMF,j) = num_detected_SeqNMF(ids_SeqNMF,j)+1;
        num_detected_FlexMF(ids_FlexMF,j) = num_detected_FlexMF(ids_FlexMF,j)+1;
    
        ids_SeqNMF_significant = ids_SeqNMF(is_significant_SeqNMF(1:length(ids_SeqNMF)));
        ids_FlexMF_significant = ids_FlexMF(is_significant_FlexMF(1:length(ids_FlexMF)));
    
        num_significant_SeqNMF(ids_SeqNMF_significant,j) = num_significant_SeqNMF(ids_SeqNMF_significant,j)+1;
        num_significant_FlexMF(ids_FlexMF_significant,j) = num_significant_FlexMF(ids_FlexMF_significant,j)+1;
    
        coeffs_SeqNMF(i,ids_SeqNMF,j) = coeff_SeqNMF;
        coeffs_FlexMF(i,ids_FlexMF,j) = coeff_FlexMF;
    
        coeffs_SeqNMF_significant(i,ids_SeqNMF_significant,j) = coeff_SeqNMF(is_significant_SeqNMF(1:length(ids_SeqNMF)));
        coeffs_FlexMF_significant(i,ids_FlexMF_significant,j) = coeff_FlexMF(is_significant_FlexMF(1:length(ids_FlexMF)));

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
figure; 
ax(1) = subplot('Position', [0.1 0.6 0.8 0.3]);
hold on
for j=1:nCond
    plot(1:K, num_significant_SeqNMF(:,j)./nSim, 'LineStyle', linestyles(j), 'Marker', markers(j), 'LineWidth',2,'Color',c1);
    plot(1:K, num_significant_FlexMF(:,j)./nSim, 'LineStyle', linestyles(j), 'Marker', markers(j), 'LineWidth',2,'Color',c2);
%     lh(j).Color(4) = alphas(j);
end

ylim([0,1])
title('Probability of being significant')
% legend({'Noise=.001', 'Noise=.01', 'Noise=.1'}, 'Location', 'bestoutside')
set(gca, 'FontSize', 12, 'XTick', [])

j = sel;
ax(2) = subplot('Position', [0.1 0.35 0.8 0.2]);
hold on
boxplot(coeffs_SeqNMF(:,:,j), 'Plotstyle', 'compact', 'Colors', c1);
title('Correlation Coeff Multiplication Rule')
set(gca, 'FontSize', 12, 'XTick', [], 'box','off', 'Position', [0.1 0.35 0.8 0.2])

subplot('Position', [0.1 0.1 0.8 0.2]);
hold on
boxplot(coeffs_FlexMF(:,:,j), 'Plotstyle', 'compact', 'Colors',c2);
title('Correlation Coeff Split Bregman Iteration')
set(gca, 'FontSize', 12, 'XTick', 1:K, 'XTickLabel', num2str([1:K]'), 'box','off', 'Position', [0.1 0.1 0.8 0.2])
% xtickangle(ax, 90)
xlabel('# Motif Occurences')
set(gcf, 'Position', [100, 100, 900, 600])


% Compare sparsity levels
figure;
subplot(131)
hold on
boxplot(sparsity_H_SeqNMF, 'Plotstyle', 'compact', 'Colors', c1, 'Positions', 2*(1:nCond)-.25);
boxplot(sparsity_H_FlexMF, 'Plotstyle', 'compact', 'Colors', c2, 'Positions', 2*(1:nCond)+.25);
title('Sparsity of H')
ylim([0,1])
set(gca, 'FontSize', 12, 'XTickLabel', legends, 'XTick', 2*(1:nCond), 'box','off', 'Position', [0.05 0.15 0.2 0.75])
xtickangle(45)

subplot(132)
hold on
boxplot(sparsity_W_SeqNMF, 'Plotstyle', 'compact', 'Colors', c1, 'Positions', 2*(1:nCond)-.25);
boxplot(sparsity_W_FlexMF, 'Plotstyle', 'compact', 'Colors', c2, 'Positions', 2*(1:nCond)+.25);
title('Sparsity of W')
ylim([0,1])
set(gca, 'FontSize', 12, 'XTickLabel', legends, 'XTick', 2*(1:nCond), 'YTick', [], 'box','off', 'Position', [0.3 0.15 0.2 0.75])
xtickangle(45)

subplot(133)
hold on
boxplot(sparsity_reg_SeqNMF, 'Plotstyle', 'compact', 'Colors', c1, 'Positions', 2*(1:nCond)-.25);
boxplot(sparsity_reg_FlexMF, 'Plotstyle', 'compact', 'Colors', c2, 'Positions', 2*(1:nCond)+.25);
title('Sparsity of Regularization')
ylim([0,1])
set(gca, 'FontSize', 12, 'XTickLabel', legends, 'XTick', 2*(1:nCond), 'YTick', [], 'box','off', 'Position', [0.55 0.15 0.2 0.75])
xtickangle(45)

boxes_SeqNMF = findobj(gca,'Tag','Box','Color',c1);
boxes_FlexMF = findobj(gca,'Tag','Box','Color',c2);
legend([boxes_SeqNMF(1), boxes_FlexMF(1)], {'Multiplication Rule', 'Split Bregman Iteration'}, 'Position', [0.75 0.75 0.2 0.1])
set(gcf, 'Position', [100, 100, 900, 600])

% Compare cost functions
figure; 
subplot(121)
hold on
boxplot(recon_errors_SeqNMF, 'Plotstyle', 'compact', 'Colors', c1, 'Positions', 2*(1:nCond)-.25);
boxplot(recon_errors_FlexMF, 'Plotstyle', 'compact', 'Colors', c2, 'Positions', 2*(1:nCond)+.25);
title('Reconstruction Costs')
set(gca, 'FontSize', 12, 'XTickLabel', legends, 'XTick', 2*(1:nCond), 'box','off', 'Position', [0.05 0.15 0.3 0.75])
xtickangle(45)

subplot(122)
hold on
boxplot(reg_costs_SeqNMF, 'Plotstyle', 'compact', 'Colors', c1, 'Positions', 2*(1:nCond)-.25);
boxplot(reg_costs_FlexMF, 'Plotstyle', 'compact', 'Colors', c2, 'Positions', 2*(1:nCond)+.25);
title('Regularization Costs')
set(gca, 'FontSize', 12, 'XTickLabel', legends, 'XTick', 2*(1:nCond), 'YTick', [], 'box','off', 'Position', [0.4 0.15 0.3 0.75])
xtickangle(45)
boxes_SeqNMF = findobj(gca,'Tag','Box','Color',c1);
boxes_FlexMF = findobj(gca,'Tag','Box','Color',c2);
legend([boxes_SeqNMF(1), boxes_FlexMF(1)], {'Multiplication Rule', 'Split Bregman Iteration'}, 'Position', [0.75 0.75 0.2 0.1])
set(gcf, 'Position', [100, 100, 900, 600])