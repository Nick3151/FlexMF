function [f1,f2,f3] = simulation_compare_algorithm_plot(file_name, noise_levels, noise_type, sel, varargin)
p = inputParser;
addOptional(p, 'nSim', 50);
addOptional(p, 'K', 10);
addOptional(p, 'c1', [0 0.4549 0.2902]);
addOptional(p, 'c2', [0.9020 0 0.7843]);
addOptional(p, "linestyles", [repmat("-",1,6), repmat("--",1,6)]);
addOptional(p, "markers", repmat(["o", "x", "s", "d", ">", "p"],1,2));
addOptional(p, 'thresh_x', 3.5)
addOptional(p, 'thresh_y', .8)
addOptional(p, 'thresh_reliable', .9:-.1:.7)
addOptional(p, 'lambda', .01)
addOptional(p, 'metric', 'corr')

parse(p, varargin{:});
nSim = p.Results.nSim;
K = p.Results.K;
c1 = p.Results.c1;
c2 = p.Results.c2;
linestyles = p.Results.linestyles;
markers = p.Results.markers;
thresh_x = p.Results.thresh_x;
thresh_y = p.Results.thresh_y;
thresh_reliable = p.Results.thresh_reliable;
lambda = p.Results.lambda;
metric = p.Results.metric;

load(file_name)
nCond = length(noise_levels);
nThresh = length(thresh_reliable);
if isnumeric(noise_levels)
    xticklabels = arrayfun(@num2str, noise_levels, 'UniformOutput', false);
else
    xticklabels = noise_levels;
end

num_success_SeqNMF = zeros(K,nCond);
num_success_SBI = zeros(K,nCond);
num_significant_SeqNMF = zeros(K,nCond);
num_significant_SBI = zeros(K,nCond);
min_occur_detect_SeqNMF = nan(nThresh,nCond);
min_occur_detect_SBI = nan(nThresh,nCond);
coeffs_W_SeqNMF = zeros(nSim, K, nCond);
coeffs_W_SBI = zeros(nSim, K, nCond);
coeffs_H_SeqNMF = zeros(nSim, K, nCond);
coeffs_H_SBI = zeros(nSim, K, nCond);
F1s_H_SeqNMF = zeros(nSim, K, nCond);
F1s_H_SBI = zeros(nSim, K, nCond);
sparsity_H_SeqNMF = zeros(nSim, nCond);
sparsity_W_SeqNMF = zeros(nSim, nCond);
sparsity_H_SBI = zeros(nSim, nCond);
sparsity_W_SBI = zeros(nSim, nCond);
sparsity_reg_SeqNMF = zeros(nSim, nCond);
sparsity_reg_SBI = zeros(nSim, nCond);
recon_errors_SeqNMF = zeros(nSim, nCond);
reg_costs_SeqNMF = zeros(nSim, nCond);
recon_errors_SBI = zeros(nSim, nCond);
reg_costs_SBI = zeros(nSim, nCond);
L = size(W_hats_SeqNMF{1,1}, 3);

for i=1:nSim
    for j=1:nCond
        is_significant_SeqNMF = logical(is_significants_SeqNMF{i,j});
        is_significant_FlexMF = logical(is_significants_FlexMF{1,j});
        W = Ws{i,j};
        H = Hs{i,j};
        W_hat_SeqNMF = W_hats_SeqNMF{i,j};
        W_hat_SBI = W_hats_FlexMF{i,j};
        H_hat_SeqNMF = H_hats_SeqNMF{i,j};
        H_hat_SBI = H_hats_FlexMF{i,j};
        % pair with ground truth factor
        [coeff_W_SeqNMF, coeff_H_SeqNMF, ids_SeqNMF] = helper.similarity_WH(W, H, W_hat_SeqNMF, H_hat_SeqNMF);
        [coeff_W_SBI, coeff_H_SBI, ids_SBI] = helper.similarity_WH(W, H, W_hat_SBI, H_hat_SBI);
        [S_W_SeqNMF, S_H_SeqNMF, ids_SeqNMF2, details_SeqNMF] = helper.similarity_WH_shape_activation(W, H, W_hat_SeqNMF, H_hat_SeqNMF);
        [S_W_SBI, S_H_SBI, ids_SBI2, details_SBI] = helper.similarity_WH_shape_activation(W, H, W_hat_SBI, H_hat_SBI);

        assert(isequal(ids_SeqNMF2, ids_SeqNMF), 'SeqNMF ids do not match!')
        assert(isequal(ids_SBI2, ids_SBI), 'SBI ids do not match!')

        ids_success_SeqNMF = ids_SeqNMF(is_significant_SeqNMF & (coeff_W_SeqNMF>.8));
        ids_success_SBI = ids_SBI(is_significant_FlexMF & (coeff_W_SBI>.8));
        ids_significant_SeqNMF = ids_SeqNMF(is_significant_SeqNMF);
        ids_significant_SBI = ids_SBI(is_significant_FlexMF);
        
        ids_success_SeqNMF(~ids_success_SeqNMF) = [];
        ids_success_SBI(~ids_success_SBI) = [];
        ids_significant_SeqNMF(~ids_significant_SeqNMF) = [];
        ids_significant_SBI(~ids_significant_SBI) = [];

        if ~isempty(ids_success_SeqNMF)
            num_success_SeqNMF(ids_success_SeqNMF,j) = num_success_SeqNMF(ids_success_SeqNMF,j)+1;
        end
        if ~isempty(ids_success_SBI)
            num_success_SBI(ids_success_SBI,j) = num_success_SBI(ids_success_SBI,j)+1;
        end
        if ~isempty(ids_significant_SeqNMF)
            num_significant_SeqNMF(ids_significant_SeqNMF,j) = num_significant_SeqNMF(ids_significant_SeqNMF,j)+1;
        end
        if ~isempty(ids_significant_SBI)
            num_significant_SBI(ids_significant_SBI,j) = num_significant_SBI(ids_significant_SBI,j)+1;
        end

        coeff_W_SeqNMF(~ids_SeqNMF) = [];
        coeff_H_SeqNMF(~ids_SeqNMF) = [];
        S_H_SeqNMF(~ids_SeqNMF) = [];
        ids_SeqNMF(~ids_SeqNMF) = [];
        coeff_W_SBI(~ids_SBI) = [];
        coeff_H_SBI(~ids_SBI) = [];
        S_H_SBI(~ids_SBI) = [];
        ids_SBI(~ids_SBI) = [];
    
        coeffs_W_SeqNMF(i,ids_SeqNMF,j) = coeff_W_SeqNMF;
        coeffs_W_SBI(i,ids_SBI,j) = coeff_W_SBI;
        coeffs_H_SeqNMF(i,ids_SeqNMF,j) = coeff_H_SeqNMF;
        coeffs_H_SBI(i,ids_SBI,j) = coeff_H_SBI;
        F1s_H_SeqNMF(i,ids_SeqNMF,j) = S_H_SeqNMF;
        F1s_H_SBI(i,ids_SBI,j) = S_H_SBI;

        sparsity_W_SeqNMF(i,j) = mean(W_hat_SeqNMF>1e-3, "all");
        sparsity_H_SeqNMF(i,j) = mean(H_hat_SeqNMF>1e-3, "all");
        sparsity_W_SBI(i,j) = mean(W_hat_SBI>1e-3, "all");
        sparsity_H_SBI(i,j) = mean(H_hat_SBI>1e-3, "all");

        TrainingData = TrainingDatas{i,j};
        
%         % Sparsity of regularization term
%         Q = ones(K);
%         Q(1:K+1:end) = 0;   % off diagonal mask
%         smoothkernel = ones(1,(2*L)-1); 
%         WTX_SeqNMF = helper.transconv(W_hat_SeqNMF, TrainingData);
%         WTXS_SeqNMF = conv2(abs(WTX_SeqNMF), smoothkernel, 'same');
%         reg_SeqNMF = WTXS_SeqNMF*H_hat_SeqNMF';
%         sparsity_reg_SeqNMF(i,j) = mean(reg_SeqNMF.*Q>1e-3, "all");
%         WTX_FlexMF = helper.transconv(W_hat_FlexMF, TrainingData);
%         WTXS_FlexMF = conv2(abs(WTX_FlexMF), smoothkernel, 'same');
%         reg_FlexMF = WTXS_FlexMF*H_hat_SeqNMF';
%         sparsity_reg_FlexMF(i,j) = mean(reg_FlexMF.*Q>1e-3, "all");

        [recon_errors_SeqNMF(i,j), reg_cross, ~, ~] = helper.get_FlexMF_cost(TrainingData,W_hat_SeqNMF,H_hat_SeqNMF);
        reg_costs_SeqNMF(i,j) = reg_cross*lambda;
        [recon_errors_SBI(i,j), reg_cross, ~, ~] = helper.get_FlexMF_cost(TrainingData,W_hat_SBI,H_hat_SBI);
        reg_costs_SBI(i,j) = reg_cross*lambda;
    end
end

% Minimum motif occurrence to be reliably detected
for i=1:nThresh
    for j=1:nCond
        can_reliably_detect = num_success_SeqNMF(:,j)/nSim>=thresh_reliable(i);
        if any(can_reliably_detect)
            min_occur_detect_SeqNMF(i,j) = find(can_reliably_detect, 1);
        end
        can_reliably_detect = num_success_SBI(:,j)/nSim>=thresh_reliable(i);
        if any(can_reliably_detect)
            min_occur_detect_SBI(i,j) = find(can_reliably_detect, 1);
        end
    end
end

%% Make plots, motif coefficients to ground truth
f1 = figure; 
legends_MR = cellfun(@(x) [noise_type, '=', x, ' MUR'], xticklabels,'UniformOutput',false);
legends_SB = cellfun(@(x) [noise_type, '=', x ' SBI'], xticklabels,'UniformOutput',false);   
ax1 = subplot('Position', [0.1 0.75 0.6 0.2]);
hold on
for j=1:nCond
    p1(j) = plot(2*(1:K), num_success_SeqNMF(:,j)./nSim, 'LineStyle', linestyles(j), 'Marker', markers(j), 'LineWidth',2,'Color',c1);
%     p2(j) = plot(2*(1:K), num_success_SBI(:,j)./nSim, 'LineStyle', linestyles(j), 'Marker', markers(j), 'LineWidth',2,'Color',c2);
%     lh(j).Color(4) = alphas(j);
end
patch(2*[thresh_x, K+.5, K+.5, thresh_x], [0, 0, 1, 1], [.5 .5 .5], 'EdgeColor', 'none', 'FaceAlpha', .3)
set(gca, 'FontSize', 12, 'XTick', [])

ax2 = subplot('Position', [0.1 0.5 0.6 0.2]);
hold on
for j=1:nCond
%     p1(j) = plot(2*(1:K), num_success_SeqNMF(:,j)./nSim, 'LineStyle', linestyles(j), 'Marker', markers(j), 'LineWidth',2,'Color',c1);
    p2(j) = plot(2*(1:K), num_success_SBI(:,j)./nSim, 'LineStyle', linestyles(j), 'Marker', markers(j), 'LineWidth',2,'Color',c2);
%     lh(j).Color(4) = alphas(j);
end

% ylabel('P(significant)')
patch(2*[thresh_x, K+.5, K+.5, thresh_x], [0, 0, 1, 1], [.5 .5 .5], 'EdgeColor', 'none', 'FaceAlpha', .3)
legend([p1 p2], [legends_MR legends_SB], 'Orientation','vertical', 'NumColumns',1, 'Position',[0.7 0.1 0.25 0.8], 'Box','off')
set(gca, 'FontSize', 12, 'XTick', [])

j = sel;
ax3 = subplot('Position', [0.1 0.3 0.6 0.15]);
% ax1 = subplot(211);
hold on
% boxplot(coeffs_W_SeqNMF(:,:,j), 'Colors', c1, 'Positions', 2*(1:K)-.25, 'Jitter', .1, 'Widths', .2);
% boxplot(coeffs_W_FlexMF(:,:,j), 'Colors', c2, 'Positions', 2*(1:K)+.25, 'Jitter', .1, 'Widths', .2);
x = repmat(2*(1:K), [nSim,1]);
swarmchart(x-.25, coeffs_W_SeqNMF(:,:,j), 12, c1, 'filled', 'XJitterWidth', .2);
swarmchart(x+.25, coeffs_W_SBI(:,:,j), 12, c2, 'filled', 'XJitterWidth', .2);
% title('Correlation Coeff Multiplicative Update Rule')
patch(2*[thresh_x, K+.5, K+.5, thresh_x], [0, 0, 1, 1], [.5 .5 .5], 'EdgeColor', 'none', 'FaceAlpha', .3)
yline(thresh_y, 'LineWidth', 2, 'LineStyle','--');
set(gca, 'FontSize', 12, 'XTick', 2*(1:K), 'XTickLabel', [], 'YTick', [0, .5, thresh_y, 1], 'box','off')

if strcmp(metric, 'corr')
    ax4 = subplot('Position', [0.1 0.1 0.6 0.15]);
    hold on
    s1 = swarmchart(x-.25, coeffs_H_SeqNMF(:,:,j), 12, c1, 'filled', 'XJitterWidth', .2);
    s2 = swarmchart(x+.25, coeffs_H_SBI(:,:,j), 12, c2, 'filled', 'XJitterWidth', .2);
    % title('Correlation Coeff Split Bregman Iteration')
    set(gca, 'FontSize', 12, 'XTick', 2*(1:K), 'XTickLabel', num2str((1:K)'), 'YTick', [0, .5, thresh_y, 1], 'box','off')
    patch(2*[thresh_x, K+.5, K+.5, thresh_x], [0, 0, 1, 1], [.5 .5 .5], 'EdgeColor', 'none', 'FaceAlpha', .3)
    yline(thresh_y, 'LineWidth', 2, 'LineStyle','--');
    % legend([s1(1), s2(1)], {'MUR', 'SBI'}, 'box','off')
    % xtickangle(ax, 90)
    % xlabel('# Motif Occurences')
elseif strcmp(metric, 'F1')
    ax4 = subplot('Position', [0.1 0.1 0.6 0.15]);
    hold on
    s1 = swarmchart(x-.25, F1s_H_SeqNMF(:,:,j), 12, c1, 'filled', 'XJitterWidth', .2);
    s2 = swarmchart(x+.25, F1s_H_SBI(:,:,j), 12, c2, 'filled', 'XJitterWidth', .2);
    % title('Correlation Coeff Split Bregman Iteration')
    set(gca, 'FontSize', 12, 'XTick', 2*(1:K), 'XTickLabel', num2str((1:K)'), 'YTick', [0, .5, thresh_y, 1], 'box','off')
    patch(2*[thresh_x, K+.5, K+.5, thresh_x], [0, 0, 1, 1], [.5 .5 .5], 'EdgeColor', 'none', 'FaceAlpha', .3)
    yline(thresh_y, 'LineWidth', 2, 'LineStyle','--');
else
    error('metric must be corr or F1')
end

set(f1, 'Position', [100, 100, 1000, 1000])
linkaxes([ax1, ax2, ax3, ax4], 'xy')
set(gca, 'XLim', 2*[0, K+.5], 'YLim', [0,1])

%% Compare sparsity levels
f2 = figure;
subplot(211)
hold on
x = repmat(2*(1:nCond), [nSim,1]);
swarmchart(x-.25, sparsity_H_SeqNMF, 12, c1, 'filled', 'XJitterWidth', .2);
swarmchart(x+.25, sparsity_H_SBI, 12, c2, 'filled', 'XJitterWidth', .2);
% boxplot(sparsity_H_SeqNMF, 'Colors', c1, 'Positions', 2*(1:nCond)-.25, 'Jitter', .2, 'Widths',.2);
% boxplot(sparsity_H_FlexMF, 'Colors', c2, 'Positions', 2*(1:nCond)+.25, 'Jitter', .2, 'Widths',.2);
title('Sparsity of H')
ylim([0,.2])
set(gca, 'FontSize', 12, 'XTickLabel', xticklabels, 'XTick', 2*(1:nCond), 'box','off')
xtickangle(45)

subplot(212)
hold on
% boxplot(sparsity_W_SeqNMF, 'Colors', c1, 'Positions', 2*(1:nCond)-.25, 'Jitter', .2, 'Widths',.2);
% boxplot(sparsity_W_FlexMF, 'Colors', c2, 'Positions', 2*(1:nCond)+.25, 'Jitter', .2, 'Widths',.2);
s1 = swarmchart(x-.25, sparsity_W_SeqNMF, 12, c1, 'filled', 'XJitterWidth', .2);
s2 = swarmchart(x+.25, sparsity_W_SBI, 12, c2, 'filled', 'XJitterWidth', .2);
xlabel(noise_type)
title('Sparsity of W')
ylim([0,.2])
set(gca, 'FontSize', 12, 'XTickLabel', xticklabels, 'XTick', 2*(1:nCond), 'box','off')
xtickangle(45)

% subplot(133)
% hold on
% boxplot(sparsity_reg_SeqNMF, 'Plotstyle', 'compact', 'Colors', c1, 'Positions', 2*(1:nCond)-.25, 'Jitter', 0);
% boxplot(sparsity_reg_FlexMF, 'Plotstyle', 'compact', 'Colors', c2, 'Positions', 2*(1:nCond)+.25, 'Jitter', 0);
% title('Sparsity of Regularization')
% ylim([0,1])
% set(gca, 'FontSize', 12, 'XTickLabel', legends, 'XTick', 2*(1:nCond), 'YTick', [], 'box','off', 'Position', [0.55 0.15 0.2 0.75])
% xtickangle(45)

% s1 = findobj(gca,'Tag','Scatter','Color',c1);
% s2 = findobj(gca,'Tag','Scatter','Color',c2);
legend([s1(1), s2(1)], {'MUR', 'SBI'}, 'box','off')
set(f2, 'Position', [100, 100, 1000, 800])

%% Compare cost functions
f3 = figure; 
ax1 = subplot(211);
hold on
swarmchart(x-.25, recon_errors_SeqNMF, 12, c1, 'filled', 'XJitterWidth', .2);
swarmchart(x+.25, recon_errors_SBI, 12, c2, 'filled', 'XJitterWidth', .2);
title('Reconstruction Costs')
set(gca, 'FontSize', 12, 'XTickLabel', xticklabels, 'XTick', 2*(1:nCond), 'box','off')

ax2 = subplot(212);
hold on
s1 = swarmchart(x-.25, reg_costs_SeqNMF, 12, c1, 'filled', 'XJitterWidth', .2);
s2 = swarmchart(x+.25, reg_costs_SBI, 12, c2, 'filled', 'XJitterWidth', .2);
xlabel(noise_type)
title('Regularization Costs')
set(gca, 'FontSize', 12, 'XTickLabel', xticklabels, 'XTick', 2*(1:nCond), 'box','off')
linkaxes([ax1, ax2], 'y')
legend([s1(1), s2(1)], {'MUR', 'SBI'}, 'box','off')
set(f3, 'Position', [100, 100, 1000, 800])

%% Minimum motif occurrence for reliable detection
% f4 = figure;
% legends_MR = arrayfun(@(x) sprintf('%d%% detected MUR', round(x*100)), thresh_reliable, 'UniformOutput', false);
% legends_SB = arrayfun(@(x) sprintf('%d%% detected SBI', round(x*100)), thresh_reliable, 'UniformOutput', false);  
% hold on 
% clear p1 p2;
% for i=1:nThresh
%     p1(i) = plot(1:nCond, min_occur_detect_SeqNMF(i,:), 'Color', c1, 'LineStyle', linestyles(i), 'Marker', markers(i));
%     p2(i) = plot(1:nCond, min_occur_detect_SBI(i,:), 'Color', c2, 'LineStyle', linestyles(i), 'Marker', markers(i));
% end
% set(gca, 'FontSize', 12, 'XTickLabel', xticklabels, 'XTick', 1:nCond, 'box','off')
% xlabel(noise_type)
% ylim([0,K])
% xlim([1,nCond])
% legend([p1 p2], [legends_MR legends_SB], 'box','off', 'Location','best')
