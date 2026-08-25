clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))

data_type = 'warpnoise';
lambda_R = 1;
lambda = 1e-1;

% load(fullfile('Simulation_Results', sprintf('FlexMF_EMD_results_%s_lambda=%0.3e.mat', data_type, lambda)))
load(fullfile('Simulation_Results', sprintf('FlexMF_EMD_results_%s_lambda=%0.3e.mat', data_type, lambda)))
[nMs, nSim] = size(H_hats);

%% Plot synthetic data
j = 5;
H = Hs{j};
W = Ws{j};
X = Xs{j};
T = size(H,2);
Htrain = H(:,1:round(T/2));
Xtrain = X(:,1:round(T/2));

plotAll = 1;
figure; SimpleWHPlot(W,Htrain,'Data',Xtrain, 'plotAll', plotAll); title('generated data raw','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
save2pdf(fullfile('Simulation_Results',sprintf('Synthetic_data_%s.pdf', data_type)))
%% Look at factors
Mi = 2;
j = 5;
plotAll = 1;

figure; SimpleWHPlot_patch(W_hats{Mi,j}, H_hats{Mi,j}, 'plotAll', plotAll); title('FlexMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
save2pdf(fullfile('Simulation_Results',sprintf('EMD_%s_results_lambda=%0.3e_lambdaM=%0.3e_lambdaR=%0.3e.pdf', data_type, lambda, lambda_Ms(Mi), lambda_R)))
figure; plot_MR(Ms{Mi,j}, Rs{Mi,j})

%% Validate EMD
Mi = 5;

W_hat = W_hats{Mi,j};
H_hat = H_hats{Mi,j};
Xhat = helper.reconstruct(W_hat, H_hat);
[emd_W, emd_H, ids] = helper.similarity_WH_EMD(W, Htrain, W_hat, H_hat);
% [d, M, R, out] = compute_EMD(Xtrain,Xhat,opts, 'continuationOptions', continue_opts, 'lambdaR', 1e2);
%% Plot EMDs relative to parameters
emd_Hs = zeros(nMs, nSim);
emd_Ws = zeros(nMs, nSim);
L1Ms = cellfun(@(x) norm(x(:),1), Ms);
L1Rs = cellfun(@(x) norm(x(:),1), Rs);
num_detected_all = cell2mat(num_detected);
num_significant_all = cell2mat(num_significant);
L1Ms_norm = L1Ms./max(L1Ms(:));
L1Rs_norm = L1Rs./max(L1Rs(:));
maxRecon = prctile(recon_costs(:),90);
minRecon = prctile(recon_costs(:),10);
maxReg = prctile(reg_costs(:),90);
minReg = prctile(reg_costs(:),10);
recon_costs_norm = (recon_costs-minRecon)./(maxRecon-minRecon);
reg_costs_norm = (reg_costs-minReg)./(maxReg-minReg);

for Mi = 1:nMs
    for n=1:nSim
        emd_Ws(Mi,n) = mean(emds_W{Mi,n});
        emd_Hs(Mi,n) =  mean(emds_H{Mi,n});
    end
end
emd_Hs_norm = emd_Hs./max(mean(emd_Hs,2));
emd_Ws_norm = emd_Ws./max(mean(emd_Ws,2));
emd_Hs_norm(isnan(emd_Hs_norm)) = 1;
emd_Ws_norm(isnan(emd_Ws_norm)) = 1;

c1 = [.1 .65 .35];
c2 = [.7 .2 .8];
figure;
ax1 = subplot('Position', [0.1 0.75 0.8 0.2]);
hold on
e1 = errorbar(lambda_Ms, median(L1Ms_norm,2), ...
    median(L1Ms_norm,2)-prctile(L1Ms_norm,25,2), ...
    prctile(L1Ms_norm,75,2)-median(L1Ms_norm,2), ...
    '-', 'Marker', '.', 'MarkerSize', 12, 'Color', c1);
e2 = errorbar(lambda_Ms, median(L1Rs_norm,2), ...
    median(L1Rs_norm,2)-prctile(L1Rs_norm,25,2), ...
    prctile(L1Rs_norm,75,2)-median(L1Rs_norm,2), ...
    '-', 'Marker', '.', 'MarkerSize', 12, 'Color', c2);
% yyaxis left
% swarmchart(lambda_Ms, L1Ms, 12, c1, 'filled')
% p1 = plot(lambda_Ms, mean(L1Ms,2), 'Color', c1);
% set(gca, 'XScale', 'log', 'YScale', 'log', 'YColor', c1)
% 
% yyaxis right
% hold on
% swarmchart(lambda_Ms, L1Rs, 12, c2, 'filled')
% p2 = plot(lambda_Ms, mean(L1Rs,2), 'Color', c2);
% set(gca, 'XScale', 'log', 'YColor', c2)
% 
% set(ax1, 'XLabel', [], 'XTicklabel', [])
% legend([p1, p2], {'L1M', 'L1R'}, 'Location','best')
set(ax1, 'XScale', 'log', 'XLabel', [], 'XTicklabel', [])
legend([e1, e2], {'L1M', 'L1R'}, 'Location','best')

ax2 = subplot('Position', [0.1 0.55 0.8 0.15]);
hold on
e1 = errorbar(lambda_Ms, median(reg_costs_norm,2), ...
    median(reg_costs_norm,2)-prctile(reg_costs_norm,25,2), ...
    prctile(reg_costs_norm,75,2)-median(reg_costs_norm,2), ...
    '-b.', 'MarkerSize', 12);
e2 = errorbar(lambda_Ms, median(recon_costs_norm,2), ...
    median(recon_costs_norm,2)-prctile(recon_costs_norm,25,2), ...
    prctile(recon_costs_norm,75,2)-median(recon_costs_norm,2), ...
    '-r.', 'MarkerSize', 12);
% swarmchart(lambda_Ms, reg_costs_norm, 12, 'b', 'filled')
% p1 = plot(lambda_Ms, mean(reg_costs_norm,2), 'Color', 'b');
% swarmchart(lambda_Ms, recon_costs_norm, 12, 'r', 'filled')
% p2 = plot(lambda_Ms, mean(recon_costs_norm,2), 'Color', 'r');
% set(ax2, 'XScale', 'log', 'XLabel', [], 'XTicklabel', [])
% legend([p1, p2], {'Reg Cost', 'Recon cost'}, 'Location','best')
set(ax2, 'XScale', 'log', 'XLabel', [], 'XTicklabel', [])
legend([e1, e2], {'Reg Cost', 'Recon cost'}, 'Location','best')

ax3 = subplot('Position', [0.1 0.35 0.8 0.15]);
hold on
% swarmchart(lambda_Ms, emd_Ws_norm, 12, 'k', 'filled')
% plot(lambda_Ms, mean(emd_Ws_norm,2), 'Color', 'k')
errorbar(lambda_Ms, median(emd_Ws_norm,2), ...
    median(emd_Ws_norm,2)-prctile(emd_Ws_norm,25,2), ...
    prctile(emd_Ws_norm,75,2)-median(emd_Ws_norm,2), ...
    '-k.', 'MarkerSize', 12);
set(ax3, 'XScale', 'log', 'XLabel', [], 'XTicklabel', [])
ylabel(ax3, 'EMD of W')

ax4 = subplot('Position', [0.1 0.1 0.8 0.2]);
plot(lambda_Ms, mean(num_significant_all,2), 'k')
set(ax4, 'XScale', 'log')
ylabel(ax4, 'Number Significant')
xlabel(ax4, 'lambda_M')
linkaxes([ax1, ax2, ax3, ax4], 'x')
set(gcf, 'Position', [100,100,400,800])

save2pdf(fullfile('Simulation_Results', sprintf('Choose_lambdaM_%s_plot.pdf', data_type)))

% % Average across simulations
% K = 3;
% Z = mean(emd_Ws, 3);   % (nLambdas x nMs)
% Z((mean_num_detected<K)) = nan;
% 
% % Make meshgrid of Lambda and M
% [LambdaGrid, MGrid] = meshgrid(lambdas, lambda_Ms);

% % Surface plot (log scaling handled by axes)
% figure;
% surf(LambdaGrid, MGrid, Z')   % transpose so dims match (M along Y)
% shading interp
% colorbar
% xlabel('\lambda')
% ylabel('M')
% zlabel('mean\_emd\_Hs')
% title('Surface plot of mean(mean\_emd\_Hs,3)')
% 
% % Set both x and y axes to log scale
% set(gca, 'XScale', 'log', 'YScale', 'log')
% view(45,30)  % adjust view angle
% grid on
% 
% %% 2-D plots of EMD
% % Average across simulations
% Zw = mean(emd_Ws, 3);   % (nLambdas x nMs)
% Zh = mean(emd_Hs, 3);   % (nLambdas x nMs)
% Zw((mean_num_detected<K)) = nan;
% Zh((mean_num_detected<K)) = nan;
% 
% graymap = gray(256); % Create 256 levels of gray
% flipped_graymap = flipud(graymap);
% 
% figure; imagesc(Zw, [min(Zw(:)), max(Zw(:))])
% colormap(flipped_graymap)
% colorbar
% title('mean EMD W')
% xticks(1:nMs);
% xticklabels(arrayfun(@(x) sprintf('%.0e', x), lambda_Ms, 'UniformOutput', false));
% yticks(1:nLambdas);
% yticklabels(arrayfun(@(x) sprintf('%.0e', x), lambdas, 'UniformOutput', false));
% xlabel('lambda_M')
% ylabel('lambda')
% set(gcf, 'Position', [100, 100, 600, 800])
% save2pdf(['EMD_W_choose_params_simulation_', data_type])
% 
% 
% figure; imagesc(Zh, [min(Zh(:)), max(Zh(:))])
% colormap(flipped_graymap)
% colorbar
% title('mean EMD H')
% xticks(1:nMs);
% xticklabels(arrayfun(@(x) sprintf('%.0e', x), lambda_Ms, 'UniformOutput', false));
% yticks(1:nLambdas);
% yticklabels(arrayfun(@(x) sprintf('%.0e', x), lambdas, 'UniformOutput', false));
% xlabel('lambda_M')
% ylabel('lambda')
% set(gcf, 'Position', [100, 100, 600, 800])
% save2pdf(['EMD_H_choose_params_simulation_', data_type])
