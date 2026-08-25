clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))

data_type = 'warpnoise';
lambda_R = 1;
lambda_M = 1e-2;

load(fullfile('Simulation_Results', sprintf('FlexMF_choose_lambda_%s_lambdaM=%0.3e.mat', data_type, lambda_M)))
[nLs, nSim] = size(H_hats);

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
% save2pdf(fullfile('Simulation_Results',sprintf('Synthetic_data_%s.pdf', data_type)))
%% Look at factors
Li = 9;
j = 5;
plotAll = 1;

figure; SimpleWHPlot_patch(W_hats{Li,j}, H_hats{Li,j}, 'plotAll', plotAll); title('FlexMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
save2pdf(fullfile('Simulation_Results',sprintf('EMD_%s_results_lambda=%0.3e_lambdaM=%0.3e_lambdaR=%0.3e.pdf', data_type, lambdas(Li), lambda_M, lambda_R)))

figure; plot_MR(Ms{Li,j}, Rs{Li,j})

%% Plot EMDs relative to parameters
emd_Hs = zeros(nLs, nSim);
emd_Ws = zeros(nLs, nSim);
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

for Mi = 1:nLs
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
e1 = errorbar(lambdas, median(L1Ms_norm,2), ...
    median(L1Ms_norm,2)-prctile(L1Ms_norm,25,2), ...
    prctile(L1Ms_norm,75,2)-median(L1Ms_norm,2), ...
    '-', 'Marker', '.', 'MarkerSize', 12, 'Color', c1);
e2 = errorbar(lambdas, median(L1Rs_norm,2), ...
    median(L1Rs_norm,2)-prctile(L1Rs_norm,25,2), ...
    prctile(L1Rs_norm,75,2)-median(L1Rs_norm,2), ...
    '-', 'Marker', '.', 'MarkerSize', 12, 'Color', c2);
% yyaxis left
% swarmchart(lambdas, L1Ms, 12, c1, 'filled')
% p1 = plot(lambdas, mean(L1Ms,2), 'Color', c1);
% set(gca, 'XScale', 'log', 'YScale', 'log', 'YColor', c1)
% yyaxis right
% swarmchart(lambdas, L1Rs, 12, c2, 'filled')
% p2 = plot(lambdas, mean(L1Rs,2), 'Color', c2);
% set(gca, 'XScale', 'log', 'YColor', c2)

set(ax1, 'XScale', 'log', 'XLabel', [], 'XTicklabel', [])
legend([e1, e2], {'L1M', 'L1R'}, 'Location','best')

ax2 = subplot('Position', [0.1 0.55 0.8 0.15]);
hold on
e1 = errorbar(lambdas, median(reg_costs_norm,2), ...
    median(reg_costs_norm,2)-prctile(reg_costs_norm,25,2), ...
    prctile(reg_costs_norm,75,2)-median(reg_costs_norm,2), ...
    '-b.', 'MarkerSize', 12);
e2 = errorbar(lambdas, median(recon_costs_norm,2), ...
    median(recon_costs_norm,2)-prctile(recon_costs_norm,25,2), ...
    prctile(recon_costs_norm,75,2)-median(recon_costs_norm,2), ...
    '-r.', 'MarkerSize', 12);
% swarmchart(lambdas, reg_costs_norm, 12, 'b', 'filled')
% p1 = plot(lambdas, mean(reg_costs_norm,2), 'Color', 'b');
% swarmchart(lambdas, recon_costs_norm, 12, 'r', 'filled')
% p2 = plot(lambdas, mean(recon_costs_norm,2), 'Color', 'r');
set(ax2, 'XScale', 'log', 'XLabel', [], 'XTicklabel', [])
legend([e1, e2], {'Reg Cost', 'Recon cost'}, 'Location','best')

ax3 = subplot('Position', [0.1 0.35 0.8 0.15]);
hold on
% swarmchart(lambdas, emd_Ws_norm, 12, 'k', 'filled')
% plot(lambdas, mean(emd_Ws_norm,2), 'Color', 'k')
errorbar(lambdas, median(emd_Ws_norm,2), ...
    median(emd_Ws_norm,2)-prctile(emd_Ws_norm,25,2), ...
    prctile(emd_Ws_norm,75,2)-median(emd_Ws_norm,2), ...
    '-k.', 'MarkerSize', 12);
set(ax3, 'XScale', 'log', 'XLabel', [], 'XTicklabel', [])
ylabel(ax3, 'EMD of W')

ax4 = subplot('Position', [0.1 0.1 0.8 0.2]);
plot(lambdas, mean(num_significant_all,2), 'k')
set(ax4, 'XScale', 'log')
ylabel(ax4, 'Number Significant')
xlabel(ax4, 'lambda')
linkaxes([ax1, ax2, ax3, ax4], 'x')
set(gcf, 'Position', [100,100,400,800])

save2pdf(fullfile('Simulation_Results', sprintf('Choose_lambda_%s_plot.pdf', data_type)))
