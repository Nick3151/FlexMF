%% Demo script showing detecing sequences with EMD as cost function
clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))
% rmpath(genpath(fullfile(root, 'seqNMF-master')));
addpath(genpath(fullfile(root, 'FlexMF')));
%% Generate some synthetic data with temporal jittering or time warping
number_of_seqences = 3;
T = 800; % length of data to generate
Nneurons = 10*ones(number_of_seqences,1); % number of neurons in each sequence
Dt = 3.*ones(number_of_seqences,1); % gap between each member of the sequence
noise = 0; % probability of added noise in each bin
jitter = 5*ones(number_of_seqences,1); % Jitter std
participation = 1.*ones(number_of_seqences,1); % Participation parameter = 100%
warp = 2; % stretch should be less than Dt
gap = 100;
neg = 0;
bin = 0;
seed = 1;
[X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',noise, 'seed', seed);
[Xwarp, Wwarp, Hwarp, ~] = generate_data(T,Nneurons,Dt, 'noise',noise, 'warp', warp, 'seed', seed);
[Xjit, Wjit, Hjit, ~] = generate_data(T,Nneurons,Dt, 'noise',noise, 'jitter', jitter, 'seed', seed);
L = size(W,3);
% range = round(L/2)-25:round(L/2)+35;

plotAll = 1;
figure; SimpleWHPlot_patch(W,H,'Data',X, 'plotAll', plotAll); title('generated data raw','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
save2pdf('EMD_simulated_data_raw.pdf')
figure; SimpleWHPlot_patch(Wwarp,Hwarp,'Data',Xwarp,'plotAll', plotAll); title('generated data warping','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
save2pdf('EMD_simulated_data_warp.pdf')
figure; SimpleWHPlot_patch(Wjit,Hjit,'Data',Xjit,'plotAll',plotAll); title('generated data jittering','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
save2pdf('EMD_simulated_data_jitter.pdf')


%% EMD for sequence detection, raw data
% Normalize data
K = 10;
frob_norm = norm(X(:));
X = X/frob_norm*K;

figure;
L = 60;
lambdaL1H = 1e-3;
lambda = 1e-4;
lambda_M = 1e-4;
lambda_R = 1e2;
[What, Hhat, cost, errors, loadings, power, M, R] = FlexMF(X, 'K', K, 'L', L, ...
    'EMD',1, 'lambda', lambda, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'lambdaL1H', lambdaL1H, 'tol', 1e-4, 'maxiter', 20, 'Reweight', 1);

% figure;
% SimpleWHPlot_patch(What, Hhat, 'Data', X, 'plotAll', 1)
% figure;
% SimpleWHPlot_patch(What, Hhat, 'plotAll', 1)
% save2pdf(sprintf('EMD_raw_lambda=%0.2e_lambdaM=%0.2e_lambdaR=%0.2e_lambdaL1H=%0.2e_results.pdf', lambda, lambda_M, lambda_R, lambdaL1H))

figure;
SimpleWHPlot_patch(What, Hhat, 'Data', X, 'plotAll', 1)
figure;
SimpleWHPlot_patch(What, Hhat, 'plotAll', 1)
% save2pdf(sprintf('EMD_raw_reweighted_lambda=%0.2e_lambdaM=%0.2e_lambdaR=%0.2e_lambdaL1H=%0.2e_results.pdf', lambda, lambda_M, lambda_R, lambdaL1H))
save2pdf(sprintf('EMD_raw_smooth_reweighted_lambda=%0.2e_lambdaM=%0.2e_lambdaR=%0.2e_lambdaL1H=%0.2e_results.pdf', lambda, lambda_M, lambda_R, lambdaL1H))
%% Plot L1W and L1H as a funciton of iterations
figure; 
yyaxis left
plot(cost(2:end))
yyaxis right
hold on
plot(errors(2:end,3))
plot(errors(2:end,4))
xlabel('Iteration #')
legend('EMD', 'L1W', 'L1H')
xline(10, '--', 'Reweighted L1')
save2pdf(sprintf('EMD_raw_reweighted_lambda=%0.2e_lambdaM=%0.2e_lambdaR=%0.2e_lambdaL1H=%0.2e_costs.pdf', lambda, lambda_M, lambda_R, lambdaL1H))
%% EMD for sequence detection, warping data
% Normalize data
K = 3;
frob_norm = norm(Xwarp(:));
Xwarp = Xwarp/frob_norm*K;

figure;
L = 60;
lambdaL1H = 1e-3;
lambda = 1e-4;
lambda_M = 1e-4;
lambda_R = 1e2;
lambda_TV = 1e-4;

% Fix W and fit H, params may be different...
% [What, Hhat, cost, errors, loadings, power, M, R] = FlexMF(Xwarp, 'K', K, 'L', L, ...
%     'W_init', Wwarp, 'H_init', Hwarp, 'W_fixed', 1, 'EMD',1, 'lambda', 1e-4, 'lambda_R', 1, 'lambda_M', 1e-6, 'lambdaL1H', lambdaL1H, 'maxiter', 10);

[What, Hhat, cost, errors, loadings, power, M, R] = FlexMF(Xwarp, 'K', K, 'L', L, ...
    'EMD',1, 'lambda', lambda, 'lambda_R', lambda_R, 'lambda_M', lambda_M, ...
    'lambda_TV', lambda_TV, 'lambdaL1H', lambdaL1H, 'maxiter', 50, 'Reweight', 1);

% figure;
% SimpleWHPlot(What, Hhat, 'Data', Xwarp, 'plotAll', 1, 'neg', 1)
% figure;
% SimpleWHPlot(What, Hhat, 'plotAll', 1, 'neg', 1)
figure;
SimpleWHPlot_patch(What, Hhat, 'Data', Xwarp, 'plotAll', 1)
figure;
SimpleWHPlot_patch(What, Hhat, 'plotAll', 1)
% save2pdf(sprintf('EMD_warp_lambda=%0.2e_lambdaM=%0.2e_lambdaR=%0.2e_lambdaL1H=%0.2e_results.pdf', lambda, lambda_M, lambda_R, lambdaL1H))
% save2pdf(sprintf('EMD_warp_reweighted_lambda=%0.2e_lambdaM=%0.2e_lambdaR=%0.2e_lambdaL1H=%0.2e_results.pdf', lambda, lambda_M, lambda_R, lambdaL1H))
save2pdf(sprintf('EMD_warp_reweighted_TVnorm_lambda=%0.2e_lambdaM=%0.2e_lambdaR=%0.2e_lambdaL1H=%0.2e_results.pdf', lambda, lambda_M, lambda_R, lambdaL1H))

tmp = helper.reconstruct(What,Hhat)-Xwarp;
% sum(tmp(:).^2/2)save2pdf(sprintf('EMD_warp_reweighted_TVnorm_lambda=%0.2e_lambdaM=%0.2e_lambdaR=%0.2e_lambdaL1H=%0.2e_results.pdf', lambda, lambda_M, lambda_R, lambdaL1H))

% Compute similarity to ground truth
[coeffs_W, coeffs_H, ids] = helper.similarity_WH(W, H, What, Hhat);

%% Plot EMD and other costs as a funciton of iterations
figure; 
yyaxis left
plot(cost(2:end))
hold on
plot(errors(2:end,2))
yyaxis right
hold on
plot(errors(2:end,3))
plot(errors(2:end,4))
xlabel('Iteration #')
legend('EMD', 'Regularization', 'L1W', 'L1H')
title('EMD')
% save2pdf(sprintf('EMD_warp_lambda=%0.2e_lambdaM=%0.2e_lambdaR=%0.2e_lambdaL1H=%0.2e_costs.pdf', lambda, lambda_M, lambda_R, lambdaL1H))

%% Plot M and R 

figure;
ax_res = subplot('Position', [0.05, 0.55, 0.8, 0.4]);
imagesc(R,[-max(abs(R(:))),max(abs(R(:)))])
title('R', 'FontSize', 16)
set(ax_res, 'XTickLabel', [], 'YTickLabel', []);
colorbar('Position', [0.9 0.55 0.05 0.4], 'FontSize', 14);

ax_flux = subplot('Position', [0.05, 0.05, 0.8, 0.4]);
imagesc(M, [-max(abs(M(:))),max(abs(M(:)))])
title('M', 'FontSize', 16)
set(ax_flux, 'XTickLabel', [], 'YTickLabel', []);
colorbar('Position', [0.9 0.05 0.05 0.4], 'FontSize', 14);
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
% save2pdf(sprintf('EMD_warp_lambda=%0.2e_lambdaM=%0.2e_lambdaR=%0.2e_lambdaL1H=%0.2e_MR', lambda, lambda_M, lambda_R, lambdaL1H))

%% EMD between X and Xwarp
opts_default = tfocs_SCD;
opts = opts_default;
opts.continuation = 1;
opts.tol = 1e-6;
opts.stopCrit = 4;
opts.maxIts = 500;

[d, M, R, out] = compute_EMD(X, Xwarp, opts);

figure;
ax_res = subplot('Position', [0.05, 0.55, 0.8, 0.4]);
imagesc(R)
title('R', 'FontSize', 16)
set(ax_res, 'XTickLabel', [], 'YTickLabel', []);
colorbar('Position', [0.9 0.55 0.05 0.4], 'FontSize', 14);
ax_flux = subplot('Position', [0.05, 0.05, 0.8, 0.4]);
imagesc(M)
title('M', 'FontSize', 16)
set(ax_flux, 'XTickLabel', [], 'YTickLabel', []);
colorbar('Position', [0.9 0.05 0.05 0.4], 'FontSize', 14);
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
% save2pdf('EMD_warp_demo.pdf')
save2pdf('EMD_shift_demo.pdf')