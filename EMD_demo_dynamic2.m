%% Demo script showing FlexMF performance on data with calcium temporal dynamics
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
T = 4000; % length of data to generate
Nneurons = 10*ones(number_of_seqences,1); % number of neurons in each sequence
Dt = 3.*ones(number_of_seqences,1); % gap between each member of the sequence
noise = 0.001; % probability of added noise in each bin
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
save2pdf('EMD_simulated_data_dynamic_raw.pdf')
figure; SimpleWHPlot_patch(Wwarp,Hwarp,'Data',Xwarp,'plotAll', plotAll); title('generated data warping','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
save2pdf('EMD_simulated_data_dynamic_warp.pdf')
figure; SimpleWHPlot_patch(Wjit,Hjit,'Data',Xjit,'plotAll',plotAll); title('generated data jittering','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
save2pdf('EMD_simulated_data_dynamic_jitter.pdf')

%% Split into training and test set
Xtrain = Xwarp(:,1:round(T/2));
Xtest = Xwarp(:,1+round(T/2):end);

figure; SimpleWHPlot(Wwarp,Hwarp,'Data',Xwarp,'plotAll', 1, 'onsets', round(T/2)); 
title('generated data warping','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.5])

%% Procedure for choosing lambda in SeqNMF
K = 3;
L = 50;
nLambdas = 20; % increase if you're patient
lambdas = logspace(-5,-1,nLambdas); 
loadings = [];
regularization = []; 
cost = []; 
for li = 1:length(lambdas)
    [What_SeqNMF, Hhat_SeqNMF, ~,~,loadings(li,:),power]= seqNMF(Xtrain,'K',K,'L',L,...
        'lambda', lambdas(li), 'maxiter', 100, 'showPlot', 0); 
    [cost(li),regularization(li),~] = helper.get_seqNMF_cost(Xtrain,What_SeqNMF,Hhat_SeqNMF);
    display(['Testing lambda ' num2str(li) '/' num2str(length(lambdas))])
end

%% plot costs as a function of lambda
windowSize = 3; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
Rs = filtfilt(b,a,regularization); 
minRs = prctile(regularization,10); maxRs= prctile(regularization,90);
Rs = (Rs-minRs)/(maxRs-minRs); 
R = (regularization-minRs)/(maxRs-minRs); 
Cs = filtfilt(b,a,cost); 
minCs =  prctile(cost,10); maxCs =  prctile(cost,90); 
Cs = (Cs -minCs)/(maxCs-minCs); 
C = (cost -minCs)/(maxCs-minCs); 

figure; hold on
plot(lambdas,Rs, 'b')
plot(lambdas,Cs,'r')
scatter(lambdas, R, 'b', 'markerfacecolor', 'flat');
scatter(lambdas, C, 'r', 'markerfacecolor', 'flat');
xlabel('Lambda'); ylabel('Cost (au)')
set(legend('Correlation cost', 'Reconstruction cost'), 'Box', 'on')
set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
% save2pdf('Simulate_warp_noise_dynamics_choose_lambda_SeqNMF')

%% Run SeqNMF on training data
K = 3;
L = 50;
lambda = .01;
lambdaL1H = 0;
lambdaL1W = 0;
lambdaOrthoH = 0;

figure;
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
[What_SeqNMF, Hhat_SeqNMF, ~, errors_SeqNMF,loadings,power]= seqNMF(Xtrain,'K',K,'L',L,...
            'lambda', lambda, 'maxiter', 50, 'showPlot', 1); 

%% Look at factors
plotAll = 1;
figure; SimpleWHPlot_patch(What_SeqNMF, Hhat_SeqNMF, 'plotAll', plotAll); title('SeqNMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
% figure; SimpleWHPlot_patch(What_SeqNMF, Hhat_SeqNMF, 'Data', Xwarp, 'plotAll', plotAll); title('SeqNMF factors, with raw data')
figure; SimpleWHPlot_patch(What_SeqNMF, Hhat_SeqNMF, 'Data', Xtrain, 'plotAll', plotAll); title('SeqNMF factors, with raw data')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
save2pdf(sprintf('Simulated_warp_noise_dynamic_result_SeqNMF'))

%% EMD for sequence detection, with calcium dynamics
figure;
K = 3;
L = 50;
lambdaL1H = 1;
lambda = 1e-2;
lambda_M = 1e-2;
lambda_R = 1;
[What, Hhat, cost, errors, loadings, power, M, R] = FlexMF(Xtrain, 'K', K, 'L', L, ...
    'EMD',1, 'lambda', lambda, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'lambdaL1H', lambdaL1H, ...
    'W_init', What_SeqNMF, 'H_init', Hhat_SeqNMF, 'maxiter', 50, 'Reweight', 1);

% figure;
% SimpleWHPlot_patch(What, Hhat, 'Data', Xtrain, 'plotAll', 1)
% figure;
% SimpleWHPlot_patch(What, Hhat, 'plotAll', 1)
% save2pdf(sprintf('Simulated_warp_noise_dynamic_result_EMD_lambda=%0.2e_lambdaM=%0.2e_lambdaR=%0.2e_lambdaL1H=%0.2e_results.pdf', lambda, lambda_M, lambda_R, lambdaL1H))

figure;
SimpleWHPlot_patch(What, Hhat, 'Data', Xtrain, 'plotAll', 1, 'compare', true)
figure;
SimpleWHPlot_patch(What, Hhat, 'plotAll', 1)
save2pdf(sprintf('Simulated_warp_noise_dynamic_result_EMD_reweightL1H_lambda=%0.2e_lambdaM=%0.2e_lambdaR=%0.2e_lambdaL1H=%0.2e_results.pdf', lambda, lambda_M, lambda_R, lambdaL1H))

%% Plot M and R 
figure;
plot_MR(M,R)

% Ttrain = size(Xtrain,2);
% D = eye(Ttrain) - diag(ones(Ttrain-1,1),-1);
% D(Ttrain,Ttrain) = 0;
% constraint = M*D'-R-helper.reconstruct(What,Hhat)+Xtrain;
% figure;
% imagesc(constraint)
% colorbar

%% Fix What, rerun FlexMF on test data
figure;
[What, Hhat_test, cost_test, errors_test, ~, ~, M_test, R_test] = FlexMF(Xtest, 'K', K, 'L', L, 'W_fixed', 1, 'W_init', What,...
    'EMD',1, 'lambda', lambda, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'lambdaL1H', lambdaL1H, 'maxiter', 10, 'Reweight', 1);

figure; SimpleWHPlot_patch(What, Hhat_test, 'plotAll', 1); title('FlexMF test recon')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; plot_MR(M_test, R_test)

%% test significance
[pvals,is_significant,is_single] = test_significance_EMD(Xtest, What, M_test, 'plot', 0);

%% Plot L1W and L1H as a funciton of iterations
figure; 
yyaxis left
p1 = plot(errors(2:end,2));
yyaxis right
hold on
p2 = plot(errors(2:end,3));
p3 = plot(errors(2:end,4));
xline(10, '--', 'Reweighted L1')
xlabel('Iteration #')
legend([p1,p2,p3], {'Reg', 'L1W', 'L1H'})

save2pdf(sprintf('Simulated_warp_noise_dynamic_result_EMD_reweightL1H_lambda=%0.2e_lambdaM=%0.2e_lambdaR=%0.2e_lambdaL1H=%0.2e_costs.pdf', lambda, lambda_M, lambda_R, lambdaL1H))
%% EMD + reweight L1H + TV norm on W
K = 3;
figure;
L = 50;
lambdaL1H = 1;
lambda = 1e-2;
lambda_M = 1e-2;
lambda_R = 1;
lambda_TV = 1e-4;

% Fix W and fit H, params may be different...
% [What, Hhat, cost, errors, loadings, power, M, R] = FlexMF(Xtrain, 'K', K, 'L', L, ...
%     'W_init', Wwarp, 'H_init', Hwarp, 'W_fixed', 1, 'EMD',1, 'lambda', 1e-4, 'lambda_R', 1, 'lambda_M', 1e-6, 'lambdaL1H', lambdaL1H, 'maxiter', 10);

[What, Hhat, cost, errors, loadings, power, M, R] = FlexMF(Xtrain, 'K', K, 'L', L, ...
    'EMD',1, 'lambda', lambda, 'lambda_R', lambda_R, 'lambda_M', lambda_M, ...
    'lambda_TV', lambda_TV, 'lambdaL1H', lambdaL1H, 'maxiter', 50, 'Reweight', 1);

figure;
SimpleWHPlot_patch(What, Hhat, 'Data', Xtrain, 'plotAll', 1)
figure;
SimpleWHPlot_patch(What, Hhat, 'plotAll', 1)

save2pdf(sprintf('Simulated_warp_noise_dynamic_result_EMD_reweighted_lambdaTV=%0.2e_lambda=%0.2e_lambdaM=%0.2e_lambdaR=%0.2e_lambdaL1H=%0.2e_results.pdf', lambda_TV, lambda, lambda_M, lambda_R, lambdaL1H))

% Compute similarity to ground truth
[emds_W, emds_H, ids] = helper.similarity_WH_EMD(W, H, What, Hhat);
% sum(W,[1,3]);
% sum(H,2)

%% Plot EMD and other costs as a funciton of iterations
figure; 
yyaxis left
plot(errors(2:end,2))
yyaxis right
hold on
plot(errors(2:end,3))
plot(errors(2:end,4))
xlabel('Iteration #')
legend('Regularization', 'L1W', 'L1H')
title('EMD')
save2pdf(sprintf('Simulated_warp_noise_dynamic_result_EMD_reweighted_lambdaTV=%0.2e_lambda=%0.2e_lambdaM=%0.2e_lambdaR=%0.2e_lambdaL1H=%0.2e_costs.pdf', lambda_TV, lambda, lambda_M, lambda_R, lambdaL1H))

%% Plot M and R 
figure;
plot_MR(M,R)
