%% Demo script: Comparing SeqNMF and FlexMF on temporal warped/jittered data
clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))
rmpath(genpath(fullfile(root, 'seqNMF-master')));
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
seed = 1;
[X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',noise, 'seed', seed, 'len_burst', 1, 'dynamic', 0);
[Xwarp, Wwarp, Hwarp, ~] = generate_data(T,Nneurons,Dt, 'noise',noise, 'warp', warp, 'seed', seed, 'len_burst', 1, 'dynamic', 0);
[Xjit, Wjit, Hjit, ~] = generate_data(T,Nneurons,Dt, 'noise',noise, 'jitter', jitter, 'seed', seed, 'len_burst', 1, 'dynamic', 0);
L = size(W,3);

plotAll = 1;
figure; SimpleWHPlot(W,H,'Data',X, 'plotAll', plotAll); title('generated data raw','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
save2pdf('EMD_simulated_data_raw.pdf')
figure; SimpleWHPlot(Wwarp,Hwarp,'Data',Xwarp,'plotAll', plotAll); title('generated data warping','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
save2pdf('EMD_simulated_data_warp.pdf')
figure; SimpleWHPlot(Wjit,Hjit,'Data',Xjit,'plotAll',plotAll); title('generated data jittering','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
save2pdf('EMD_simulated_data_jitter.pdf')

%% Procedure for choosing lambda in SeqNMF
% Normalize data
K = 3;
frob_norm = norm(Xwarp(:));
Xwarp = Xwarp/frob_norm*K;
Wwarp = Wwarp/frob_norm*K;

nLambdas = 20; % increase if you're patient
lambdas = sort(logspace(-1,-5,nLambdas), 'ascend'); 
loadings = [];
regularization = []; 
cost = []; 
for li = 1:length(lambdas)
    [What_SeqNMF, Hhat_SeqNMF, ~,~,loadings(li,:),power]= seqNMF(Xwarp,'K',K,'L',L,...
        'lambda', lambdas(li), 'maxiter', 100, 'showPlot', 0); 
    [cost(li),regularization(li),~] = helper.get_seqNMF_cost(X,What_SeqNMF,Hhat_SeqNMF);
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
save2pdf('EMD_Simulate_choose_lambda_SeqNMF')

%% Run SeqNMF
lambda = .05;
lambdaL1H = 0;
lambdaL1W = 0;
lambdaOrthoH = 0;

figure;
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
[What_SeqNMF, Hhat_SeqNMF, ~, errors_SeqNMF,loadings,power]= seqNMF(Xwarp,'K',K,'L',L,...
            'lambda', lambda, 'maxiter', 50, 'showPlot', 1); 

% plot, sorting neurons by latency within each factor
[max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(What_SeqNMF(:,:,:),1);
indSort = hybrid(:,3);

%% Look at factors
plotAll = 1;
figure; SimpleWHPlot(What_SeqNMF, Hhat_SeqNMF, 'plotAll', plotAll); title('SeqNMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot(What_SeqNMF, Hhat_SeqNMF, 'Data', Xwarp, 'plotAll', plotAll); title('SeqNMF factors, with raw data')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

save2pdf('EMD_Simulated_warp_data_SeqNMF.pdf', gcf)

%% Run FlexMF with EMD
lambda_M = 1e-3;
lambda_R = 1e2;
tic
figure;
[What_FlexMF, Hhat_FlexMF, cost, errors_FlexMF, loadings, power, M, R] = FlexMF(Xwarp, 'K', K, 'L', L, ...
    'EMD',1, 'lambda', lambda, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'maxiter', 50);
toc

%% Look at factors
plotAll = 1;
figure; SimpleWHPlot(What_FlexMF, Hhat_FlexMF, 'plotAll', plotAll); title('FlexMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot(What_FlexMF, Hhat_FlexMF, 'Data', Xwarp, 'plotAll', plotAll); title('FlexMF factors, with raw data')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
save2pdf('EMD_Simulated_warp_data_FlexMF.pdf', gcf)

%% Plot M, R
figure;
plot_MR(M,R)
save2pdf('FlexMF_warp_demo.pdf')

%% Compare algorithms
[emds_W_SeqNMF, emds_H_SeqNMF, ids_SeqNMF] = helper.similarity_WH_EMD(Wwarp, Hwarp, What_SeqNMF, Hhat_SeqNMF);
[emds_W_FlexMF, emds_H_FlexMF, ids_FlexMF] = helper.similarity_WH_EMD(Wwarp, Hwarp, What_FlexMF, Hhat_FlexMF);
[coeffs_W_SeqNMF, coeffs_H_SeqNMF, ~] = helper.similarity_WH(Wwarp, Hwarp, What_SeqNMF, Hhat_SeqNMF);
tic
[coeffs_W_FlexMF, coeffs_H_FlexMF, ~] = helper.similarity_WH(Wwarp, Hwarp, What_FlexMF, Hhat_FlexMF);
toc

emds_W_all = zeros(2,K);
emds_W_all(1, ids_SeqNMF) = emds_W_SeqNMF;
emds_W_all(2, ids_FlexMF) = emds_W_FlexMF;

figure; bar(K:-1:1, emds_W_all);
legend({'SeqNMF', 'FlexMF'}, 'Location', 'north')
set(gca, 'FontSize', 14)
title('EMDs of W', 'FontSize', 16)
save2pdf('EMD_Simulated_warp_data_compare_W.pdf', gcf)

emds_H_all = zeros(2,K);
emds_H_all(1, ids_SeqNMF) = emds_H_SeqNMF;
emds_H_all(2, ids_FlexMF) = emds_H_FlexMF;

figure; bar(K:-1:1, emds_H_all);
legend({'SeqNMF', 'FlexMF'}, 'Location', 'north')
set(gca, 'FontSize', 14)
title('EMDs of H', 'FontSize', 16)
save2pdf('EMD_Simulated_warp_data_compare_H.pdf', gcf)