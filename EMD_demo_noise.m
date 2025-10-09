%% Demo script: FlexMF on warped/jittered data with noise
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
noise = .001; % probability of added noise in each bin
jitter = 5*ones(number_of_seqences,1); % Jitter std
participation = 1.*ones(number_of_seqences,1); % Participation parameter = 100%
warp = 2; % stretch should be less than Dt
gap = 100;
neg = 0;
seed = 1;
[X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',0, 'seed', seed, 'len_burst', 1, 'dynamic', 0);
[Xnoise, Wnoise, Hnoise, ~] = generate_data(T,Nneurons,Dt, 'noise',noise, 'seed', seed, 'len_burst', 1, 'dynamic', 0);
[Xwarp, Wwarp, Hwarp, ~] = generate_data(T,Nneurons,Dt, 'noise',noise, 'warp', warp, 'seed', seed, 'len_burst', 1, 'dynamic', 0);
[Xjit, Wjit, Hjit, ~] = generate_data(T,Nneurons,Dt, 'noise',noise, 'jitter', jitter, 'seed', seed, 'len_burst', 1, 'dynamic', 0);
[Xwarp_noise, Wwarp_noise, Hwarp_noise, ~] = generate_data(T,Nneurons,Dt, 'noise',noise, 'warp', warp, 'seed', seed, 'len_burst', 1, 'dynamic', 0);
L = size(W,3);

plotAll = 1;
figure; SimpleWHPlot(W,H,'Data',X, 'plotAll', plotAll); title('generated data raw','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot(W,H,'Data',Xnoise, 'plotAll', plotAll); title('generated data noise','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot(Wwarp,Hwarp,'Data',Xwarp,'plotAll', plotAll); title('generated data warping','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot(Wwarp_noise,Hwarp_noise,'Data',Xwarp_noise,'plotAll', plotAll); title('generated data warping w noise','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
save2pdf('EMD_simulated_data_warp_noise.pdf')
figure; SimpleWHPlot(Wjit,Hjit,'Data',Xjit,'plotAll',plotAll); title('generated data jittering','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

%% Normalize data
K = 3;
% frob_norm = norm(Xwarp(:));
% Xwarp = Xwarp/frob_norm*K;
% Wwarp = Wwarp/frob_norm*K;
% frob_norm = norm(Xnoise(:));
% Xnoise = Xnoise/frob_norm*K;
% Wnoise = Wnoise/frob_norm*K;

%% Procedure for choosing lambda in SeqNMF
nLambdas = 20; % increase if you're patient
lambdas = sort(logspace(-1,-5,nLambdas), 'ascend'); 
loadings = [];
regularization = []; 
cost = []; 
for li = 1:length(lambdas)
    [What_SeqNMF, Hhat_SeqNMF, ~,~,loadings(li,:),power]= seqNMF(Xwarp,'K',K,'L',L,...
        'lambda', lambdas(li), 'maxiter', 100, 'showPlot', 0); 
    [cost(li),regularization(li),~] = helper.get_seqNMF_cost(Xwarp,What_SeqNMF,Hhat_SeqNMF);
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

save2pdf('EMD_Simulated_warp_noise_data_SeqNMF.pdf', gcf)

%% Run FlexMF with EMD
lambda = 1e-5;
lambda_M = 1;
lambda_R = 1;
tic
figure;
% [What_FlexMF, Hhat_FlexMF, cost, errors_FlexMF, loadings, power, M, R] = FlexMF(Xnoise, 'K', K, 'L', L, ...
%     'EMD',1, 'lambda', lambda, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'maxiter', 50, 'tolerance', 1e-4);
[What_FlexMF, Hhat_FlexMF, cost, errors_FlexMF, loadings, power, M, R] = FlexMF(Xwarp_noise, 'K', K, 'L', L, ...
    'EMD',1, 'lambda', lambda, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'maxiter', 50, 'tolerance', 1e-4);
toc

%% Look at factors
plotAll = 1;
figure; SimpleWHPlot(What_FlexMF, Hhat_FlexMF, 'plotAll', plotAll); title('FlexMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot(What_FlexMF, Hhat_FlexMF, 'Data', Xwarp, 'plotAll', plotAll); title('FlexMF factors, with raw data')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
% save2pdf('EMD_Simulated_noise_data_FlexMF.pdf', gcf)
% save2pdf('EMD_Simulated_warp_noise_data_FlexMF.pdf', gcf)

%% Plot M, R
figure;
plot_MR(M,R)
% save2pdf('FlexMF_noise_demo.pdf')
% save2pdf('FlexMF_warp_noise_demo.pdf')