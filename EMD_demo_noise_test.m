%% Test script: FlexMF on warped/jittered data with noise
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
T = 2000; % length of data to generate
Nneurons = 10*ones(number_of_seqences,1); % number of neurons in each sequence
Dt = 3.*ones(number_of_seqences,1); % gap between each member of the sequence
noise = .001; % probability of added noise in each bin
jitter = 5*ones(number_of_seqences,1); % Jitter std
participation = 1.*ones(number_of_seqences,1); % Participation parameter = 100%
warp = 2; % stretch should be less than Dt
gap = 100;
neg = 0;
seed = 1;

[Xwarp, Wwarp, Hwarp, ~] = generate_data(T,Nneurons,Dt, 'noise',noise, 'warp', warp, 'seed', seed, 'len_burst', 1, 'dynamic', 0);

L = size(Wwarp,3);

plotAll = 1;
figure; SimpleWHPlot(Wwarp,Hwarp,'Data',Xwarp,'plotAll', plotAll); title('generated data warping','Fontsize',16)

%% Normalize data
K = 3;

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

%% Run SeqNMF
lambda_SeqNMF = .05;
lambdaL1H = 0;
lambdaL1W = 0;
lambdaOrthoH = 0;

figure;
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
[What_SeqNMF, Hhat_SeqNMF, ~, errors_SeqNMF,loadings,power]= seqNMF(Xwarp,'K',K,'L',L,...
            'lambda', lambda_SeqNMF, 'maxiter', 50, 'showPlot', 1); 

% plot, sorting neurons by latency within each factor
[max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(What_SeqNMF(:,:,:),1);
indSort = hybrid(:,3);

%% Look at factors
plotAll = 1;
figure; SimpleWHPlot_patch(What_SeqNMF, Hhat_SeqNMF, 'plotAll', plotAll); title('SeqNMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot_patch(What_SeqNMF, Hhat_SeqNMF, 'Data', Xwarp, 'plotAll', plotAll); title('SeqNMF factors, with raw data')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

%% Consistency
lambdaL1H = 0;
lambdaL1W = 0;
lambda_FlexMF = 1e-4;
lambda_M = .1;
lambda_R = 1;

nRuns = 10;
Whats = cell(nRuns,1);
Hhats = cell(nRuns,1);
errors = zeros(nRuns, 4);
Ms = cell(nRuns,1);
Rs = cell(nRuns,1);

parfor n=1:nRuns
%     lambda_FlexMF = lambdas(n);
    disp(lambda_FlexMF)
    
%     [Whats{n}, Hhats{n}, cost, error, loadings, power, Ms{n}, Rs{n}] = FlexMF(Xwarp, 'K', K, 'L', L, ...
%         'EMD',1, 'lambda', lambda_FlexMF, ...
%         'lambdaL1H', lambdaL1H, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'maxiter', 1, 'tolerance', 1e-3, ...
%         'W_init', Wwarp, 'W_fixed', 1, 'showPlot', 0, 'verbal', 0);
    figure;     
    [Whats{n}, Hhats{n}, cost, error, loadings, power, Ms{n}, Rs{n}] = FlexMF(Xwarp, 'K', K, 'L', L, ...
            'EMD',1, 'lambda', lambda_FlexMF, 'lambdaL1W', lambdaL1W, ...
            'lambdaL1H', lambdaL1H, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'maxiter', 50, 'tolerance', 1e-6, ...
            'W_init', What_SeqNMF, 'H_init', Hhat_SeqNMF, 'showPlot', 0, 'verbal', 0);
        
    errors(n,:) = error(end,:);
end

%% Look at factors
plotAll = 1;
n = 7;
figure; SimpleWHPlot_patch(Whats{n}, Hhats{n}, 'plotAll', plotAll); title('FlexMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot_patch(Whats{n}, Hhats{n}, 'Data', Xwarp, 'plotAll', plotAll); title('FlexMF factors, with raw data')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
% save2pdf(sprintf('EMD_Simulated_warp_noise_data_FlexMF_T=%d_lambda=%1.1e_lambdaM=%1.1e_lambdaR=%1.1e.pdf', T, lambda, lambda_M, lambda_R), gcf)
figure;
plot_MR(Ms{n},Rs{n})
% save2pdf(sprintf('FlexMF_warp_noise_demo_MR_T=%d_lambda=%1.1e_lambdaM=%1.1e_lambdaR=%1.1e.pdf', T, lambda, lambda_M, lambda_R))

L1Ms = cellfun(@(x) norm(x(:),1), Ms);
L1Rs = cellfun(@(x) norm(x(:),1), Rs);
L1Ws = cellfun(@(x) norm(x(:),1), Whats);
L1Hs = cellfun(@(x) norm(x(:),1), Hhats);
%% Errors vs lambdaL1H
figure;
plot(lambdas, errors(:,1), lambdas, L1Hs)
legend('Recon', 'L1H')
xlabel('Lambda')
export_vector_pdf('Error vs lambda');
