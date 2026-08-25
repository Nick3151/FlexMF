clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))

%% Generate synthetic data, with noise/probabilistic participation/warping/jittering
data_type = 'noise';
K = 3;
T = 4000; % length of data to generate
Nneurons = 5*ones(K,1); % number of neurons in each sequence
Dt = 5.*ones(K,1); % gap between each member of the sequence
neg = 0;
noise_levels = 0:.005:.02;
participation_levels = 1:-0.1:.5; 
jitter_levels = 0:.5:4.5;
warp_levels = 0:4;
gap = 100;

nSim = 10;
rng(1)
seeds = randperm(1000, nSim);

n = 1;
switch data_type
    case 'warp'
        [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',0, 'warp', warp_levels(end), 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
    case 'jitter'
        [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',0, 'jitter', jitter_levels(end).*ones(K,1), 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
    case 'participation'
        [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',0, 'participation', participation_levels(end).*ones(K,1), 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
    case 'noise'
        [X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise', noise_levels(2), 'seed', seeds(n), 'len_burst', 1, 'dynamic', 0);
end

L = size(W,3);
X_train = X(:,1:round(T/2));
X_test = X(:,1+round(T/2):end);
H_train = H(:,1:round(T/2));
H_test = H(:,1+round(T/2):end);

plotAll = 1;
figure; SimpleWHPlot(W,H_train,'Data',X_train, 'plotAll', plotAll); title('generated data raw','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

%% Split into training and test set, normalize data
frob_norm = norm(X_train(:));
X_train = X_train/frob_norm*K;
W = W/frob_norm*K;
frob_norm = norm(X_test(:));
X_test = X_test/frob_norm*K;

%% Procedure for choosing lambda in SeqNMF
nLambdas = 20; % increase if you're patient
lambdas = sort(logspace(-1,-5,nLambdas), 'ascend'); 
loadings = [];
regularization = []; 
cost = []; 
for li = 1:length(lambdas)
    [What_SeqNMF, Hhat_SeqNMF, ~,~,loadings(li,:),power]= seqNMF(X_train,'K',K,'L',L,...
        'lambda', lambdas(li), 'maxiter', 100, 'showPlot', 0); 
    [cost(li),regularization(li),~] = helper.get_seqNMF_cost(X_train,What_SeqNMF,Hhat_SeqNMF);
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
lambda = .05;
lambdaL1H = 0;
lambdaL1W = 0;
lambdaOrthoH = 0;

figure;
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
[What_SeqNMF, Hhat_SeqNMF, ~, errors_SeqNMF,loadings,power]= seqNMF(X_train,'K',K,'L',L,...
            'lambda', lambda, 'maxiter', 50, 'showPlot', 1); 

% plot, sorting neurons by latency within each factor
[max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(What_SeqNMF(:,:,:),1);
indSort = hybrid(:,3);

plotAll = 1;
figure; SimpleWHPlot(What_SeqNMF, Hhat_SeqNMF, 'plotAll', plotAll); title('SeqNMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot(What_SeqNMF, Hhat_SeqNMF, 'Data', X_train, 'plotAll', plotAll); title('SeqNMF factors, with raw data')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])