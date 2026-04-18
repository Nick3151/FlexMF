%% Robustness of FlexMF to initialization
clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))
rmpath(genpath(fullfile(root, 'seqNMF-master')));
addpath(genpath(fullfile(root, 'FlexMF')));

%% Generate some synthetic data with various noise
K = 3;
T = 2000; % length of data to generate
Nneurons = 5*ones(K,1); % number of neurons in each sequence
Dt = 3.*ones(K,1); % gap between each member of the sequence
noise = .001; % probability of added noise in each bin
jitter = 2*ones(K,1); % Jitter std
participation = .8.*ones(K,1); % Participation parameter = 100%
warp = 2; % stretch should be less than Dt
gap = 100;
neg = 0;
seed = 1;
[X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',0, 'seed', seed, 'len_burst', 1, 'dynamic', 0);
[Xpart, Wpart, Hpart, ~] = generate_data(T,Nneurons,Dt, 'noise', 0, 'participation', participation, 'seed', seed, 'len_burst', 1, 'dynamic', 0);
L = size(W,3);

plotAll = 1;
figure; SimpleWHPlot(W,H,'Data',X, 'plotAll', plotAll); title('generated data raw','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot(Wpart,Hpart,'Data',Xpart, 'plotAll', plotAll); title('generated data participation','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

%% Normalize data
Khat = 5;
frob_norm = norm(Xpart(:));
Xpart = Xpart/frob_norm*Khat;
Wpart = Wpart/frob_norm*Khat;

%% Procedure for choosing lambda in SeqNMF
nLambdas = 20; % increase if you're patient
lambdas = sort(logspace(-1,-5,nLambdas), 'ascend'); 
loadings = [];
regularization = []; 
cost = []; 
for li = 1:length(lambdas)
    [What_SeqNMF, Hhat_SeqNMF, ~,~,loadings(li,:),power]= seqNMF(Xpart,'K',Khat,'L',L,...
        'lambda', lambdas(li), 'maxiter', 100, 'showPlot', 0); 
    [cost(li),regularization(li),~] = helper.get_seqNMF_cost(Xpart,What_SeqNMF,Hhat_SeqNMF);
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

% save2pdf('Simulate_jitter_noise_choose_lambda_SeqNMF')

%% Run SeqNMF or FlexMF multiple times
lambdaL1H = 0;
lambdaL1W = 0;
lambdaOrthoH = 0;
lambda_M = .1;
lambda_R = 1;

% Running nSim times
nSim = 10;
Whats_SeqNMF = cell(nSim,1);
Hhats_SeqNMF = cell(nSim,1);
% First column: random initialization
% second column: initialize with seqNMF results
Whats_FlexMF = cell(nSim,2);
Hhats_FlexMF = cell(nSim,2);
Ms = cell(nSim,2);
Rs = cell(nSim,2);
objs = zeros(nSim,2);
emds_W_SeqNMF = cell(nSim,1);
emds_H_SeqNMF = cell(nSim,1);
ids_SeqNMF = cell(nSim,1);
emds_W_FlexMF = cell(nSim,2);
emds_H_FlexMF = cell(nSim,2);
ids_FlexMF = cell(nSim,2);

% figure;
% set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
for n=1:nSim
    fprintf('n=%d\n', n);
    tic
    lambda = .1;
    [What_SeqNMF, Hhat_SeqNMF]= seqNMF(Xpart,'K',Khat,'L',L,...
            'lambda', lambda, 'maxiter', 50, 'showPlot', 0); 
    Whats_SeqNMF{n} = What_SeqNMF;
    Hhats_SeqNMF{n} = Hhat_SeqNMF;
    [emds_W_SeqNMF{n}, emds_H_SeqNMF{n}, ids_SeqNMF{n}] = helper.similarity_WH_EMD(Wpart, Hpart, What_SeqNMF, Hhat_SeqNMF);
    toc
    tic
    lambda = .01;
    [What_FlexMF, Hhat_FlexMF, cost, errors_FlexMF, loadings, power, M, R] = FlexMF(Xpart, 'K', Khat, 'L', L, ...
        'EMD',1, 'lambda', lambda, ...
        'lambdaL1H', lambdaL1H, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'maxiter', 50, 'tolerance', 1e-3, ...
        'neg_prop', 0, 'Reweight', 1, 'showPlot', 0, 'verbal', 0);
    Whats_FlexMF{n,1} = What_FlexMF;
    Hhats_FlexMF{n,1} = Hhat_FlexMF;
    Ms{n,1} = M;
    Rs{n,1} = R;
    objs(n,1) = lambda*errors_FlexMF(end,2)+lambda_M*norm(M(:),1)+lambda_R*norm(R(:),1);
    [emds_W_FlexMF{n,1}, emds_H_FlexMF{n,1}, ids_FlexMF{n,1}] = helper.similarity_WH_EMD(Wpart, Hpart, What_FlexMF, Hhat_FlexMF);
    toc

    tic
    [What_FlexMF, Hhat_FlexMF, cost, errors_FlexMF, loadings, power, M, R] = FlexMF(Xpart, 'K', Khat, 'L', L, ...
        'EMD',1, 'lambda', lambda, ...
        'lambdaL1H', lambdaL1H, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'maxiter', 50, 'tolerance', 1e-3, ...
        'W_init', What_SeqNMF, 'H_init', Hhat_SeqNMF, 'showPlot', 0, 'verbal', 0);
    Whats_FlexMF{n,2} = What_FlexMF;
    Hhats_FlexMF{n,2} = Hhat_FlexMF;
    Ms{n,2} = M;
    Rs{n,2} = R;
    objs(n,2) = lambda*errors_FlexMF(end,2)+lambda_M*norm(M(:),1)+lambda_R*norm(R(:),1);
    [emds_W_FlexMF{n,2}, emds_H_FlexMF{n,2}, ids_FlexMF{n,2}] = helper.similarity_WH_EMD(Wpart, Hpart, What_FlexMF, Hhat_FlexMF);
    toc
end

% % plot, sorting neurons by latency within each factor
% [max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(What_SeqNMF(:,:,:),1);
% indSort = hybrid(:,3);

%% Look at factors seqNMF
n = 2;
What_SeqNMF = Whats_SeqNMF{n};
Hhat_SeqNMF = Hhats_SeqNMF{n};
plotAll = 1;
figure; SimpleWHPlot_patch(What_SeqNMF, Hhat_SeqNMF, 'plotAll', plotAll); title('SeqNMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

figure; SimpleWHPlot_patch(What_SeqNMF, Hhat_SeqNMF, 'Data', Xpart, 'plotAll', plotAll); title('SeqNMF factors, with raw data')
% figure; SimpleWHPlot_patch(What_SeqNMF, Hhat_SeqNMF, 'Data', Xjit, 'plotAll', plotAll); title('SeqNMF factors, with raw data')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

%% Run FlexMF with EMD
% lambda = .01;

%% Look at factors FlexMF
What_FlexMF = Whats_FlexMF{n,1};
Hhat_FlexMF = Hhats_FlexMF{n,1};
plotAll = 1;
figure; SimpleWHPlot_patch(What_FlexMF, Hhat_FlexMF, 'plotAll', plotAll); title('FlexMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

What_FlexMF = Whats_FlexMF{n,2};
Hhat_FlexMF = Hhats_FlexMF{n,2};
plotAll = 1;
figure; SimpleWHPlot_patch(What_FlexMF, Hhat_FlexMF, 'plotAll', plotAll); title('FlexMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

%% Plot M, R
M = Ms{n,1};
R = Rs{n,1};
figure;
plot_MR(M,R)

M = Ms{n,2};
R = Rs{n,2};
figure;
plot_MR(M,R)

%% Constraints for FlexMF
constraints = zeros(nSim,2);
for n=1:nSim
    What_FlexMF = Whats_FlexMF{n,1};
    Hhat_FlexMF = Hhats_FlexMF{n,1};
    M = Ms{n,1};
    R = Rs{n,1};
    Xcorr = helper.correct_warp(Xpart,M);
    constraint = Xcorr-R-helper.reconstruct(What_FlexMF,Hhat_FlexMF);
    constraints(n,1) = norm(constraint(:),1);

    What_FlexMF = Whats_FlexMF{n,2};
    Hhat_FlexMF = Hhats_FlexMF{n,2};
    M = Ms{n,2};
    R = Rs{n,2};
    Xcorr = helper.correct_warp(Xpart,M);
    constraint = Xcorr-R-helper.reconstruct(What_FlexMF,Hhat_FlexMF);
    constraints(n,2) = norm(constraint(:),1);
end

%% Plot resutls for differnt algorithms
emds_W_all = zeros(nSim, 3);
emds_H_all = zeros(nSim, 3);
num_seq_all = zeros(nSim, 3);
emds_W_all(:,1) = cellfun(@nanmean, emds_W_SeqNMF);
emds_H_all(:,1) = cellfun(@nanmean, emds_H_SeqNMF);
emds_W_all(:,2:3) = cellfun(@nanmean, emds_W_FlexMF);
emds_H_all(:,2:3) = cellfun(@nanmean, emds_H_FlexMF);
num_seq_all(:,1) = cellfun(@nnz, ids_SeqNMF);
num_seq_all(:,2:3) = cellfun(@nnz, ids_FlexMF);

figure;
ax1 = subplot('Position', [0.1 0.65 0.8 0.3]);
hold on
swarmchart(ones(nSim,1), emds_W_all(:,1))
swarmchart(2*ones(nSim,1), emds_W_all(:,2))
swarmchart(3*ones(nSim,1), emds_W_all(:,3))
set(gca, 'XLabel', [], 'XTicklabel', [])
ylabel(ax1, 'EMDs W')

ax2 = subplot('Position', [0.1 0.3 0.8 0.3]);
hold on
swarmchart(ones(nSim,1), emds_H_all(:,1))
swarmchart(2*ones(nSim,1), emds_H_all(:,2))
swarmchart(3*ones(nSim,1), emds_H_all(:,3))
set(gca, 'XLabel', [], 'XTicklabel', [])
ylabel(ax2, 'EMDs H')

ax3 = subplot('Position', [0.1 0.1 0.8 0.15]);
hold on
errorbar(1:3, median(num_seq_all), ...
median(num_seq_all)-prctile(num_seq_all,25), ...
prctile(num_seq_all,75)-median(num_seq_all), ...
'-', 'Marker', '.', 'MarkerSize', 12, 'Color', 'k');
set(gca, 'XTick', 1:3, 'XTickLabel', {'SeqNMF', 'FlexMF Rand Inti', 'FlexMF SeqNMF Init'})
ylabel('#Sequences')
ylim([0,Khat])
linkaxes([ax1, ax2, ax3], 'x')
set(gcf, 'Position', [100,100,600,800])