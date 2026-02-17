%% Test FlexMF parameters (Fix lambda_R, sweep lambda, lambda_M)
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
noise = 0.001; % probability of added noise in each bin
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

%% Visualize data
figure; SimpleWHPlot(Wwarp,Hwarp,'Data',Xwarp,'plotAll', 1); title('generated data warping','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

%% Normalize data
K = 3;
% frob_norm = norm(Xwarp(:));
% Xwarp = Xwarp/frob_norm*K;
% Wwarp = Wwarp/frob_norm*K;

%% Procedure for choosing lambda in SeqNMF
nLambdas = 20; % increase if you're patient
lambdas = sort(logspace(0,-4,nLambdas), 'ascend'); 
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
lambda = .05;
lambdaL1H = 0;
lambdaL1W = 0;
lambdaOrthoH = 0;

figure;
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
[What_SeqNMF, Hhat_SeqNMF, ~, errors_SeqNMF,loadings,power]= seqNMF(Xwarp,'K',K,'L',L,...
            'lambda', lambda, 'maxiter', 50, 'showPlot', 1); 

%% Run FlexMF with different lambda and lambda_M
n = 9;
lambda_R = 1;
lambdas = logspace(-4,0,n);
lambda_Ms = logspace(-3,1,n);

Whats = cell(n,n);
Hhats = cell(n,n);
Ms = cell(n,n);
Rs = cell(n,n);
emds_W = cell(n, n);
emds_H = cell(n, n);
num_detected = cell(n, n);
ids_match = cell(n, n);

for Li = 1:n
    lambda = lambdas(Li);
    parfor Mi = 1:n       
        lambda_M = lambda_Ms(Mi);
        display(['Testing lambda=' num2str(lambda) ' lambda_M=' num2str(lambda_M)])
        tic
        [What, Hhat, cost, errors, loadings, power, M, R] = FlexMF(Xwarp, 'K', K, 'L', L, ...
            'EMD',1, 'lambda', lambda, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'maxiter', 50, 'tolerance', 1e-4,...
            'W_init', What_SeqNMF, 'H_init', Hhat_SeqNMF, 'showPlot', 0, 'verbal', 0);
        toc
        Whats{Li,Mi} = What;
        Hhats{Li,Mi} = Hhat;
        Ms{Li,Mi} = M;
        Rs{Li,Mi} = R;
        
        disp('Evaluate EMDs of results')
        tic
        [emds_W{Li,Mi}, emds_H{Li,Mi}, ids] = helper.similarity_WH_EMD(Wwarp, Hwarp, What, Hhat);
        num_detected{Li,Mi} = length(ids);
        ids_match{Li,Mi} = ids;
        toc
    end
end

% save(sprintf('Simulation_Results/EMD_params_warp.mat'), ...
%     "Hhats", "Whats", "Ms", "Rs", "lambda_Ms", "lambdas", "Xwarp")
save(sprintf('Simulation_Results/EMD_params_warp_noise.mat'), ...
    "Hhats", "Whats", "Ms", "Rs", "lambda_Ms", "lambdas", "Xwarp")
%% Look at factors
i = 6;
j = 4;
plotAll = 1;
figure; SimpleWHPlot_patch(Whats{i,j}, Hhats{i,j}, 'plotAll', plotAll); title('FlexMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

figure;
plot_MR(Ms{i,j},Rs{i,j})

%% Visualize EMDs to ground truth
mean_emd_Hs = zeros(n, n);
mean_emd_Ws = zeros(n, n);
num_detected_all = zeros(n, n);

for Li=1:n
    for Mi=1:n
        mean_emd_Hs(Li,Mi) = mean(emds_H{Li,Mi});
        mean_emd_Ws(Li,Mi) = mean(emds_W{Li,Mi});
        num_detected_all(Li,Mi)  = num_detected{Li,Mi};
    end
end

M_norms = cellfun(@(x) norm(x(:),1), Ms);
R_norms = cellfun(@(x) norm(x(:),1), Rs);