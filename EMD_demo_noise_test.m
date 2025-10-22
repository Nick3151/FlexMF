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

[Xwarp, Wwarp, Hwarp, ~] = generate_data(T,Nneurons,Dt, 'noise',noise, 'warp', warp, 'seed', seed, 'len_burst', 1, 'dynamic', 0);

L = size(Wwarp,3);

plotAll = 1;
figure; SimpleWHPlot(Wwarp,Hwarp,'Data',Xwarp,'plotAll', plotAll); title('generated data warping','Fontsize',16)

%% Normalize data
K = 3;

%% Fix W with ground truth, fit HMR, parameter sweep on lambdaL1H
lambdaL1Hs = .01:.01:.1;
lambda = 0;
lambda_M = .5;
lambda_R = 1;

nLambdas = length(lambdaL1Hs);
Whats = cell(nLambdas,1);
Hhats = cell(nLambdas,1);
errors = zeros(nLambdas, 4);
Ms = cell(nLambdas,1);
Rs = cell(nLambdas,1);

for n=1:nLambdas
    lambdaL1H = lambdaL1Hs(n);
    disp(lambdaL1H)
    
    [Whats{n}, Hhats{n}, cost, error, loadings, power, Ms{n}, Rs{n}] = FlexMF(Xwarp, 'K', K, 'L', L, ...
        'EMD',1, 'lambda', lambda, ...
        'lambdaL1H', lambdaL1H, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'maxiter', 1, 'tolerance', 1e-4, ...
        'W_init', Wwarp, 'W_fixed', 1, 'showPlot', 0, 'verbal', 0);
    errors(n,:) = error(end,:);
end

%% Look at factors
plotAll = 1;
n = 3;
figure; SimpleWHPlot(Whats{n}, Hhats{n}, 'plotAll', plotAll); title('FlexMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot(Whats{n}, Hhats{n}, 'Data', Xwarp, 'plotAll', plotAll); title('FlexMF factors, with raw data')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure;
plot_MR(Ms{n},Rs{n})

%% Errors vs lambdaL1H
L1Hs = cellfun(@(x) norm(x(:),1), Hhats);
figure;
plot(lambdaL1Hs, errors(:,1), lambdaL1Hs, L1Hs)
legend('Recon', 'L1H')
xlabel('LambdaL1H')
save2pdf('Error vs lambdaL1H')
