clear all
close all
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
%% Generate some synthetic data with temporal jittering or time warping
number_of_seqences = 3;
T = 800; % length of data to generate
Nneurons = 10*ones(number_of_seqences,1); % number of neurons in each sequence
Dt = 3.*ones(number_of_seqences,1); % gap between each member of the sequence
NeuronNoise = 0.001; % probability of added noise in each bin
SeqNoiseTime = 0.25*ones(number_of_seqences,1); % Jitter parameter = 25%
SeqNoiseNeuron = 1.*ones(number_of_seqences,1); % Participation parameter = 100%
stretch = 2; % stretch should be less than Dt
neg = 0;
bin = 0;
[X, W, H, ~] = generate_data(T,Nneurons,Dt,NeuronNoise,zeros(number_of_seqences,1),SeqNoiseNeuron,0,bin,neg,1);
[Xwarp, Wwarp, Hwarp, ~] = generate_data(T,Nneurons,Dt,NeuronNoise,zeros(number_of_seqences,1),SeqNoiseNeuron,stretch,bin,neg,1);
[Xjit, Wjit, Hjit, ~] = generate_data(T,Nneurons,Dt,NeuronNoise,SeqNoiseNeuron,SeqNoiseNeuron,0,bin,neg,1);
L = size(W,3);
range = round(L/2)-25:round(L/2)+35;

figure; SimpleWHPlot(W(:,:,range),H,X); title('generated data raw','Fontsize',16)
set(gcf,'position',[200,200,1600,900])
figure; SimpleWHPlot(Wwarp(:,:,range),Hwarp,Xwarp); title('generated data warping','Fontsize',16)
set(gcf,'position',[200,200,1600,900])
figure; SimpleWHPlot(Wjit(:,:,range),Hjit,Xjit); title('generated data jittering','Fontsize',16)
set(gcf,'position',[200,200,1600,900])

%% Compute EMD
mu = 1;
N = size(X,1);
A = @(Y, mode)Beckmann_UOT_constraint(N+1, T, Y, mode);
W = @(Y, mode)Beckmann_UOT_obj(N+1, T, mu, Y, mode);
b = Xwarp-X;
[Y, out] = solver_sBPDN_W(A,W,b,0,0);