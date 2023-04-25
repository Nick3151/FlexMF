clear all
close all
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
rmpath(genpath(fullfile(root, 'seqNMF-master')));
addpath(genpath(fullfile(root, 'FlexMF')));
%% Generate some synthetic data with temporal jittering or time warping
number_of_seqences = 3;
T = 800; % length of data to generate
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
[X, W, H, ~] = generate_data(T,Nneurons,Dt, 'seed', seed);
[Xwarp, Wwarp, Hwarp, ~] = generate_data(T,Nneurons,Dt, 'warp', warp, 'seed', seed);
[Xjit, Wjit, Hjit, ~] = generate_data(T,Nneurons,Dt, 'jitter', jitter, 'seed', seed);
L = size(W,3);
% range = round(L/2)-25:round(L/2)+35;

plotAll = 1;
figure; SimpleWHPlot_patch(W,H,[],[],[],X,plotAll); title('generated data raw','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot_patch(Wwarp,Hwarp,[],[],[],Xwarp,plotAll); title('generated data warping','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot_patch(Wjit,Hjit,[],[],[],Xjit,plotAll); title('generated data jittering','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

%% Compute EMD
opts_default = tfocs_SCD;
opts = opts_default;
opts.continuation = 1;
opts.tol = 1e-6;
opts.stopCrit = 4;
opts.maxIts = 5000;

mu = 1e6;
N = size(X,1);
A = @(Y, mode)Beckmann_UOT_constraint(N, T, Y, mode);
W = @(Y, mode)Beckmann_UOT_obj(N, T, mu, Y, mode);
% b = Xwarp-X;
b = Xjit-X;
[Y, out] = solver_sBPDN_W(A,W,b,0,1,[],[],opts);

M = Y(1:N,:);
R = Y(N+1:2*N,:);
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

figure;
plot(out.f)
title('out.f')

emd = norm(M(:),1)+mu*norm(R(:),1)

% Check constraint
D = eye(T) - diag(ones(T-1,1),-1);
C = M*D'-R-b;
figure;
imagesc(C)
set(gca, 'XTickLabel', [], 'YTickLabel', []);
colorbar

%% EMD vs different levels of jittering
number_of_seqences = 3;
T = 800; % length of data to generate
Nneurons = 10*ones(number_of_seqences,1); % number of neurons in each sequence
Dt = 3.*ones(number_of_seqences,1); % gap between each member of the sequence
noise = 0.001; % probability of added noise in each bin
participation = 1.*ones(number_of_seqences,1); % Participation parameter = 100%
gap = 100;
neg = 0;
bin = 0;
[X, W, H, ~] = generate_data(T,Nneurons,Dt,noise,zeros(number_of_seqences,1),participation,gap,0,bin,neg,1);
L = size(W,3);

opts_default = tfocs_SCD;
opts = opts_default;
opts.continuation = 1;
opts.tol = 1e-6;
opts.stopCrit = 4;
opts.maxIts = 5000;

mu = 1e6;
N = size(X,1);
A = @(Y, mode)Beckmann_UOT_constraint(N, T, Y, mode);
W = @(Y, mode)Beckmann_UOT_obj(N, T, mu, Y, mode);

for j=1:5
    jitter = j*ones(number_of_seqences,1); % Jitter std
    [Xjit, Wjit, Hjit, ~] = generate_data(T,Nneurons,Dt,noise,jitter,participation,gap,0,bin,neg,1);
    b = Xjit-X;
    [Y, out] = solver_sBPDN_W(A,W,b,0,1,[],[],opts);
    M = Y(1:N,:);
    R = Y(N+1:2*N,:);
    EMDs(j) = norm(M(:),1)+mu*norm(R(:),1);
end

figure;
plot(1:5, EMDs, 'LineWidth',2)
ylabel('EMD')
xlabel('Jitter SD')

%% EMD vs different levels of warping
number_of_seqences = 3;
T = 800; % length of data to generate
Nneurons = 10*ones(number_of_seqences,1); % number of neurons in each sequence
Dt = 3.*ones(number_of_seqences,1); % gap between each member of the sequence
noise = 0.001; % probability of added noise in each bin
jitter = 0*ones(number_of_seqences,1); % Jitter std
participation = 1.*ones(number_of_seqences,1); % Participation parameter = 100%
gap = 100;
neg = 0;
bin = 0;
[X, W, H, ~] = generate_data(T,Nneurons,Dt,noise,jitter,participation,gap,0,bin,neg,1);
L = size(W,3);

opts_default = tfocs_SCD;
opts = opts_default;
opts.continuation = 1;
opts.tol = 1e-6;
opts.stopCrit = 4;
opts.maxIts = 5000;

mu = 1e6;
N = size(X,1);
A = @(Y, mode)Beckmann_UOT_constraint(N, T, Y, mode);
W = @(Y, mode)Beckmann_UOT_obj(N, T, mu, Y, mode);

for warp=0:2
    [Xwarp, Wwarp, Hwarp, ~] = generate_data(T,Nneurons,Dt,noise,jitter,participation,gap,warp,bin,neg,1);
    b = Xwarp-X;
    [Y, out] = solver_sBPDN_W(A,W,b,0,1,[],[],opts);
    M = Y(1:N,:);
    R = Y(N+1:2*N,:);
    EMDs(warp) = norm(M(:),1)+mu*norm(R(:),1);
end

figure;
plot(0:2, EMDs, 'LineWidth',2)
ylabel('EMD')
xlabel('Warping level')