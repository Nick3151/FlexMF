%% Demo script showing EMD between two sequences
clear all
close all
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(fullfile(root, 'Utils'))
rmpath(genpath(fullfile(root, 'seqNMF-master')));
addpath(genpath(fullfile(root, 'FlexMF')));
%% Generate some synthetic sequence with temporal shifting or time warping
T = 50;
N = 10;
X = generate_sequence(T,N,3);
Xwarp = generate_sequence(T,N,3,'warp',1);
Xwarp_noise = Xwarp;
Xwarp_noise(1:5, end-3) = 1;
Xshift = generate_sequence(T,N,3,'shift',2);
Xshift_noise = Xshift;
Xshift_noise(1:5, end-3) = 1;

figure; 
subplot('Position', [0.05 0.65 0.9 0.3])
SimpleXplot(X)
set(gca, 'Position', [0.05 0.65 0.9 0.3], 'XTick', [], 'YTick', [], 'xlabel', [], 'ylabel', [])
subplot('Position', [0.05 0.35 0.9 0.3])
SimpleXplot(Xwarp)
set(gca, 'Position', [0.05 0.35 0.9 0.3], 'XTick', [], 'YTick', [], 'xlabel', [], 'ylabel', [])
subplot('Position', [0.05 0.05 0.9 0.3])
SimpleXplot(Xshift)
set(gca, 'Position', [0.05 0.05 0.9 0.3], 'XTick', [], 'YTick', [], 'xlabel', [], 'ylabel', [])
save2pdf('Sequence_simulated.pdf')

figure; 
subplot('Position', [0.05 0.65 0.9 0.3])
SimpleXplot(X)
set(gca, 'Position', [0.05 0.65 0.9 0.3], 'XTick', [], 'YTick', [], 'xlabel', [], 'ylabel', [])
subplot('Position', [0.05 0.35 0.9 0.3])
SimpleXplot(Xwarp_noise)
set(gca, 'Position', [0.05 0.35 0.9 0.3], 'XTick', [], 'YTick', [], 'xlabel', [], 'ylabel', [])
subplot('Position', [0.05 0.05 0.9 0.3])
SimpleXplot(Xshift_noise)
set(gca, 'Position', [0.05 0.05 0.9 0.3], 'XTick', [], 'YTick', [], 'xlabel', [], 'ylabel', [])
save2pdf('Sequence_simulated_noise.pdf')

%% Compute EMD demo
opts_default = tfocs_SCD;
opts = opts_default;
opts.continuation = 1;
opts.tol = 1e-6;
opts.stopCrit = 4;
opts.maxIts = 500;
% opts.alg = 'N83';

tic
[d, M, R, out] = compute_EMD(X, Xwarp, opts);
% [d, M, R, out] = compute_EMD(X, Xshift, opts);
toc

figure;
plot_MR(M,R)
save2pdf('EMD_warp_demo_MR.pdf')
% save2pdf('EMD_shift_demo_MR.pdf')

figure;
plot(out.f)
title('out.f')

% Check constraint
D = eye(T) - diag(ones(T-1,1),-1);
D(T,T) = 0;
% C = M*D'-R-(Xshift-X);
C = M*D'-R-(Xwarp-X);
figure;
imagesc(C)
set(gca, 'XTickLabel', [], 'YTickLabel', []);
title('Constraint error')
colorbar

%% Choose lambdaR
nlambdaRs = 20;
lambdaRs = logspace(-1, 3, nlambdaRs);
Ms = cell(nlambdaRs,1);
Rs = cell(nlambdaRs,1);
ds = cell(nlambdaRs,1);
for n=1:nlambdaRs
    disp(n)
    tic
%     [ds{n}, Ms{n}, Rs{n}, out] = compute_EMD(X, Xshift_noise, opts, 'lambdaR', lambdaRs(n));
    [ds{n}, Ms{n}, Rs{n}, out] = compute_EMD(X, Xwarp_noise, opts, 'lambdaR', lambdaRs(n));
    toc
end
M_norms = cellfun(@(x) norm(x(:),1), Ms);
R_norms = cellfun(@(x) norm(x(:),1), Rs);
% Err = X-Xshift_noise;
Err = X-Xwarp_noise;

figure;
plot(lambdaRs, M_norms, 'r', lambdaRs, R_norms, 'b')
set(gca, 'XScale', 'log')
hold on
yline(norm(Err(:),1), 'k')
xlabel('lambdaR')
legend('L1M', 'L1R', 'L1Err', 'Location', 'best')
save2pdf('EMD_Choose_lambdaR_warp_noise')
% save2pdf('EMD_Choose_lambdaR_shift_noise')

figure;
n = 10;
plot_MR(Ms{n},Rs{n})
% save2pdf(sprintf('EMD_shift_noise_demo_lambdaR=%0.1f_MR.pdf', lambdaRs(n)))
save2pdf(sprintf('EMD_warp_noise_demo_lambdaR=%0.1f_MR.pdf', lambdaRs(n)))

%% EMD vs different levels of warping/shift
T = 100;
N = 10;
X = generate_sequence(T,N,3);

EMDs_shift = zeros(5,1);
L2_shift = zeros(5,1);
EMDs_warp = zeros(5,1);
L2_warp = zeros(5,1);

for j=0:4
    Xshift = generate_sequence(T,N,3,'shift',j);
    Xwarp = generate_sequence(T,N,3,'warp',j);

    [d, M, R, out] = compute_EMD(X, Xshift, opts, 'lambdaR', 1e2);
    EMDs_shift(j+1) = d;
    b = Xshift-X;
    L2_shift(j+1) = norm(b(:))^2;

    [d, M, R, out] = compute_EMD(X, Xwarp, opts, 'lambdaR', 1e2);
    EMDs_warp(j+1) = d;
    b = Xwarp-X;
    L2_warp(j+1) = norm(b(:))^2;
end

figure;
plot(0:4, EMDs_shift, 'LineWidth',2)
hold on
plot(0:4, L2_shift, 'LineWidth',2)
ylabel('Distance')
xlabel('Shift step')
legend({'EMD', 'L2-square'})
save2pdf('EMD_vs_shift.pdf')

figure;
plot(0:4, EMDs_warp, 'LineWidth',2)
hold on
plot(0:4, L2_warp, 'LineWidth',2)
ylabel('Distance')
xlabel('Warp step')
legend({'EMD', 'L2-square'})
save2pdf('EMD_vs_warp.pdf')
