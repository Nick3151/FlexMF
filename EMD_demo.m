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
Xshift = generate_sequence(T,N,3,'shift',2);

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

%% Compute EMD demo
opts_default = tfocs_SCD;
opts = opts_default;
opts.continuation = 1;
opts.tol = 1e-6;
opts.stopCrit = 4;
opts.maxIts = 500;
% opts.alg = 'N83';

[d, M, R, out] = compute_EMD(X, Xwarp, opts);
% [d, M, R, out] = compute_EMD(X, Xshift, opts);

figure;
plot_MR(M,R)
save2pdf('EMD_warp_demo.pdf')
% save2pdf('EMD_shift_demo.pdf')

figure;
plot(out.f)
title('out.f')

% Check constraint
D = eye(T) - diag(ones(T-1,1),-1);
% C = M*D'-R-(Xshift-X);
C = M*D'-R-(Xwarp-X);
figure;
imagesc(C)
set(gca, 'XTickLabel', [], 'YTickLabel', []);
title('Constraint error')
colorbar
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

    [d, M, R, out] = compute_EMD(X, Xshift, opts);
    EMDs_shift(j+1) = d;
    b = Xshift-X;
    L2_shift(j+1) = norm(b(:))^2;

    [d, M, R, out] = compute_EMD(X, Xwarp, opts);
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

%% Compute EMD under noise
opts_default = tfocs_SCD;
opts = opts_default;
opts.continuation = 1;
opts.tol = 1e-6;
opts.stopCrit = 4;
opts.maxIts = 500;
% opts.alg = 'N83';

X = generate_sequence(T,N,3, 'len_burst', 5, 'dynamic', 1);
X_noise = generate_sequence(T,N,3, 'noise', .01, 'len_burst', 5, 'dynamic', 1);
X_warp = generate_sequence(T,N,3, 'noise', .01, 'warp', 1, 'len_burst', 5, 'dynamic', 1);
figure; SimpleXplot_patch(X)
save2pdf('Sequence_dynamic.pdf')
figure; SimpleXplot_patch(X_noise)
save2pdf('Sequence_dynamic_noise.pdf')
figure; SimpleXplot_patch(X_warp)
save2pdf('Sequence_dynamic_warp.pdf')

[d, M, R, out] = compute_EMD(X_noise, X_warp, opts);

figure;
plot_MR(M,R)
save2pdf('EMD_noise_dynamic_demo.pdf')

% Check constraint
[N,T] = size(X_warp);
D = eye(T) - diag(ones(T-1,1),-1);
C = M*D'-R-(X_warp-X_noise);
figure;
imagesc(C)
set(gca, 'XTickLabel', [], 'YTickLabel', []);
title('Constraint error')
colorbar
