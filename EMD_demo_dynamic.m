clear all
close all
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(fullfile(root, 'Utils'))
rmpath(genpath(fullfile(root, 'seqNMF-master')));
addpath(genpath(fullfile(root, 'FlexMF')));

opts_default = tfocs_SCD;
opts = opts_default;
opts.continuation = 1;
opts.tol = 1e-6;
opts.stopCrit = 4;
opts.maxIts = 500;
% opts.alg = 'N83';

%% Generate sequences
T = 100;
N = 10;
X = generate_sequence(T,N,3, 'len_burst', 5, 'dynamic', 1);
X_noise = generate_sequence(T,N,3, 'noise', .01, 'len_burst', 5, 'dynamic', 1);
Xwarp_noise = generate_sequence(T,N,3, 'noise', .01, 'warp', 1, 'len_burst', 5, 'dynamic', 1);
figure; SimpleXplot_patch(X)
save2pdf('Sequence_dynamic.pdf')
figure; SimpleXplot_patch(X_noise)
save2pdf('Sequence_dynamic_noise.pdf')
figure; SimpleXplot_patch(Xwarp_noise)
save2pdf('Sequence_dynamic_warp.pdf')

%% Choose lambdaR
nlambdaRs = 20;
lambdaRs = logspace(-1, 3, nlambdaRs);
Ms = cell(nlambdaRs,1);
Rs = cell(nlambdaRs,1);
ds = cell(nlambdaRs,1);
for n=1:nlambdaRs
    disp(n)
    tic
    [ds{n}, Ms{n}, Rs{n}, out] = compute_EMD(X, Xwarp_noise, opts, 'lambdaR', lambdaRs(n));
    toc
end
M_norms = cellfun(@(x) norm(x(:),1), Ms);
R_norms = cellfun(@(x) norm(x(:),1), Rs);
Err = X-Xwarp_noise;

figure;
plot(lambdaRs, M_norms, 'r', lambdaRs, R_norms, 'b')
set(gca, 'XScale', 'log')
hold on
yline(norm(Err(:),1), 'k')
xlabel('lambdaR')
legend('L1M', 'L1R', 'L1Err', 'Location', 'best')
save2pdf('EMD_Choose_lambdaR_dynamic_warp_noise')

figure;
n = 12;
plot_MR(Ms{n},Rs{n})
save2pdf(sprintf('EMD_dynamic_warp_noise_demo_lambdaR=%0.1f_MR.pdf', lambdaRs(n)))

%% Compute EMD under noise and calcium dynamics
[d, M, R, out] = compute_EMD(X_noise, Xwarp_noise, opts, 'lambdaR', 10);

figure;
plot_MR(M,R)
% save2pdf('EMD_dynamic_warp_noise_demo_MR.pdf')

% Check constraint
[N,T] = size(Xwarp_noise);
D = eye(T) - diag(ones(T-1,1),-1);
D(T,T) = 0;
C = M*D'-R-(Xwarp_noise-X_noise);
figure;
imagesc(C)
set(gca, 'XTickLabel', [], 'YTickLabel', []);
title('Constraint error')
colorbar
