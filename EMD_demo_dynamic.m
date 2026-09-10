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
export_vector_pdf('Sequence_dynamic.pdf');
figure; SimpleXplot_patch(X_noise)
export_vector_pdf('Sequence_dynamic_noise.pdf');
figure; SimpleXplot_patch(Xwarp_noise)
export_vector_pdf('Sequence_dynamic_warp.pdf');

%% Choose lambdaR
nlambdaRs = 20;
lambdaRs = logspace(-1, 3, nlambdaRs);
Ms = cell(nlambdaRs,1);
Rs = cell(nlambdaRs,1);
ds = cell(nlambdaRs,1);
constraint_rel = zeros(nlambdaRs,1);
Ddiff = eye(T) - diag(ones(T-1,1), -1);
Ddiff(T,T) = 0;
b = Xwarp_noise - X;
for n=1:nlambdaRs
    disp(n)
    tic
    [ds{n}, Ms{n}, Rs{n}, out] = compute_EMD(X, Xwarp_noise, opts, 'lambdaR', lambdaRs(n));
    toc
    C = Ms{n}*Ddiff' - Rs{n} - b;
    constraint_rel(n) = norm(C(:),1) / max(norm(b(:),1), eps);
end
M_norms = cellfun(@(x) norm(x(:),1), Ms);
R_norms = cellfun(@(x) norm(x(:),1), Rs);
Err = X-Xwarp_noise;

figure;
yyaxis left
plot(lambdaRs, M_norms, 'r-', lambdaRs, R_norms, 'b-')
hold on
yline(norm(Err(:),1), 'k--')
ylabel('L1')
yyaxis right
plot(lambdaRs, constraint_rel, 'g-', 'LineWidth', 1.5)
ylabel('||constraint||_1 / ||X2-X1||_1')
set(gca, 'XScale', 'log')
xlabel('lambdaR')
legend('L1M', 'L1R', 'L1Err', 'constraint_{rel}', 'Location', 'best')
export_vector_pdf('EMD_Choose_lambdaR_dynamic_warp_noise');

figure;
n = 12;
plot_MR(Ms{n},Rs{n})
export_vector_pdf(sprintf('EMD_dynamic_warp_noise_demo_lambdaR=%0.1f_MR.pdf', lambdaRs(n)));

%% Compute EMD under noise and calcium dynamics
[d, M, R, out] = compute_EMD(X_noise, Xwarp_noise, opts, 'lambdaR', 10);

figure;
plot_MR(M,R)
% export_vector_pdf('EMD_dynamic_warp_noise_demo_MR.pdf')

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
