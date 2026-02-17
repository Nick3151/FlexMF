%% Demo script: EMD between two completely different sequences
clear all
close all
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(fullfile(root, 'Utils'))
addpath(genpath(fullfile(root, 'FlexMF')));
%% Generate some complete different synthetic sequence 
T = 50;
N = 10;
X = generate_sequence(T,N,3);
X1 = zeros(2*N,T);
X2 = zeros(2*N,T);
X1(1:N,:) = X;
X2(N+1:2*N,:) = X;

figure;
SimpleXplot(X1)
save2pdf('Sequence_simulated1.pdf')

figure;
SimpleXplot(X2)
save2pdf('Sequence_simulated2.pdf')

%% Choose lambdaR different sequences
opts_default = tfocs_SCD;
opts = opts_default;
opts.continuation = 1;
opts.tol = 1e-6;
opts.stopCrit = 4;
opts.maxIts = 500;
% opts.alg = 'N83';

nlambdaRs = 20;
lambdaRs = logspace(-1, 3, nlambdaRs);
Ms = cell(nlambdaRs,1);
Rs = cell(nlambdaRs,1);
ds = cell(nlambdaRs,1);
for n=1:nlambdaRs
    disp(n)
    tic
    [ds{n}, Ms{n}, Rs{n}, out] = compute_EMD(X1, X2, opts, 'lambdaR', lambdaRs(n));
    toc
end
M_norms = cellfun(@(x) norm(x(:),1), Ms);
R_norms = cellfun(@(x) norm(x(:),1), Rs);
Err = X1-X2;

figure;
plot(lambdaRs, M_norms, 'r', lambdaRs, R_norms, 'b')
set(gca, 'XScale', 'log')
hold on
yline(norm(Err(:),1), 'k')
xlabel('lambdaR')
legend('L1M', 'L1R', 'L1Err', 'Location', 'best')
save2pdf('EMD_Choose_lambdaR_diff_seq')

figure;
n = 10;
plot_MR(Ms{n},Rs{n})
save2pdf(sprintf('EMD_diff_seq_demo_lambdaR=%0.1f_MR.pdf', lambdaRs(n)))

%% Choose lambdaR reverse sequences
Xrev = flip(X,2);

figure;
SimpleXplot(X)
save2pdf('Sequence_simulated_X.pdf')
figure;
SimpleXplot(Xrev)
save2pdf('Sequence_simulated_Xrev.pdf')

nlambdaRs = 20;
lambdaRs = logspace(-1, 3, nlambdaRs);
Ms = cell(nlambdaRs,1);
Rs = cell(nlambdaRs,1);
ds = cell(nlambdaRs,1);
for n=1:nlambdaRs
    disp(n)
    tic
    [ds{n}, Ms{n}, Rs{n}, out] = compute_EMD(X, Xrev, opts, 'lambdaR', lambdaRs(n));
    toc
end
M_norms = cellfun(@(x) norm(x(:),1), Ms);
R_norms = cellfun(@(x) norm(x(:),1), Rs);
Err = X-Xrev;

figure;
plot(lambdaRs, M_norms, 'r', lambdaRs, R_norms, 'b')
set(gca, 'XScale', 'log')
hold on
yline(norm(Err(:),1), 'k')
xlabel('lambdaR')
legend('L1M', 'L1R', 'L1Err', 'Location', 'best')
save2pdf('EMD_Choose_lambdaR_rev_seq')

figure;
n = 6;
plot_MR(Ms{n},Rs{n})
save2pdf(sprintf('EMD_rev_seq_demo_lambdaR=%0.1f_MR.pdf', lambdaRs(n)))

%% Choose lambdaR same sequences with noise
X_noise = X+(rand(size(X))<.01);
nlambdaRs = 20;
lambdaRs = logspace(-1, 3, nlambdaRs);
Ms = cell(nlambdaRs,1);
Rs = cell(nlambdaRs,1);
ds = cell(nlambdaRs,1);
for n=1:nlambdaRs
    disp(n)
    tic
    [ds{n}, Ms{n}, Rs{n}, out] = compute_EMD(X, X_noise, opts, 'lambdaR', lambdaRs(n));
    toc
end
M_norms = cellfun(@(x) norm(x(:),1), Ms);
R_norms = cellfun(@(x) norm(x(:),1), Rs);
Err = X-X_noise;

figure;
plot(lambdaRs, M_norms, 'r', lambdaRs, R_norms, 'b')
set(gca, 'XScale', 'log')
hold on
yline(norm(Err(:),1), 'k')
xlabel('lambdaR')
legend('L1M', 'L1R', 'L1Err', 'Location', 'best')
save2pdf('EMD_Choose_lambdaR_same_seq_noise')