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
export_vector_pdf('Sequence_simulated1.pdf');

figure;
SimpleXplot(X2)
export_vector_pdf('Sequence_simulated2.pdf');

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
constraint_rel = zeros(nlambdaRs,1);
Ddiff = eye(T) - diag(ones(T-1,1), -1);
Ddiff(T,T) = 0;
b = X2 - X1;
for n=1:nlambdaRs
    disp(n)
    tic
    [ds{n}, Ms{n}, Rs{n}, out] = compute_EMD(X1, X2, opts, 'lambdaR', lambdaRs(n));
    toc
    C = Ms{n}*Ddiff' - Rs{n} - b;
    constraint_rel(n) = norm(C(:),1) / max(norm(b(:),1), eps);
end
M_norms = cellfun(@(x) norm(x(:),1), Ms);
R_norms = cellfun(@(x) norm(x(:),1), Rs);
Err = X1-X2;

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
export_vector_pdf('EMD_Choose_lambdaR_diff_seq');

figure;
n = 10;
plot_MR(Ms{n},Rs{n}, [], 'imagesc')
export_vector_pdf(sprintf('EMD_diff_seq_demo_lambdaR=%0.1f_MR.pdf', lambdaRs(n)));

%% Choose lambdaR reverse sequences
Xrev = flip(X,2);

figure;
SimpleXplot(X)
export_vector_pdf('Sequence_simulated_X.pdf');
figure;
SimpleXplot(Xrev)
export_vector_pdf('Sequence_simulated_Xrev.pdf');

nlambdaRs = 20;
lambdaRs = logspace(-1, 3, nlambdaRs);
Ms = cell(nlambdaRs,1);
Rs = cell(nlambdaRs,1);
ds = cell(nlambdaRs,1);
constraint_rel = zeros(nlambdaRs,1);
b = Xrev - X;
for n=1:nlambdaRs
    disp(n)
    tic
    [ds{n}, Ms{n}, Rs{n}, out] = compute_EMD(X, Xrev, opts, 'lambdaR', lambdaRs(n));
    toc
    C = Ms{n}*Ddiff' - Rs{n} - b;
    constraint_rel(n) = norm(C(:),1) / max(norm(b(:),1), eps);
end
M_norms = cellfun(@(x) norm(x(:),1), Ms);
R_norms = cellfun(@(x) norm(x(:),1), Rs);
Err = X-Xrev;

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
export_vector_pdf('EMD_Choose_lambdaR_rev_seq');

figure;
n = 6;
plot_MR(Ms{n},Rs{n}, [], 'imagesc')
export_vector_pdf(sprintf('EMD_rev_seq_demo_lambdaR=%0.1f_MR.pdf', lambdaRs(n)));

%% Choose lambdaR same sequences with noise
X_noise = X+(rand(size(X))<.05);
figure;
SimpleXplot(X_noise)
export_vector_pdf('Sequence_simulated_Xnoise.pdf');

nlambdaRs = 20;
lambdaRs = logspace(-1, 3, nlambdaRs);
Ms = cell(nlambdaRs,1);
Rs = cell(nlambdaRs,1);
ds = cell(nlambdaRs,1);
constraint_rel = zeros(nlambdaRs,1);
b = X_noise - X;
for n=1:nlambdaRs
    disp(n)
    tic
    [ds{n}, Ms{n}, Rs{n}, out] = compute_EMD(X, X_noise, opts, 'lambdaR', lambdaRs(n));
    toc
    C = Ms{n}*Ddiff' - Rs{n} - b;
    constraint_rel(n) = norm(C(:),1) / max(norm(b(:),1), eps);
end
M_norms = cellfun(@(x) norm(x(:),1), Ms);
R_norms = cellfun(@(x) norm(x(:),1), Rs);
Err = X-X_noise;

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
export_vector_pdf('EMD_Choose_lambdaR_same_seq_noise');

figure;
n = 10;
plot_MR(Ms{n},Rs{n}, [], 'imagesc')
export_vector_pdf(sprintf('EMD_same_seq_noise_demo_lambdaR=%0.1f_MR.pdf', lambdaRs(n)));
