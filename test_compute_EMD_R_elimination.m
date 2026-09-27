%% Compare pre- and post-R-elimination compute_EMD: accuracy and speed
% compute_EMD.m now solves for M directly (R is eliminated and recovered
% algebraically); the legacy stacked-[M;R] solver is kept only as a local
% function below, purely for this before/after comparison.
clear
close all

root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(fullfile(root, 'Utils'))
addpath(genpath(fullfile(root, 'FlexMF')))

opts = tfocs_SCD;
opts.continuation = 1;
opts.tol = 1e-6;
opts.stopCrit = 4;
opts.maxIts = 500;
opts.printEvery = 0;
opts.alg = 'N83';
copts = continuation();
copts.verbose = 0;
lambdaR = 1e2;

%% Test cases: small shift/warp/noise (as in EMD_demo.m) plus large true/wrong pairs (as in diagnose_uot_scale.m)
T = 50;
N = 10;
X = generate_sequence(T,N,3);
Xshift = generate_sequence(T,N,3,'shift',2);
Xwarp = generate_sequence(T,N,3,'warp',1);
Xwarp_noise = Xwarp;
Xwarp_noise(1:5, end-3) = 1;

K = 3;
Tbig = 2000;
[~, Wbig] = generate_data(Tbig, 5*ones(K,1), 3*ones(K,1), 'noise', 0, ...
    'participation', 0.8*ones(K,1), 'seed', 1, 'len_burst', 1, 'dynamic', 0);
Lbig = size(Wbig, 3);
w1 = reshape(Wbig(:,1,:), size(Wbig,1), Lbig);
w2 = reshape(Wbig(:,2,:), size(Wbig,1), Lbig);

cases = struct('name', {'shift', 'warp', 'warp+noise', 'large (true match)', 'large (disjoint)'}, ...
    'X1', {X, X, X, w1, w1}, ...
    'X2', {Xshift, Xwarp, Xwarp_noise, helper.shift_profiles(w1,3,Lbig), w2});

%% Run old vs new on every case
% d_old is the legacy solver's own (possibly infeasible) reported cost;
% d_old_honest recomputes it from M_old using the true constraint-derived R,
% which is the fair comparison against d_new since d_new can never cheat via
% infeasibility (R is always derived algebraically from M).
fprintf('%-20s %10s %10s %8s %12s %12s %12s %10s %10s %11s %11s\n', ...
    'case', 't_old(s)', 't_new(s)', 'speedup', 'd_old', 'd_old_honest', 'd_new', 'reldiff_M', 'infeas_old', 'infeas_new', 'improved');

for i = 1:numel(cases)
    X1 = cases(i).X1;
    X2 = cases(i).X2;
    Ti = size(X1,2);
    D = eye(Ti) - diag(ones(Ti-1,1),-1);
    D(Ti,Ti) = 0;
    b = X2 - X1;

    tic
    [d_old, M_old, R_old] = compute_EMD_old(X1, X2, opts, lambdaR, copts);
    t_old = toc;

    tic
    [d_new, M_new, R_new] = compute_EMD(X1, X2, opts, 'lambdaR', lambdaR, 'continuationOptions', copts);
    t_new = toc;

    infeas_old = norm(M_old*D' - R_old - b, 'fro') / max(norm(b(:)), eps);
    infeas_new = norm(R_new - (M_new*D' + X1 - X2), 'fro');

    R_old_true = M_old*D' + X1 - X2;
    d_old_honest = norm(M_old(:),1) + lambdaR*norm(R_old_true(:),1);

    reldiff_M = norm(M_new(:) - M_old(:)) / max(norm(M_old(:)), eps);

    fprintf('%-20s %10.4f %10.4f %8.2f %12.4f %12.4f %12.4f %10.4f %11.2e %11.2e %11s\n', ...
        cases(i).name, t_old, t_new, t_old / max(t_new, eps), d_old, d_old_honest, d_new, ...
        reldiff_M, infeas_old, infeas_new, mat2str(d_new <= d_old_honest*(1+1e-6)));
end

%% Legacy pre-elimination solver (stacked [M;R], hard equality constraint via proj_Rn)
% Kept only for this comparison; production code now uses compute_EMD.m.
function [d, M, R, out] = compute_EMD_old(X1, X2, opts, lambdaR, continuationOptions)
[N,T] = size(X1);
if isequal(X1,X2)
    d = 0;
    M = zeros(N,T);
    R = zeros(N,T);
    out = [];
    return
end
A = @(Y, mode)Beckmann_UOT_constraint(N, T, Y, mode);
W = @(Y, mode)Beckmann_UOT_obj(N, T, lambdaR, Y, mode);
b = X2-X1;
[Y, out] = solver_sBPDN_W(A,W,b,0,.1,[],[],opts, continuationOptions);
M = Y(1:N,:);
R = Y(N+1:2*N,:);
d = norm(M(:),1) + lambdaR*norm(R(:),1);
end
