% DIAGNOSE_UOT_SCALE
%
% Why does matching on the raw UOT cost fail, and what exactly does the
% unit-mass normalisation change?
%
% Beckmann_UOT_constraint builds A([M; R]) = M*D' - R with
% D = I - subdiag(1), D(T,T) = 0, so the divergence acts along time only and
% never couples neurons. Summing that divergence over t cancels, so at
% feasibility every neuron satisfies sum_t R(n,t) = sum_t b(n,t), giving the
% exact lower bound
%
%   ||R||_1 >= sum_n |mass2(n) - mass1(n)|
%
% For two motifs on disjoint neurons that bound is mass1 + mass2, i.e. the cost
% of a hopeless pair is extensive in the masses involved. That is the suspected
% reason the raw cost is not comparable across pairs.
%
% This script checks the claim on controlled inputs instead of on seqNMF
% output: a true pair (a motif against a shifted copy of itself) and a wrong
% pair (two motifs on disjoint neurons), each at several overall scales, and
% each at several solver tolerances. It reports whether ||R||_1 reaches its
% lower bound and how far the returned point is from feasibility, because
% compute_EMD calls solver_sBPDN_W with a fixed smoothing parameter mu = 0.1
% that is not scale invariant.

clear

root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(fullfile(root, 'Utils'))
addpath(genpath(fullfile(root, 'FlexMF')))

K = 3;
T = 2000;
[~, W, H, ~] = generate_data(T, 5 * ones(K, 1), 3 * ones(K, 1), 'noise', 0, ...
    'participation', 0.8 * ones(K, 1), 'seed', 1, ...
    'len_burst', 1, 'dynamic', 0);
L = size(W, 3);
N = size(W, 1);

w1 = reshape(W(:, 1, :), N, L);
w2 = reshape(W(:, 2, :), N, L);
w1 = w1 / sum(w1(:));            % unit mass reference profiles
w2 = w2 / sum(w2(:));

fprintf('neuron support: motif 1 uses %s\n', mat2str(find(any(w1, 2))'));
fprintf('neuron support: motif 2 uses %s\n', mat2str(find(any(w2, 2))'));
fprintf('overlap: %d neurons\n\n', nnz(any(w1, 2) & any(w2, 2)));

lambdaR = 1e3;
scales = [1, 5, 20];             % total mass fed to compute_EMD
tols = [1e-4, 1e-6];

D = eye(L) - diag(ones(L - 1, 1), -1);
D(L, L) = 0;

fprintf('%-14s %7s %6s %9s %9s %9s %9s %10s %11s\n', ...
    'pair', 'scale', 'tol', '||b||_1', '||b||_2', '||M||_1', '||R||_1', ...
    'R bound', 'infeas rel');

for pairIdx = 1:2
    if pairIdx == 1
        X1 = helper.shift_profiles(w1, 3, L);   % true pair, 3-bin shift
        X2 = w1;
        pairName = 'true (shift 3)';
    else
        X1 = w1;                                % wrong pair, disjoint neurons
        X2 = w2;
        pairName = 'wrong (disjoint)';
    end

    for sc = scales
        A1 = X1 * sc;
        A2 = X2 * sc;
        % Exact lower bound on ||R||_1 implied by per-neuron mass conservation.
        Rbound = sum(abs(sum(A2, 2) - sum(A1, 2)));
        b = A2 - A1;

        for tol = tols
            opts = tfocs_SCD;
            opts.continuation = 1;
            opts.tol = tol;
            opts.stopCrit = 4;
            opts.maxIts = 500;
            opts.printEvery = 0;
            opts.alg = 'N83';
            copts = continuation();
            copts.verbose = 0;

            [~, M, R] = compute_EMD(A1, A2, opts, ...
                'continuationOptions', copts, 'lambdaR', lambdaR);

            infeas = norm(M * D' - R - b, 'fro') / max(norm(b, 'fro'), eps);

            fprintf('%-14s %7.4g %6.0e %9.4f %9.4f %9.4f %9.4f %10.4f %11.4f\n', ...
                pairName, sc, tol, norm(b(:), 1), norm(b(:), 2), ...
                norm(M(:), 1), norm(R(:), 1), Rbound, infeas);
        end
    end
    fprintf('\n');
end

%% Does the cost of a hopeless pair simply scale with mass?
fprintf('cost of the wrong pair against overall scale, tol=1e-4:\n');
opts = tfocs_SCD;
opts.continuation = 1;
opts.tol = 1e-4;
opts.stopCrit = 4;
opts.maxIts = 500;
opts.printEvery = 0;
opts.alg = 'N83';
copts = continuation();
copts.verbose = 0;

massScales = [0.5, 1, 2, 4, 8, 16];
costs = nan(size(massScales));
for k = 1:numel(massScales)
    % Only the estimate is rescaled, as happens when one factor explains less
    % of the data than the other.
    costs(k) = compute_EMD(w1, w2 * massScales(k), opts, ...
        'continuationOptions', copts, 'lambdaR', lambdaR);
    fprintf('  mass(estimate)=%5.1f  cost=%12.4f  cost/mass=%10.4f\n', ...
        massScales(k), costs(k), costs(k) / massScales(k));
end

%% Which lags in the sweep still contain the ground-truth motif at all?
% generate_data leaves the motif in the first columns of the L-wide window and
% pads the rest, so shift_profiles zero-fills the motif out of existence well
% before the sweep reaches its end. Those lags compare the estimate against an
% empty profile, which is cheaper than any real overlap, so they set the floor
% of the curve and the flanks that the retired dip-ratio classifier averaged
% over.
fprintf('\nprofile window L = %d\n', L);
for k = 1:K
    active = find(any(reshape(W(:, k, :), N, L), 1));
    fprintf('motif %d: mass lives in columns %d..%d\n', k, min(active), max(active));
end

fprintf('\n%-6s %s\n', 'lag', 'fraction of motif-1 mass still inside the window');
for s = -L:5:L
    kept = sum(helper.shift_profiles(w1, s, L), 'all') / sum(w1, 'all');
    fprintf('%-6d %.3f\n', s, kept);
end
