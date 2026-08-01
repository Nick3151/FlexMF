% BENCHMARK_LAG_RANGE
%
% How much runtime is there to gain in the lag sweep of similarity_WH_EMD?
%
% The sweep currently costs K*Khat*(2*max(L,Lhat)+1) calls to compute_EMD plus
% Khat calls for the loadings. Three ways to shrink that are compared:
%
%   1  restrict the lags to those that keep the whole motif inside the window,
%      which is also what removes the annihilation floor
%   2  search the lags coarse to fine instead of exhaustively, which is only
%      safe if the cost curve has a single minimum
%   3  both together
%
% Option 2 is only available if the curve shape is not needed, which is now the
% case: the dip-ratio classifier that needed the flanks has been retired.
%
% Per-solve times are measured at both scales and both tolerances, since the
% unit-mass normalisation changes the conditioning of the problem.

clear

root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(fullfile(root, 'Utils'))
addpath(genpath(fullfile(root, 'FlexMF')))

K = 3;
T = 2000;
[~, W, H, ~] = generate_data(T, 5 * ones(K, 1), 3 * ones(K, 1), 'noise', 0, ...
    'participation', 0.8 * ones(K, 1), 'seed', 1, 'len_burst', 1, 'dynamic', 0);
[N, ~, L] = size(W);
Khat = 5;
lambdaR = 1e3;

%% Lag counts under each scheme
support = find(any(reshape(W(:, 1, :), N, L), 1));
a = min(support);
b = max(support);

nLagsNow = 2 * max(L, L) + 1;
lagsKeep = (1 - a):(L - b);          % motif stays wholly inside the window
nLagsKeep = numel(lagsKeep);

fprintf('L = %d, motif support columns %d..%d\n', L, a, b);
fprintf('lags now                      : %d  (-%d:%d)\n', nLagsNow, L, L);
fprintf('lags keeping the motif intact : %d  (%d:%d)\n', ...
    nLagsKeep, lagsKeep(1), lagsKeep(end));

%% Is the cost curve unimodal enough for a coarse-to-fine search?
if isfile('emd_cost_vs_lag_unitmass_results.mat')
    S = load('emd_cost_vs_lag_unitmass_results.mat');
    res = S.res;
    fprintf('\ncoarse-to-fine argmin check on the unit-mass curves\n');
    for s = 1:numel(res)
        lags = res(s).lags;
        nLocalMin = 0;
        nRecovered = 0;
        nCurves = 0;
        excess = nan(1, res(s).K * res(s).Khat);
        for ii = 1:res(s).K
            for jj = 1:res(s).Khat
                c = squeeze(res(s).cost(ii, jj, :))';
                nCurves = nCurves + 1;
                nLocalMin = nLocalMin + ...
                    nnz([false, c(2:end-1) < c(1:end-2) & c(2:end-1) < c(3:end), false]);

                % Coarse pass on every second sample, then refine locally.
                coarseIdx = 1:2:numel(c);
                [~, ci] = min(c(coarseIdx));
                centre = coarseIdx(ci);
                win = max(1, centre - 2):min(numel(c), centre + 2);
                [~, ri] = min(c(win));
                if win(ri) == find(c == min(c), 1)
                    nRecovered = nRecovered + 1;
                end
                % Whether the argmin index is recovered matters less than how
                % much cost is given up by landing on a neighbouring lag.
                excess(nCurves) = (c(win(ri)) - min(c)) / max(min(c), eps);
            end
        end
        fprintf('  %-46s %d curves, %.1f local minima each, argmin %d/%d, excess cost median %.1f%% worst %.1f%%\n', ...
            res(s).name, nCurves, nLocalMin / nCurves, nRecovered, nCurves, ...
            100 * median(excess), 100 * max(excess));
    end
end

%% Per-solve timing
w1 = reshape(W(:, 1, :), N, L);
w2 = reshape(W(:, 2, :), N, L);
scaleCases = {'raw (mass ~20)', 20 / sum(w1(:)); 'unit mass', 1 / sum(w1(:))};
tols = [1e-4, 1e-6];

fprintf('\n%-16s %6s %14s %14s\n', 'scale', 'tol', 'W solve (s)', 'H solve (s)');
tW = nan(2, 2);
tH = nan(2, 2);
for sc = 1:2
    A1 = w1 * scaleCases{sc, 2};
    A2 = helper.shift_profiles(w1, 5, L) * scaleCases{sc, 2};
    A3 = w2 * scaleCases{sc, 2};
    h1 = H(1, :) / sum(H(1, :)) * scaleCases{sc, 2} * sum(w1(:));
    h2 = helper.shift_loadings(H(1, :), 5) / sum(H(1, :)) * scaleCases{sc, 2} * sum(w1(:));

    for ti = 1:2
        opts = tfocs_SCD;
        opts.continuation = 1;
        opts.tol = tols(ti);
        opts.stopCrit = 4;
        opts.maxIts = 500;
        opts.printEvery = 0;
        opts.alg = 'N83';
        copts = continuation();
        copts.verbose = 0;

        % Average a matching and a non-matching pair: the sweep does both.
        t0 = tic;
        compute_EMD(A1, A2, opts, 'continuationOptions', copts, 'lambdaR', lambdaR);
        compute_EMD(A1, A3, opts, 'continuationOptions', copts, 'lambdaR', lambdaR);
        tW(sc, ti) = toc(t0) / 2;

        t0 = tic;
        compute_EMD(h1, h2, opts, 'continuationOptions', copts, 'lambdaR', lambdaR);
        tH(sc, ti) = toc(t0);

        fprintf('%-16s %6.0e %14.3f %14.3f\n', ...
            scaleCases{sc, 1}, tols(ti), tW(sc, ti), tH(sc, ti));
    end
end

%% Projected cost of one similarity_WH_EMD call
fprintf('\nprojected wall time for one call, K=%d Khat=%d L=%d, unit mass, tol=1e-6\n', ...
    K, Khat, L);
tw = tW(2, 2);
th = tH(2, 2);
% Loadings are only solved for pairs that actually got matched, so at most
% min(K, Khat) of them, not Khat.
nH = min(K, Khat);
nCoarse = ceil(nLagsNow / 4) + 6;        % coarse step 4 plus a local refinement
nCoarseKeep = ceil(nLagsKeep / 4) + 6;

schemes = { ...
    'exhaustive, full lag range (today)', K * Khat * nLagsNow; ...
    'exhaustive, motif-preserving lags', K * Khat * nLagsKeep; ...
    'coarse-to-fine, full lag range',    K * Khat * nCoarse; ...
    'coarse-to-fine, motif-preserving',  K * Khat * nCoarseKeep};

base = schemes{1, 2} * tw + nH * th;
fprintf('  %-36s %6s %10s %8s %8s\n', '', 'W solves', 'W time', 'H time', 'speedup');
for s = 1:size(schemes, 1)
    total = schemes{s, 2} * tw + nH * th;
    fprintf('  %-36s %6d %8.1f m %6.1f m %7.1fx\n', ...
        schemes{s, 1}, schemes{s, 2}, schemes{s, 2} * tw / 60, ...
        nH * th / 60, base / total);
end
