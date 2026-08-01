% ANALYZE_EMD_COST_VS_LAG_NORMALIZED
%
% Follow-up to analyze_emd_cost_vs_lag. That script showed that on simulated
% data neither the dip ratio nor the raw UOT cost separates a true motif pair
% from a wrong one, and that the cost of a pair is nearly independent of which
% ground-truth motif is involved. The suspicion is that lambdaR*||R||_1 is
% paying for a mismatch in total mass between a ground-truth factor and an
% estimated one, which swamps the shape information.
%
% This script tests that. Every profile is scaled to unit total mass before the
% EMD, so no cost can come from a global scale difference and ||R||_1 is
% bounded by 2. lambdaR is then the price of abandoning a unit of mass measured
% in bins of displacement, so it should be of the order of the profile width,
% not 1e3. Both values are swept.
%
% Under this normalisation the cost is comparable across pairs: it is the mean
% displacement per unit of mass, saturating at 2*lambdaR when nothing can be
% matched. The question is whether the minimum over lags then separates true
% pairs from wrong ones without any curve-shape test.
%
% Lags are subsampled by 2 to keep the runtime near half an hour; the curves
% are smooth on that scale.

clear
close all

root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(fullfile(root, 'Utils'))
addpath(genpath(fullfile(root, 'FlexMF')))

lambdaRs = [70, 1000];
lagStep = 2;
cosineThreshold = 0.8;

opts = tfocs_SCD;
opts.continuation = 1;
opts.tol = 1e-4;
opts.stopCrit = 4;
opts.maxIts = 500;
opts.printEvery = 0;
opts.alg = 'N83';
continue_opts = continuation();
continue_opts.verbose = 0;

%% Fixtures, identical to analyze_emd_cost_vs_lag
K = 3;
T = 2000;
Nneurons = 5 * ones(K, 1);
Dt = 3 * ones(K, 1);
participation = 0.8 * ones(K, 1);

[X, W, H, ~] = generate_data(T, Nneurons, Dt, 'noise', 0, ...
    'participation', participation, 'seed', 1, ...
    'len_burst', 1, 'dynamic', 0);
L = size(W, 3);
N = size(W, 1);

Khat = 5;
frob = norm(X(:));
X = X / frob * Khat;
W = W / frob * Khat;

rng(0);
[WhatA, HhatA] = seqNMF(X, 'K', Khat, 'L', L, 'lambda', 0.005, ...
    'maxiter', 100, 'showPlot', 0);

WhatB = zeros(N, 3, L);
HhatB = zeros(3, T);
WhatB(:, 1, :) = W(:, 1, :);
HhatB(1, :) = H(1, :);
WhatB(:, 2, :) = helper.shift_profiles(reshape(W(:, 1, :), N, L), 1, L);
HhatB(2, :) = helper.shift_loadings(H(1, :), -1);
rng(1);
noiseProfile = rand(N, L) .* (rand(N, L) > 0.7);
WhatB(:, 3, :) = noiseProfile * mean(sum(W, [1 3])) / max(sum(noiseProfile(:)), eps);
HhatB(3, :) = double(rand(1, T) > 0.98);

%% Where does the current cost come from?
report_masses('A: seqNMF, Khat=5 > K=3', W, H, WhatA, HhatA);
report_masses('B: copy, duplicate, noise', W, H, WhatB, HhatB);

%% Mass-normalised sweeps
res = struct([]);
for r = 1:numel(lambdaRs)
    for s = 1:2
        if s == 1
            [Wh, Hh, tag] = deal(WhatA, HhatA, 'A: seqNMF, Khat=5 > K=3');
        else
            [Wh, Hh, tag] = deal(WhatB, HhatB, 'B: copy, duplicate, noise');
        end
        out = sweep_pairs(W, H, Wh, Hh, opts, continue_opts, lambdaRs(r), lagStep);
        out.name = sprintf('%s  |  lambdaR=%g, unit mass', tag, lambdaRs(r));
        out.lambdaR = lambdaRs(r);
        out.setting = s;
        print_table(out);
        if isempty(res)
            res = out;
        else
            res(end + 1) = out; %#ok<SAGROW>
        end
    end
end

%% Plots and verdict
for r = 1:numel(lambdaRs)
    sel = find([res.lambdaR] == lambdaRs(r));
    for s = sel
        plot_curves(res(s), sprintf('emd_cost_vs_lag_unitmass_lambdaR%g_%s.png', ...
            res(s).lambdaR, char('A' + res(s).setting - 1)));
    end
end

plot_separation(res, cosineThreshold, 'emd_rejection_signals_unitmass.png');

for r = 1:numel(lambdaRs)
    sel = [res.lambdaR] == lambdaRs(r);
    summarize(res(sel), cosineThreshold, lambdaRs(r));
end

save('emd_cost_vs_lag_unitmass_results.mat', 'res', 'cosineThreshold', 'lambdaRs');

%% ===================================================================== local
function [isMatch, ratio] = dip_ratio(curve, threshold)
% The retired helper.emd_match_classifier, kept here as a local function so
% that this analysis, which is the evidence for retiring it, still runs. The
% minimum over the middle 40% of the curve divided by the mean over the flanks,
% called a match when that falls below the threshold.
if nargin < 2 || isempty(threshold)
    threshold = 0.85;
end

curve = curve(:);
n = numel(curve);
mid = ceil(n / 2);
halfw = max(1, round(n * 0.4 / 2));
cStart = max(1, mid - halfw);
cEnd = min(n, mid + halfw);

flanks = curve([1:cStart - 1, cEnd + 1:n]);
if isempty(flanks) || mean(flanks) == 0
    ratio = 1;
    isMatch = false;
    return
end

ratio = min(curve(cStart:cEnd)) / mean(flanks);
isMatch = ratio < threshold;
end

function report_masses(name, W, H, W_hat, H_hat)
% Total mass of each factor the way similarity_WH_EMD currently sees it, and
% the pairwise mass mismatch that lambdaR*||R||_1 has to pay for.
rowMass = sum(H, 2)';
rowMassHat = sum(H_hat, 2)';
massW = sum(W, [1, 3]) .* rowMass;
massWhat = sum(W_hat, [1, 3]) .* rowMassHat;

fprintf('\n=== factor masses, %s ===\n', name);
fprintf('ground truth : %s\n', num2str(massW, '%12.4f'));
fprintf('estimates    : %s\n', num2str(massWhat, '%12.4f'));
fprintf('|mass difference| per pair (lower bound on ||R||_1):\n');
for ii = 1:numel(massW)
    fprintf('   truth %d : %s\n', ii, ...
        num2str(abs(massW(ii) - massWhat), '%12.4f'));
end
end

function res = sweep_pairs(W, H, W_hat, H_hat, opts, continue_opts, lambdaR, lagStep)
[N, K, L] = size(W);
[~, Khat, Lhat] = size(W_hat);

% Unit total mass per profile and per loading row: matching becomes invariant
% to how a factorisation splits scale between W and H, and to the overall
% amplitude of an estimated factor.
Hn = H ./ (sum(H, 2) + eps);
Wn = W ./ (sum(W, [1, 3]) + eps);
Wn_hat = W_hat ./ (sum(W_hat, [1, 3]) + eps);

maxLag = max(L, Lhat);
lags = -maxLag:lagStep:maxLag;
nLags = numel(lags);

cost = nan(K, Khat, nLags);
transport = nan(K, Khat, nLags);
residualMass = nan(K, Khat, nLags);

for ii = 1:K
    wk = reshape(Wn(:, ii, :), N, L);
    for jj = 1:Khat
        wk_hat = reshape(Wn_hat(:, jj, :), N, Lhat);
        for li = 1:nLags
            wtmp = helper.shift_profiles(wk, lags(li), Lhat);
            [c, M, R] = compute_EMD(wtmp, wk_hat, opts, ...
                'continuationOptions', continue_opts, 'lambdaR', lambdaR);
            cost(ii, jj, li) = c;
            transport(ii, jj, li) = norm(M(:), 1);
            residualMass(ii, jj, li) = norm(R(:), 1);
        end
        fprintf('  lambdaR=%g pair (%d,%d) done, min cost %.4f\n', ...
            lambdaR, ii, jj, min(cost(ii, jj, :)));
    end
end

[~, ~, cosIds, cosDetails] = helper.similarity_WH(W, Hn, W_hat, H_hat);

minCost = nan(K, Khat);
bestLag = nan(K, Khat);
bestTransport = nan(K, Khat);
bestResidualMass = nan(K, Khat);
dipRatio = nan(K, Khat);
isMatch = false(K, Khat);
for ii = 1:K
    for jj = 1:Khat
        curve = squeeze(cost(ii, jj, :))';
        [minCost(ii, jj), li] = min(curve);
        bestLag(ii, jj) = lags(li);
        bestTransport(ii, jj) = transport(ii, jj, li);
        bestResidualMass(ii, jj) = residualMass(ii, jj, li);
        [isMatch(ii, jj), dipRatio(ii, jj)] = dip_ratio(curve);
    end
end

res = struct();
res.lags = lags;
res.cost = cost;
res.minCost = minCost;
res.bestLag = bestLag;
res.bestTransport = bestTransport;
res.bestResidualMass = bestResidualMass;
res.dipRatio = dipRatio;
res.isMatch = isMatch;
res.cosine = cosDetails.S_pair;
res.cosIds = cosIds;
res.K = K;
res.Khat = Khat;
res.name = '';
res.lambdaR = lambdaR;
res.setting = 0;
end

function print_table(res)
fprintf('\n=== %s ===\n', res.name);
fprintf('%4s %4s %8s %12s %12s %10s %8s %9s %8s\n', ...
    'i', 'j', 'cosine', 'min cost', 'transport', '||R||_1', 'lag', 'dip', 'classif');
for ii = 1:res.K
    for jj = 1:res.Khat
        fprintf('%4d %4d %8.3f %12.4f %12.4f %10.4f %8d %9.3f %8d\n', ...
            ii, jj, res.cosine(ii, jj), res.minCost(ii, jj), ...
            res.bestTransport(ii, jj), res.bestResidualMass(ii, jj), ...
            res.bestLag(ii, jj), res.dipRatio(ii, jj), res.isMatch(ii, jj));
    end
end
end

function plot_curves(res, fileName)
fig = figure('Visible', 'off', 'Position', [100 100 1500 420], 'Color', 'w');
for ii = 1:res.K
    subplot(1, res.K, ii);
    hold on
    for jj = 1:res.Khat
        curve = squeeze(res.cost(ii, jj, :));
        if res.cosine(ii, jj) > 0.8
            plot(res.lags, curve, 'LineWidth', 2.5, 'Color', [0.85 0.2 0.15], ...
                'DisplayName', sprintf('est %d (cos %.2f) TRUE', jj, res.cosine(ii, jj)));
        else
            plot(res.lags, curve, 'LineWidth', 1.2, 'Color', [0.5 0.55 0.6 0.9], ...
                'DisplayName', sprintf('est %d (cos %.2f)', jj, res.cosine(ii, jj)));
        end
    end
    yline(2 * res.lambdaR, ':', 'no mass matched', 'Color', [0.2 0.2 0.2]);
    set(gca, 'YScale', 'log');
    xlabel('lag (bins)');
    ylabel('EMD cost, unit mass');
    title(sprintf('ground-truth motif %d', ii));
    legend('Location', 'southoutside', 'Box', 'off', 'FontSize', 7);
    grid on
    box off
end
sgtitle(sprintf('Cost vs lag  --  %s', res.name));
print(fig, fileName, '-dpng', '-r120');
close(fig);
fprintf('wrote %s\n', fileName);
end

function plot_separation(res, cosineThreshold, fileName)
lambdaRs = unique([res.lambdaR]);
fig = figure('Visible', 'off', 'Position', [100 100 1000 420], 'Color', 'w');
for r = 1:numel(lambdaRs)
    sel = find([res.lambdaR] == lambdaRs(r));
    cosVals = [];
    costVals = [];
    for s = sel
        cosVals = [cosVals; res(s).cosine(:)]; %#ok<AGROW>
        costVals = [costVals; res(s).minCost(:)]; %#ok<AGROW>
    end
    isTrue = cosVals > cosineThreshold;

    subplot(1, numel(lambdaRs), r);
    hold on
    scatter(cosVals(~isTrue), costVals(~isTrue), 45, [0.45 0.5 0.55], 'filled', ...
        'DisplayName', 'wrong pair');
    scatter(cosVals(isTrue), costVals(isTrue), 70, [0.85 0.2 0.15], 'filled', ...
        'DisplayName', 'true pair');
    yline(2 * lambdaRs(r), ':', 'no mass matched', 'Color', [0.2 0.2 0.2]);
    yline(lambdaRs(r), '--', 'half the mass unmatched', 'Color', [0.2 0.2 0.2]);
    set(gca, 'YScale', 'log');
    xlabel('cosine similarity (independent label)');
    ylabel('minimum cost, unit mass');
    title(sprintf('lambdaR = %g', lambdaRs(r)));
    legend('Location', 'east', 'Box', 'off', 'FontSize', 8);
    grid on
    box off
end
sgtitle('Minimum cost after unit-mass normalisation');
print(fig, fileName, '-dpng', '-r120');
close(fig);
fprintf('wrote %s\n', fileName);
end

function summarize(res, cosineThreshold, lambdaR)
cosVals = [];
costVals = [];
dipVals = [];
for s = 1:numel(res)
    cosVals = [cosVals; res(s).cosine(:)]; %#ok<AGROW>
    costVals = [costVals; res(s).minCost(:)]; %#ok<AGROW>
    dipVals = [dipVals; res(s).dipRatio(:)]; %#ok<AGROW>
end
isTrue = cosVals > cosineThreshold;

fprintf('\n=== separation, lambdaR=%g, unit mass, %d pairs (%d true, %d wrong) ===\n', ...
    lambdaR, numel(cosVals), nnz(isTrue), nnz(~isTrue));
report_gap('minimum cost', costVals, isTrue);
report_gap('dip ratio', dipVals, isTrue);
fprintf('  a threshold at lambdaR=%g (half the mass unmatched) would give:\n', lambdaR);
fprintf('    true pairs kept     : %d / %d\n', nnz(isTrue & costVals < lambdaR), nnz(isTrue));
fprintf('    wrong pairs rejected: %d / %d\n', nnz(~isTrue & costVals >= lambdaR), nnz(~isTrue));
end

function report_gap(label, vals, isTrue)
maxTrue = max(vals(isTrue));
minWrong = min(vals(~isTrue));
if maxTrue < minWrong
    verdict = sprintf('SEPARATES, gap %.4g to %.4g (factor %.1f)', ...
        maxTrue, minWrong, minWrong / max(maxTrue, eps));
else
    verdict = sprintf('OVERLAPS, true up to %.4g, wrong down to %.4g', maxTrue, minWrong);
end
fprintf('  %-14s true [%.4g, %.4g]  wrong [%.4g, %.4g]  -> %s\n', ...
    label, min(vals(isTrue)), maxTrue, min(vals(~isTrue)), max(vals(~isTrue)), verdict);
end
