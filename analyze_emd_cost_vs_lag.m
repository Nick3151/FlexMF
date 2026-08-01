% ANALYZE_EMD_COST_VS_LAG
%
% Does helper.similarity_WH_EMD still need helper.emd_match_classifier now
% that matching is one-to-one greedy, or can the minimum cost be taken
% directly?
%
% "Take the minimum" and "reject a pair outright" are separate jobs. Greedy
% assignment already decides which pair to take; a rejection rule decides
% whether a pair should be matched at all. This script measures three
% candidate rejection signals on simulated data:
%
%   dip ratio        what the since-retired emd_match_classifier used: the
%                    minimum over the central lags divided by the mean over
%                    the flanks
%   minimum cost     the full UOT objective at the best lag
%   residual share   lambdaR*||R||_1 / cost at the best lag, the fraction of
%                    the cost that comes from discarded rather than
%                    transported mass
%
% Two settings are run:
%   A  seqNMF on data from generate_data with Khat > K, as in
%      FlexMF_test_init, so two estimated factors have no true counterpart
%   B  a controlled set of estimates built from the ground truth: one faithful
%      copy of motif 1, one jittered duplicate of motif 1, and one pure noise
%      factor.  Motifs 2 and 3 exist in the truth but nothing resembles them,
%      so a matcher without a rejection rule is forced into two wrong pairs.
%
% The lag sweep uses opts.tol = 1e-4 rather than the 1e-6 that
% similarity_WH_EMD now uses. That costs about 20% accuracy on the smallest
% values and saves most of the runtime; the quantities compared here differ by
% orders of magnitude, so it does not affect any conclusion.

clear
close all

root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(fullfile(root, 'Utils'))
addpath(genpath(fullfile(root, 'FlexMF')))

lambdaR = 1e3;
cosineThreshold = 0.8;   % above this, call the pair a true correspondence

opts = tfocs_SCD;
opts.continuation = 1;
opts.tol = 1e-4;
opts.stopCrit = 4;
opts.maxIts = 500;
opts.printEvery = 0;
opts.alg = 'N83';
continue_opts = continuation();
continue_opts.verbose = 0;

%% ---------------------------------------------------------------- setting A
% Same generator and normalisation as FlexMF_test_init.
K = 3;
T = 2000;
Nneurons = 5 * ones(K, 1);
Dt = 3 * ones(K, 1);
participation = 0.8 * ones(K, 1);
seed = 1;

[X, W, H, ~] = generate_data(T, Nneurons, Dt, 'noise', 0, ...
    'participation', participation, 'seed', seed, ...
    'len_burst', 1, 'dynamic', 0);
L = size(W, 3);

Khat = 5;
frob = norm(X(:));
X = X / frob * Khat;
W = W / frob * Khat;

fprintf('setting A: N=%d K=%d L=%d T=%d Khat=%d\n', size(W, 1), K, L, T, Khat);

rng(0);
[What, Hhat] = seqNMF(X, 'K', Khat, 'L', L, 'lambda', 0.005, ...
    'maxiter', 100, 'showPlot', 0);

resA = sweep_pairs(W, H, What, Hhat, opts, continue_opts, lambdaR);
resA.name = 'A: seqNMF, Khat=5 > K=3';

%% ---------------------------------------------------------------- setting B
% Faithful copy of motif 1, jittered duplicate of motif 1, noise factor.
WhatB = zeros(size(W, 1), 3, L);
HhatB = zeros(3, T);

WhatB(:, 1, :) = W(:, 1, :);
HhatB(1, :) = H(1, :);

WhatB(:, 2, :) = helper.shift_profiles(reshape(W(:, 1, :), size(W, 1), L), 1, L);
HhatB(2, :) = helper.shift_loadings(H(1, :), -1);

rng(1);
noiseProfile = rand(size(W, 1), L) .* (rand(size(W, 1), L) > 0.7);
WhatB(:, 3, :) = noiseProfile * mean(sum(W, [1 3])) / max(sum(noiseProfile(:)), eps);
HhatB(3, :) = double(rand(1, T) > 0.98);

resB = sweep_pairs(W, H, WhatB, HhatB, opts, continue_opts, lambdaR);
resB.name = 'B: copy, duplicate, noise (Khat=3 = K)';

%% ------------------------------------------------------------------- report
print_table(resA);
print_table(resB);

plot_curves(resA, 'emd_cost_vs_lag_seqNMF.png');
plot_curves(resB, 'emd_cost_vs_lag_controlled.png');
plot_separation(resA, resB, cosineThreshold, 'emd_rejection_signals.png');

summarize(resA, resB, cosineThreshold);

save('emd_cost_vs_lag_results.mat', 'resA', 'resB', 'cosineThreshold');

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

function res = sweep_pairs(W, H, W_hat, H_hat, opts, continue_opts, lambdaR)
% Replicate the normalisation and lag sweep of similarity_WH_EMD, keeping the
% whole cost curve for every pair instead of only its minimum.

[N, K, L] = size(W);
[~, Khat, Lhat] = size(W_hat);

rowMass = sum(H, 2)';
rowMassHat = sum(H_hat, 2)';
Hn = H ./ (rowMass(:) + eps);
Wn = W .* reshape(rowMass, [1, K, 1]);
Wn_hat = W_hat .* reshape(rowMassHat, [1, Khat, 1]);

maxLag = max(L, Lhat);
lags = -maxLag:maxLag;
nLags = numel(lags);

cost = nan(K, Khat, nLags);
transport = nan(K, Khat, nLags);
residual = nan(K, Khat, nLags);

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
            residual(ii, jj, li) = lambdaR * norm(R(:), 1);
        end
    end
end

% Independent labels: shift-tolerant cosine similarity, no EMD involved.
[~, ~, cosIds, cosDetails] = helper.similarity_WH(W, Hn, W_hat, H_hat);

minCost = nan(K, Khat);
bestLag = nan(K, Khat);
residualShare = nan(K, Khat);
dipRatio = nan(K, Khat);
isMatch = false(K, Khat);

for ii = 1:K
    for jj = 1:Khat
        curve = squeeze(cost(ii, jj, :))';
        [minCost(ii, jj), li] = min(curve);
        bestLag(ii, jj) = lags(li);
        residualShare(ii, jj) = residual(ii, jj, li) / max(curve(li), eps);
        [isMatch(ii, jj), dipRatio(ii, jj)] = dip_ratio(curve);
    end
end

res = struct();
res.lags = lags;
res.cost = cost;
res.transport = transport;
res.residual = residual;
res.minCost = minCost;
res.bestLag = bestLag;
res.residualShare = residualShare;
res.dipRatio = dipRatio;
res.isMatch = isMatch;
res.cosine = cosDetails.S_pair;
res.cosIds = cosIds;
res.K = K;
res.Khat = Khat;
end

function print_table(res)
fprintf('\n=== %s ===\n', res.name);
fprintf('%4s %4s %8s %12s %12s %10s %9s %8s\n', ...
    'i', 'j', 'cosine', 'min cost', 'transport', 'resid %', 'dip ratio', 'classif');
for ii = 1:res.K
    for jj = 1:res.Khat
        li = find(res.lags == res.bestLag(ii, jj), 1);
        fprintf('%4d %4d %8.3f %12.4f %12.4f %10.1f %9.3f %8d\n', ...
            ii, jj, res.cosine(ii, jj), res.minCost(ii, jj), ...
            res.transport(ii, jj, li), 100 * res.residualShare(ii, jj), ...
            res.dipRatio(ii, jj), res.isMatch(ii, jj));
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
        isTrue = res.cosine(ii, jj) > 0.8;
        if isTrue
            plot(res.lags, curve, 'LineWidth', 2.5, 'Color', [0.85 0.2 0.15], ...
                'DisplayName', sprintf('est %d (cos %.2f) TRUE', jj, res.cosine(ii, jj)));
        else
            plot(res.lags, curve, 'LineWidth', 1.2, 'Color', [0.5 0.55 0.6 0.9], ...
                'DisplayName', sprintf('est %d (cos %.2f)', jj, res.cosine(ii, jj)));
        end
    end
    set(gca, 'YScale', 'log');
    xlabel('lag (bins)');
    ylabel('EMD cost (full UOT objective)');
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

function plot_separation(resA, resB, cosineThreshold, fileName)
cos_all = [resA.cosine(:); resB.cosine(:)];
minCost_all = [resA.minCost(:); resB.minCost(:)];
share_all = [resA.residualShare(:); resB.residualShare(:)];
dip_all = [resA.dipRatio(:); resB.dipRatio(:)];
isTrue = cos_all > cosineThreshold;

fig = figure('Visible', 'off', 'Position', [100 100 1400 400], 'Color', 'w');

subplot(1, 3, 1);
plot_signal(cos_all, minCost_all, isTrue, 'minimum cost', true, []);
subplot(1, 3, 2);
plot_signal(cos_all, 100 * share_all, isTrue, 'residual share (%)', false, []);
subplot(1, 3, 3);
plot_signal(cos_all, dip_all, isTrue, 'dip ratio', false, 0.85);

sgtitle('Candidate rejection signals vs an independent cosine label');
print(fig, fileName, '-dpng', '-r120');
close(fig);
fprintf('wrote %s\n', fileName);
end

function plot_signal(cosVals, yVals, isTrue, label, useLog, hline)
hold on
scatter(cosVals(~isTrue), yVals(~isTrue), 45, [0.45 0.5 0.55], 'filled', ...
    'DisplayName', 'wrong pair');
scatter(cosVals(isTrue), yVals(isTrue), 70, [0.85 0.2 0.15], 'filled', ...
    'DisplayName', 'true pair');
if ~isempty(hline)
    yline(hline, '--', 'threshold 0.85', 'Color', [0.2 0.2 0.2]);
end
if useLog
    set(gca, 'YScale', 'log');
end
xlabel('cosine similarity (independent label)');
ylabel(label);
legend('Location', 'best', 'Box', 'off', 'FontSize', 8);
grid on
box off
end

function summarize(resA, resB, cosineThreshold)
cos_all = [resA.cosine(:); resB.cosine(:)];
minCost_all = [resA.minCost(:); resB.minCost(:)];
share_all = [resA.residualShare(:); resB.residualShare(:)];
dip_all = [resA.dipRatio(:); resB.dipRatio(:)];
match_all = [resA.isMatch(:); resB.isMatch(:)];
isTrue = cos_all > cosineThreshold;

fprintf('\n=== separation over %d pairs (%d true, %d wrong) ===\n', ...
    numel(cos_all), nnz(isTrue), nnz(~isTrue));

report_gap('minimum cost', minCost_all, isTrue);
report_gap('residual share', share_all, isTrue);
report_gap('dip ratio', dip_all, isTrue);

fprintf('\ndip-ratio classifier verdicts:\n');
fprintf('  true pairs accepted : %d / %d\n', nnz(match_all & isTrue), nnz(isTrue));
fprintf('  wrong pairs rejected: %d / %d\n', nnz(~match_all & ~isTrue), nnz(~isTrue));
end

function report_gap(label, vals, isTrue)
maxTrue = max(vals(isTrue));
minWrong = min(vals(~isTrue));
if maxTrue < minWrong
    verdict = sprintf('SEPARATES, gap %.4g to %.4g', maxTrue, minWrong);
else
    verdict = sprintf('OVERLAPS, true up to %.4g, wrong down to %.4g', ...
        maxTrue, minWrong);
end
fprintf('  %-16s true [%.4g, %.4g]  wrong [%.4g, %.4g]  -> %s\n', ...
    label, min(vals(isTrue)), maxTrue, ...
    min(vals(~isTrue)), max(vals(~isTrue)), verdict);
end
