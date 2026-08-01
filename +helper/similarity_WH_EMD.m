function [emds_W, emds_H, ids, details] = similarity_WH_EMD(W, H, W_hat, H_hat, varargin)
% SIMILARITY_WH_EMD  Match estimated motifs to ground-truth motifs by EMD.
%
%   [emds_W, emds_H, ids] = helper.similarity_WH_EMD(W, H, W_hat, H_hat)
%   [...] = helper.similarity_WH_EMD(..., 'Name', value)
%   [emds_W, emds_H, ids, details] = helper.similarity_WH_EMD(...)
%
% The pairwise EMD is computed between every ground-truth motif and every
% estimated motif at the lag that aligns them best, pairs that are too
% expensive to be a match are discarded, and the rest are claimed one-to-one,
% greedily from the cheapest pair. An estimate left with no unclaimed motif, or
% whose every pairing was discarded, is a failed match.
%
% Inputs
%   W      : N x K x L, or N x L for a single motif
%   H      : K x T
%   W_hat  : N x Khat x Lhat, or N x Lhat
%   H_hat  : Khat x T
%
% Options
%   LambdaR    price of unmatched mass inside compute_EMD (default 1e3)
%   Tol        opts.tol handed to the TFOCS solve (default 1e-6)
%   Rejection  'reference' to discard a pair costing more than RejectFrac of
%              what it costs to explain the estimate from an empty profile,
%              or 'none' to keep every pair and let greedy assign blindly
%              (default 'reference')
%   RejectFrac fraction used by 'reference' (default 0.7)
%   MaxShift   optional cap on |lag|, on top of the support-derived range
%   SupportTol fraction of profile mass allowed outside the lag window
%              (default 0.01)
%
% Outputs
%   emds_W : 1 x Khat, transport cost between each estimated profile and its
%            matched ground-truth profile, NaN where the match failed
%   emds_H : 1 x Khat, same for the loadings once the profile lag is undone
%   ids    : 1 x Khat, matched ground-truth index, 0 where the match failed
%   details: struct of per-pair diagnostics, see the end of this file
%
% SCALE. Every profile is scaled to unit total mass and every loading row to
% unit sum, so the result is invariant both to how a factorisation splits scale
% between W and H and to the overall amplitude of a factor. That second
% invariance matters: an estimated factor generally explains a different amount
% of mass than the motif it corresponds to, and lambdaR*||R||_1 charges for
% that difference, which on simulated data swamped the shape information
% entirely. With unit mass the reported transport cost reads directly as a mean
% displacement in bins.
%
% LAG RANGE. Shifting a profile zero-fills it, so a lag large enough to push
% the motif out of the window replaces it with an empty profile. Comparing an
% estimate against nothing is cheaper than overlapping it with an incompatible
% motif, so those lags used to form the floor of the cost curve and the minimum
% landed on them: the cost of a wrong pair became a function of the estimate
% alone, and the cheapest estimate won the greedy auction regardless of shape.
% The lag range is therefore derived from each motif's own support so that all
% but SupportTol of its mass stays inside the window.
%
% REJECTION. Greedy assignment decides which pair to take, never whether a
% pair is a match at all, so without a rejection rule this returns
% min(K, Khat) matches however poor they are. A curve-shape classifier used to
% do that job, comparing the middle of the cost curve against its flanks, but
% those flanks were exactly the lags that erase the motif, and on simulated
% data it rejected 1 of 18 wrong pairs while discarding 1 of 6 true ones.
%
% The 'reference' rule is the same idea measured against a well-defined
% baseline: one extra solve per estimate against an empty profile gives the
% cost of simply discarding it, and a pairing has to come in below RejectFrac
% of that. Since both profiles carry unit mass, a pairing that salvages
% nothing costs about twice the reference, so the ratio is roughly twice the
% fraction of mass abandoned. Measured over the two simulated settings in
% validate_similarity_WH_EMD, true pairs scored 0 to 0.43 and wrong pairs 1.18
% to 2.00; the default sits in that gap, nearer the wrong pairs because a
% false rejection surfaces as NaN and reads downstream as total failure.
%
% compute_EMD returns ||M||_1 + lambdaR*||R||_1. Matching uses that whole
% objective, because the unmatched-mass term is what makes an incompatible
% pairing expensive. The reported EMDs are the transport term ||M||_1 on its
% own, since lambdaR multiplies the solver's residual as well as any real
% unmatched mass. The unmatched-mass part of each matched pair is in details.

p = inputParser;
p.FunctionName = 'helper.similarity_WH_EMD';
addParameter(p, 'LambdaR', 1e3, @(x) isscalar(x) && isnumeric(x) && x > 0);
addParameter(p, 'Tol', 1e-6, @(x) isscalar(x) && isnumeric(x) && x > 0);
addParameter(p, 'Rejection', 'reference', ...
    @(x) any(strcmpi(x, {'reference', 'none'})));
addParameter(p, 'RejectFrac', 0.7, @(x) isscalar(x) && isnumeric(x) && x > 0);
addParameter(p, 'MaxShift', [], ...
    @(x) isempty(x) || (isscalar(x) && isnumeric(x) && x >= 0));
addParameter(p, 'SupportTol', 0.01, ...
    @(x) isscalar(x) && isnumeric(x) && x >= 0 && x < 0.5);
parse(p, varargin{:});

lambdaR = p.Results.LambdaR;
rejectFrac = p.Results.RejectFrac;
maxShift = p.Results.MaxShift;
supportTol = p.Results.SupportTol;
useReference = strcmpi(p.Results.Rejection, 'reference');
massTol = 1e-3;

%% Shapes and validation
szW = size(W);
if numel(szW) == 2
    N = szW(1);
    K = 1;
    L = szW(2);
    W = reshape(W, [N, 1, L]);
elseif numel(szW) == 3
    N = szW(1);
    K = szW(2);
    L = szW(3);
else
    error('helper:similarity_WH_EMD:InvalidW', ...
        'W must be N x L or N x K x L.');
end

szWhat = size(W_hat);
if numel(szWhat) == 2
    Khat = 1;
    Lhat = szWhat(2);
    W_hat = reshape(W_hat, [szWhat(1), 1, Lhat]);
elseif numel(szWhat) == 3
    Khat = szWhat(2);
    Lhat = szWhat(3);
else
    error('helper:similarity_WH_EMD:InvalidWhat', ...
        'W_hat must be N x Lhat or N x Khat x Lhat.');
end

assert(szWhat(1) == N, 'helper:similarity_WH_EMD:NeuronMismatch', ...
    'W and W_hat must have the same number of neurons.');
assert(size(H, 1) == K, 'helper:similarity_WH_EMD:InvalidH', ...
    'H must have one row per ground-truth motif.');
assert(size(H_hat, 1) == Khat, 'helper:similarity_WH_EMD:InvalidHhat', ...
    'H_hat must have one row per estimated motif.');
assert(size(H, 2) == size(H_hat, 2), 'helper:similarity_WH_EMD:TimeMismatch', ...
    'H and H_hat must have the same number of time bins.');

%% Zero factors, found before normalising while the masses are comparable
massW = sum(W, [1, 3]);
massWhat = sum(W_hat, [1, 3]);
massH = sum(H, 2)';
massHhat = sum(H_hat, 2)';

isZero = (massW == 0) | (massW < massTol * max(massW)) | (massH == 0);
isZeroHat = (massWhat == 0) | (massWhat < massTol * max(massWhat)) | ...
    (massHhat == 0);

%% Unit total mass per profile, unit sum per loading row
H = H ./ (massH(:) + eps);
H_hat = H_hat ./ (massHhat(:) + eps);
W = W ./ (reshape(massW, [1, K, 1]) + eps);
W_hat = W_hat ./ (reshape(massWhat, [1, Khat, 1]) + eps);

%% EMD options
opts = tfocs_SCD;
opts.continuation = 1;
opts.tol = p.Results.Tol;
opts.stopCrit = 4;
opts.maxIts = 500;
opts.printEvery = 0;
opts.alg = 'N83';
continue_opts = continuation();
continue_opts.verbose = 0;

%% Cost of explaining each estimate from an empty profile
% This is the reference level a pairing has to beat. It depends on lambdaR, on
% the window width and on the solver settings, so it is measured rather than
% assumed.
refCost = nan(1, Khat);
if useReference
    empty = zeros(N, Lhat);
    for jj = 1:Khat
        if isZeroHat(jj)
            continue
        end
        refCost(jj) = compute_EMD(empty, reshape(W_hat(:, jj, :), N, Lhat), ...
            opts, 'continuationOptions', continue_opts, 'lambdaR', lambdaR);
    end
end

%% Pairwise EMD at the best lag
cost = nan(K, Khat);
transport = nan(K, Khat);
residual = nan(K, Khat);
bestLag = nan(K, Khat);
rejected = false(K, Khat);
lagRange = cell(K, 1);
curves = cell(K, Khat);

for ii = 1:K
    if isZero(ii)
        continue
    end
    wk = reshape(W(:, ii, :), N, L);

    % Lags that keep all but supportTol of this motif's mass inside a window
    % of width Lhat.
    [a, b] = mass_support(wk, supportTol);
    lo = 1 - a;
    hi = Lhat - b;
    if ~isempty(maxShift)
        lo = max(lo, -maxShift);
        hi = min(hi, maxShift);
    end
    if hi < lo
        lags = 0;   % motif wider than the window, nothing else is available
    else
        lags = lo:hi;
    end
    lagRange{ii} = lags;
    nLags = numel(lags);

    for jj = 1:Khat
        if isZeroHat(jj)
            continue
        end
        wk_hat = reshape(W_hat(:, jj, :), N, Lhat);

        costLag = nan(1, nLags);
        transportLag = nan(1, nLags);
        residualLag = nan(1, nLags);

        for li = 1:nLags
            wtmp = helper.shift_profiles(wk, lags(li), Lhat);
            [costLag(li), M, R] = compute_EMD(wtmp, wk_hat, opts, ...
                'continuationOptions', continue_opts, 'lambdaR', lambdaR);
            transportLag(li) = norm(M(:), 1);
            residualLag(li) = lambdaR * norm(R(:), 1);
        end

        curves{ii, jj} = costLag;

        % Recorded whether or not the pair survives, so that details carries
        % the margin by which a pair was rejected and RejectFrac can be tuned.
        [~, li] = min(costLag);
        cost(ii, jj) = costLag(li);
        transport(ii, jj) = transportLag(li);
        residual(ii, jj) = residualLag(li);
        bestLag(ii, jj) = lags(li);

        rejected(ii, jj) = useReference && ...
            costLag(li) > rejectFrac * refCost(jj);
    end
end

%% Costs the assignment is allowed to use
matchCost = cost;
matchCost(rejected) = NaN;

%% Claim motifs one-to-one, cheapest pair first
emds_W = nan(1, Khat);
emds_H = nan(1, Khat);
ids = zeros(1, Khat);
matchedLag = nan(1, Khat);
residual_W = nan(1, Khat);
residual_H = nan(1, Khat);

unclaimed = matchCost;
for step = 1:min(K, Khat)
    if ~any(isfinite(unclaimed(:)))
        break
    end

    [~, idx] = min(unclaimed(:));
    [ii, jj] = ind2sub([K, Khat], idx);

    ids(jj) = ii;
    matchedLag(jj) = bestLag(ii, jj);
    emds_W(jj) = transport(ii, jj);
    residual_W(jj) = residual(ii, jj);

    unclaimed(ii, :) = NaN;
    unclaimed(:, jj) = NaN;
end

%% Loadings, with the profile lag undone
% Moving a profile s bins later and its loadings s bins earlier leaves the
% reconstruction unchanged, so the matched lag is applied negated here.
for jj = find(ids > 0)
    Hk = helper.shift_loadings(H(ids(jj), :), -matchedLag(jj));
    [~, M, R] = compute_EMD(Hk, H_hat(jj, :), opts, ...
        'continuationOptions', continue_opts, 'lambdaR', lambdaR);
    emds_H(jj) = norm(M(:), 1);
    residual_H(jj) = lambdaR * norm(R(:), 1);
end

%% Diagnostics
if nargout > 3
    details.lambdaR = lambdaR;
    details.rejection = p.Results.Rejection;
    details.rejectFrac = rejectFrac;
    details.refCost = refCost;
    details.lagRange = lagRange;
    details.curves = curves;
    details.cost = cost;
    details.costRatio = cost ./ refCost;
    details.matchCost = matchCost;
    details.transport = transport;
    details.residual = residual;
    details.bestLag = bestLag;
    details.rejected = rejected;
    details.matchedLag = matchedLag;
    details.residual_W = residual_W;
    details.residual_H = residual_H;
    details.isZeroFactor = isZero;
    details.isZeroEstimate = isZeroHat;
end

end

% ========================================================================

function [a, b] = mass_support(profile, tol)
% Narrowest column range holding all but tol of the profile mass at each end.
% A hard nonzero test would be useless on estimated profiles, which are dense
% with small values and would report the whole window as support.
colMass = sum(profile, 1);
total = sum(colMass);
if total <= 0
    a = 1;
    b = numel(colMass);
    return
end

c = cumsum(colMass) / total;
a = find(c > tol, 1, 'first');
b = find(c >= 1 - tol, 1, 'first');
if isempty(a)
    a = 1;
end
if isempty(b)
    b = numel(colMass);
end
b = max(a, b);
end
