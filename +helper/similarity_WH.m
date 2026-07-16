function [S_W, S_H, ids, details] = similarity_WH(W, H, W_hat, H_hat, opts)
% similarity_WH
%
% Match motifs across W and W_hat using a shift-tolerant cosine similarity,
% then evaluate the corresponding activations H after greedy motif matching.
%
% The similarity of W is always the shift-tolerant cosine similarity.
% The similarity of H is evaluated after applying the opposite of the W
% shift, using a metric selected by opts.hMetric:
%   'cosine' : cosine similarity between the shifted H and H_hat (default)
%   'f1'     : event-based F_beta score after greedy event matching
%
% Inputs:
%   W      : N x K x L or N x L
%   H      : K x T
%   W_hat  : N x Khat x Lhat or N x Lhat
%   H_hat  : Khat x That
%   opts   : optional struct
%
% Key opts:
%   opts.hMetric         : 'cosine' (default) or 'f1'
%   opts.maxShift        : maximum temporal shift for W matching
%                          default = max(L, Lhat) - 1
%   opts.eventThreshold  : scalar threshold for event detection ('f1' only)
%                          default = []
%                          If empty, threshold is set separately for each H
%                          as mean(H) + std(H).
%   opts.eventTol        : event matching tolerance in time bins ('f1' only)
%                          default = 3
%   opts.minEventDistance: minimum distance between detected events ('f1' only)
%                          default = 2
%   opts.beta            : beta for F_beta ('f1' only)
%                          default = 1
%                          beta > 1 penalizes missing events more strongly.
%
% Outputs:
%   S_W     : 1 x Khat vector, matched shift-tolerant W similarity for each
%             W_hat motif
%   S_H     : 1 x Khat vector, matched H similarity for each H_hat motif,
%             using the metric selected by opts.hMetric
%   ids     : 1 x Khat vector, ids(j) is the matched W/H index for W_hat/H_hat j
%   details : struct with additional matching and activation diagnostics,
%             including both cosine and event-based H metrics regardless of
%             the selected opts.hMetric

if nargin < 5
    opts = struct();
end

if ~isfield(opts, 'hMetric'),          opts.hMetric = 'cosine'; end
if ~isfield(opts, 'maxShift'),         opts.maxShift = []; end
if ~isfield(opts, 'eventThreshold'),   opts.eventThreshold = []; end
if ~isfield(opts, 'eventTol'),         opts.eventTol = 3; end
if ~isfield(opts, 'minEventDistance'), opts.minEventDistance = 2; end
if ~isfield(opts, 'beta'),             opts.beta = 1; end

opts.hMetric = validatestring(opts.hMetric, {'cosine', 'f1'}, ...
    'similarity_WH', 'opts.hMetric');

% ---------- reshape W ----------
szW = size(W);
if ndims(W) == 2
    N = szW(1);
    K = 1;
    L = szW(2);
    W = reshape(W, [N, 1, L]);
elseif ndims(W) == 3
    N = szW(1);
    K = szW(2);
    L = szW(3);
else
    error('W must be N x L or N x K x L.');
end

% ---------- reshape W_hat ----------
szWhat = size(W_hat);
assert(N == szWhat(1), 'W and W_hat should have the same N.');

if ndims(W_hat) == 2
    Khat = 1;
    Lhat = szWhat(2);
    W_hat = reshape(W_hat, [N, 1, Lhat]);
elseif ndims(W_hat) == 3
    Khat = szWhat(2);
    Lhat = szWhat(3);
else
    error('W_hat must be N x Lhat or N x Khat x Lhat.');
end

% ---------- check H dimensions ----------
assert(size(H, 1) == K, 'H must be K x T.');
assert(size(H_hat, 1) == Khat, 'H_hat must be Khat x That.');

That = size(H_hat, 2);

if isempty(opts.maxShift)
    maxShift = max(L, Lhat) - 1;
else
    maxShift = opts.maxShift;
end

% ============================================================
% Step 1. Pairwise shift-tolerant W similarity
% ============================================================

S_pair = zeros(K, Khat);
shift_pair = zeros(K, Khat);

for ii = 1:K
    wk = squeeze(W(:, ii, :));

    for jj = 1:Khat
        wk_hat = squeeze(W_hat(:, jj, :));

        bestS = -Inf;
        bestShift = 0;

        for s = -maxShift:maxShift
            wk_shifted = shift_template_zero(wk, Lhat, s);
            Stmp = cosine_sim(wk_shifted(:), wk_hat(:));

            if Stmp > bestS
                bestS = Stmp;
                bestShift = s;
            end
        end

        S_pair(ii, jj) = bestS;
        shift_pair(ii, jj) = bestShift;
    end
end

% ============================================================
% Step 2. Greedy motif matching based on S_W
% ============================================================

S_W = zeros(1, Khat);
S_H = zeros(1, Khat);
ids = zeros(1, Khat);

precision = nan(1, Khat);
recall = nan(1, Khat);
Fbeta = nan(1, Khat);
medianLag = nan(1, Khat);
nRefEvents = zeros(1, Khat);
nHatEvents = zeros(1, Khat);
nMatchedEvents = zeros(1, Khat);
matchedShift = zeros(1, Khat);
H_cosine_shifted = nan(1, Khat);
details.eventInfo = cell(1, Khat);

temp = S_pair;

for mm = 1:min(K, Khat)

    if ~any(isfinite(temp(:))) || max(temp(:)) == -Inf
        break
    end

    [~, linearIdx] = max(temp(:));
    [r, c] = ind2sub(size(temp), linearIdx);

    ids(c) = r;
    S_W(c) = S_pair(r, c);
    matchedShift(c) = shift_pair(r, c);

    % Shift H in the opposite direction of the W shift.
    Hk = H(r, :);
    Hk_hat = H_hat(c, :);

    Hk_shifted = shift_vector_zero(Hk, That, -matchedShift(c));

    % Cosine similarity of activations.
    H_cosine_shifted(c) = cosine_sim(Hk_shifted(:), Hk_hat(:));

    % Event-based activation similarity.
    [metrics, eventInfo] = activation_event_similarity( ...
        Hk_shifted, ...
        Hk_hat, ...
        opts.eventThreshold, ...
        opts.eventTol, ...
        opts.minEventDistance, ...
        opts.beta);

    precision(c) = metrics.precision;
    recall(c) = metrics.recall;
    Fbeta(c) = metrics.Fbeta;
    medianLag(c) = metrics.medianLag;
    nRefEvents(c) = metrics.nRefEvents;
    nHatEvents(c) = metrics.nHatEvents;
    nMatchedEvents(c) = metrics.nMatchedEvents;

    details.eventInfo{c} = eventInfo;

    % Selected H similarity metric.
    switch opts.hMetric
        case 'cosine'
            S_H(c) = H_cosine_shifted(c);
        case 'f1'
            S_H(c) = metrics.Fbeta;
    end

    % Remove matched row and column.
    temp(r, :) = -Inf;
    temp(:, c) = -Inf;
end

% ============================================================
% Output details
% ============================================================

details.hMetric = opts.hMetric;
details.S_pair = S_pair;
details.shift_pair = shift_pair;
details.matchedShift = matchedShift;
details.H_cosine_shifted = H_cosine_shifted;

details.precision = precision;
details.recall = recall;
details.Fbeta = Fbeta;
details.medianLag = medianLag;
details.nRefEvents = nRefEvents;
details.nHatEvents = nHatEvents;
details.nMatchedEvents = nMatchedEvents;

details.opts = opts;

end

% ========================================================================
% Helper: zero-padded template shift
% ========================================================================

function Xout = shift_template_zero(X, targetL, shift)
% X is N x L.
% Output is N x targetL.
% Positive shift moves X to the right.

[N, L] = size(X);
Xout = zeros(N, targetL);

for t = 1:L
    tt = t + shift;
    if tt >= 1 && tt <= targetL
        Xout(:, tt) = X(:, t);
    end
end

end

% ========================================================================
% Helper: zero-padded vector shift
% ========================================================================

function yout = shift_vector_zero(y, targetT, shift)
% y is 1 x T.
% Output is 1 x targetT.
% Positive shift moves y to the right.

y = y(:)';
T = length(y);
yout = zeros(1, targetT);

for t = 1:T
    tt = t + shift;
    if tt >= 1 && tt <= targetT
        yout(tt) = y(t);
    end
end

end

% ========================================================================
% Helper: cosine similarity
% ========================================================================

function s = cosine_sim(a, b)
a = a(:);
b = b(:);

denom = sqrt(sum(a.^2)) * sqrt(sum(b.^2)) + eps;
s = sum(a .* b) / denom;

if isnan(s)
    s = 0;
end

end

% ========================================================================
% Helper: event-based activation similarity
% ========================================================================

function [metrics, eventInfo] = activation_event_similarity(hRef, hHat, threshold, eventTol, minDist, beta)

hRef = hRef(:)';
hHat = hHat(:)';

if isempty(threshold)
    thrRef = mean(hRef) + std(hRef);
    thrHat = mean(hHat) + std(hHat);
else
    thrRef = threshold;
    thrHat = threshold;
end

[tRef, aRef] = detect_events_localmax(hRef, thrRef, minDist);
[tHat, aHat] = detect_events_localmax(hHat, thrHat, minDist);

matches = match_events_one_to_one(tRef, tHat, eventTol);

nRef = length(tRef);
nHat = length(tHat);
nMatch = size(matches, 1);

if nMatch > 0
    refMatchedIdx = matches(:, 1);
    hatMatchedIdx = matches(:, 2);

    refMatchedAmp = sum(aRef(refMatchedIdx));
    hatMatchedAmp = sum(aHat(hatMatchedIdx));

    totalRefAmp = sum(aRef) + eps;
    totalHatAmp = sum(aHat) + eps;

    weightedRecall = refMatchedAmp / totalRefAmp;
    weightedPrecision = hatMatchedAmp / totalHatAmp;

    lags = tHat(hatMatchedIdx) - tRef(refMatchedIdx);
    medLag = median(lags);
else
    weightedRecall = double(nRef == 0);
    weightedPrecision = double(nHat == 0);
    medLag = NaN;
    lags = [];
end

if weightedPrecision == 0 && weightedRecall == 0
    Fbeta = 0;
else
    Fbeta = (1 + beta^2) * weightedPrecision * weightedRecall / ...
            (beta^2 * weightedPrecision + weightedRecall + eps);
end

metrics.precision = weightedPrecision;
metrics.recall = weightedRecall;
metrics.Fbeta = Fbeta;
metrics.medianLag = medLag;
metrics.nRefEvents = nRef;
metrics.nHatEvents = nHat;
metrics.nMatchedEvents = nMatch;

eventInfo.tRef = tRef;
eventInfo.aRef = aRef;
eventInfo.tHat = tHat;
eventInfo.aHat = aHat;
eventInfo.matches = matches;
eventInfo.lags = lags;
eventInfo.thresholdRef = thrRef;
eventInfo.thresholdHat = thrHat;

end

% ========================================================================
% Helper: detect local-max events without requiring findpeaks()
% ========================================================================

function [tEvents, aEvents] = detect_events_localmax(x, threshold, minDist)

x = x(:)';
T = length(x);

candidateTimes = [];
candidateAmps = [];

for t = 1:T
    leftOK = (t == 1) || (x(t) >= x(t - 1));
    rightOK = (t == T) || (x(t) > x(t + 1));

    if x(t) >= threshold && leftOK && rightOK
        candidateTimes(end + 1) = t; %#ok<AGROW>
        candidateAmps(end + 1) = x(t); %#ok<AGROW>
    end
end

if isempty(candidateTimes)
    tEvents = [];
    aEvents = [];
    return;
end

% Enforce minimum event distance by keeping larger peaks first.
[~, order] = sort(candidateAmps, 'descend');

kept = false(size(candidateTimes));
blocked = false(size(candidateTimes));

for ii = 1:length(order)
    idx = order(ii);

    if blocked(idx)
        continue;
    end

    kept(idx) = true;

    tooClose = abs(candidateTimes - candidateTimes(idx)) < minDist;
    blocked = blocked | tooClose;
    blocked(idx) = false;
end

tEvents = candidateTimes(kept);
aEvents = candidateAmps(kept);

[tEvents, sortIdx] = sort(tEvents);
aEvents = aEvents(sortIdx);

end

% ========================================================================
% Helper: one-to-one event matching within tolerance
% ========================================================================

function matches = match_events_one_to_one(tRef, tHat, eventTol)

matches = [];

if isempty(tRef) || isempty(tHat)
    return;
end

candidates = [];

for i = 1:length(tRef)
    for j = 1:length(tHat)
        d = abs(tRef(i) - tHat(j));
        if d <= eventTol
            candidates(end + 1, :) = [i, j, d]; %#ok<AGROW>
        end
    end
end

if isempty(candidates)
    return;
end

% Greedy nearest-time matching.
[~, order] = sort(candidates(:, 3), 'ascend');
candidates = candidates(order, :);

usedRef = false(1, length(tRef));
usedHat = false(1, length(tHat));

for cc = 1:size(candidates, 1)
    i = candidates(cc, 1);
    j = candidates(cc, 2);

    if ~usedRef(i) && ~usedHat(j)
        matches(end + 1, :) = [i, j]; %#ok<AGROW>
        usedRef(i) = true;
        usedHat(j) = true;
    end
end

end
