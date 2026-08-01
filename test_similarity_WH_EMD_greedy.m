function tests = test_similarity_WH_EMD_greedy
% Specification tests for one-to-one greedy motif matching in
% helper.similarity_WH_EMD.
%
% These tests describe the intended behaviour rather than the current
% behaviour, so failures and errors are expected until the function is
% updated.  The specification under test is:
%   * motifs are matched one-to-one, greedily on the pairwise EMD matrix,
%   * a failed match returns an EMD of NaN and an id of 0,
%   * two identical motifs matched together return an EMD of exactly 0.
%
% Every expected EMD is derived from the fixture construction.  The
% derivation rests on two properties of compute_EMD:
%
%   1. d = ||M||_1 + lambdaR*||R||_1, and the divergence operator M*D' sums
%      to zero along each neuron, so ||R||_1 is at least the per-neuron mass
%      imbalance between the two inputs.  A perturbation that *moves* mass
%      leaves every row mass unchanged, so R = 0 is optimal and d is exactly
%      the transport cost, amplitude times distance in bins.
%   2. compute_EMD short-circuits to d = 0 when its two inputs are equal, so
%      an exact match is exactly 0 and needs no tolerance at all.
%
% The only number here that is not derived from a fixture is
% solver_rel_tol, which covers the accuracy of the TFOCS solve itself.
%
% Run with:
%   results = runtests('test_similarity_WH_EMD_greedy.m');
%   table(results)

tests = functiontests(localfunctions);
end

function setupOnce(testCase)
import matlab.unittest.fixtures.PathFixture

repoDir = fileparts(mfilename('fullpath'));
tfocsDir = fullfile(fileparts(repoDir), 'TFOCS');

assert(isfolder(tfocsDir), ...
    'Expected a TFOCS folder next to the FlexMF repository, found none.');

% PathFixture restores the original path on teardown, including the case
% where these folders were already on it.
testCase.applyFixture(PathFixture(repoDir));
testCase.applyFixture(PathFixture(tfocsDir));
end

% ========================================================================
% Exact matches
% ========================================================================

function testIdenticalMotifsGiveZeroEMD(testCase)
% A permutation of the ground truth is bitwise identical after the internal
% normalisation, so compute_EMD short-circuits and both EMDs are exactly 0.
fx = disjoint_fixture();
permutation = [3, 1, 2];

W_hat = fx.W(:, permutation, :);
H_hat = fx.H(permutation, :);

[emdW, emdH, ids] = helper.similarity_WH_EMD(fx.W, fx.H, W_hat, H_hat);

verifyEqual(testCase, ids, permutation);
verifyEqual(testCase, emdW, zeros(1, numel(permutation)));
verifyEqual(testCase, emdH, zeros(1, numel(permutation)));
end

function testSingleMotifPerfectMatch(testCase)
% A 2-D W with a single motif compared against itself.  A valid zero-cost
% match must not be confused with an absent match.
fx = disjoint_fixture();

W_single = reshape(fx.W(:, 1, :), size(fx.W, 1), size(fx.W, 3));
H_single = fx.H(1, :);

[emdW, emdH, ids] = ...
    helper.similarity_WH_EMD(W_single, H_single, W_single, H_single);

verifyEqual(testCase, ids, 1);
verifyEqual(testCase, emdW, 0);
verifyEqual(testCase, emdH, 0);
end

function testTemporalShiftIsUndoneInLoadings(testCase)
% Moving a profile s bins later and its loadings s bins earlier leaves the
% reconstruction W*H unchanged.  The lag search must find s, and the
% opposite shift applied to H must line the loadings back up, so both EMDs
% come back as exactly 0.  Nothing else in this file exercises the sign of
% that opposite shift.
fx = disjoint_fixture();
s = 2;
permutation = [3, 1, 2];

W_hat = helper.shift_profiles(fx.W(:, permutation, :), s);
H_hat = helper.shift_loadings(fx.H(permutation, :), -s);

[emdW, emdH, ids] = helper.similarity_WH_EMD(fx.W, fx.H, W_hat, H_hat);

verifyEqual(testCase, ids, permutation);
verifyEqual(testCase, emdW, zeros(1, numel(permutation)));
verifyEqual(testCase, emdH, zeros(1, numel(permutation)));
end

% ========================================================================
% Matches with a residual
% ========================================================================

function testMatchedPairResidualEqualsTransportCost(testCase)
% A residual of 10% of peak amplitude, displaced by a known number of bins,
% must come back as exactly that transport cost.
fx = disjoint_fixture();
permutation = [2, 3, 1];

[W_hat, H_hat, expectedW, expectedH] = ...
    perturb_estimates(fx.W(:, permutation, :), fx.H(permutation, :));

[emdW, emdH, ids] = helper.similarity_WH_EMD(fx.W, fx.H, W_hat, H_hat);

verifyEqual(testCase, ids, permutation);
verifyEqual(testCase, emdW, expectedW, 'RelTol', solver_rel_tol());
verifyEqual(testCase, emdH, expectedH, 'RelTol', solver_rel_tol());
end

function testFactorizationScalingIsInvariant(testCase)
% W(:,k,:)*c together with H(k,:)/c is the same factorisation, and
% similarity_WH_EMD normalises the profile total and the loading row sum
% separately, so the answer must not move.  Powers of two keep the rescaling
% exact in binary, leaving the solver as the only source of difference.
fx = disjoint_fixture();
permutation = [2, 3, 1];

[W_hat, H_hat, expectedW, expectedH] = ...
    perturb_estimates(fx.W(:, permutation, :), fx.H(permutation, :));
[W_scaled, H_scaled] = rescale_factors(W_hat, H_hat, [0.25, 4, 8]);

[emdW, emdH, ids] = ...
    helper.similarity_WH_EMD(fx.W, fx.H, W_scaled, H_scaled);

verifyEqual(testCase, ids, permutation);
verifyEqual(testCase, emdW, expectedW, 'RelTol', solver_rel_tol());
verifyEqual(testCase, emdH, expectedH, 'RelTol', solver_rel_tol());
end

% ========================================================================
% Unequal motif counts
% ========================================================================

function testFewerEstimatesAreAllMatched(testCase)
% Khat < K.  Every estimate should find its own ground-truth motif.
fx = disjoint_fixture();
permutation = [3, 1];

[W_hat, H_hat, expectedW, expectedH] = ...
    perturb_estimates(fx.W(:, permutation, :), fx.H(permutation, :));

[emdW, emdH, ids] = helper.similarity_WH_EMD(fx.W, fx.H, W_hat, H_hat);

verifyEqual(testCase, ids, permutation);
verifyEqual(testCase, numel(unique(ids)), numel(ids));
verifyEqual(testCase, emdW, expectedW, 'RelTol', solver_rel_tol());
verifyEqual(testCase, emdH, expectedH, 'RelTol', solver_rel_tol());
end

function testExtraEstimatesReportFailedMatches(testCase)
% Khat > K, with one duplicate and one all-zero profile among the estimates.
% The first three estimates are exact, so they win their motifs outright and
% the last two have nothing left to claim: NaN EMD and id 0.
fx = disjoint_fixture();
K = size(fx.W, 2);

[W_dup, H_dup] = perturb_estimates(fx.W(:, 1, :), fx.H(1, :));
W_hat = cat(2, fx.W, W_dup, zeros(size(fx.W(:, 1, :))));
H_hat = [fx.H; H_dup; fx.H(2, :)];

[emdW, emdH, ids] = helper.similarity_WH_EMD(fx.W, fx.H, W_hat, H_hat);

verifyEqual(testCase, ids, [1, 2, 3, 0, 0]);
verifyEqual(testCase, emdW(1:K), zeros(1, K));
verifyEqual(testCase, emdH(1:K), zeros(1, K));
verifyTrue(testCase, all(isnan(emdW(K+1:end))));
verifyTrue(testCase, all(isnan(emdH(K+1:end))));
end

function testTwoDimensionalTruthAgainstThreeDimensionalEstimates(testCase)
% Khat has to be read from size(W_hat), not size(W).  The second estimate
% has no counterpart in the single-motif ground truth.
fx = disjoint_fixture();

W_single = reshape(fx.W(:, 1, :), size(fx.W, 1), size(fx.W, 3));
H_single = fx.H(1, :);

[emdW, emdH, ids] = helper.similarity_WH_EMD( ...
    W_single, H_single, fx.W(:, [1, 2], :), fx.H([1, 2], :));

verifyEqual(testCase, ids, [1, 0]);
verifyEqual(testCase, emdW(1), 0);
verifyEqual(testCase, emdH(1), 0);
verifyTrue(testCase, isnan(emdW(2)));
verifyTrue(testCase, isnan(emdH(2)));
end

% ========================================================================
% Degenerate ground truth
% ========================================================================

function testZeroGroundTruthFactorCannotBeMatched(testCase)
% Ground-truth motif 3 is empty, so it is not available to match against.
% The third estimate is a degraded copy of motif 1, which the exact copy has
% already claimed, leaving it unmatched rather than paired with the empty
% motif or with motif 1 a second time.
fx = disjoint_fixture();

W = fx.W;
H = fx.H;
W(:, 3, :) = 0;
H(3, :) = 0;

[W_dup, H_dup] = perturb_estimates(fx.W(:, 1, :), fx.H(1, :));
W_hat = cat(2, fx.W(:, 1:2, :), W_dup);
H_hat = [fx.H(1:2, :); H_dup];

[emdW, emdH, ids] = helper.similarity_WH_EMD(W, H, W_hat, H_hat);

verifyEqual(testCase, ids, [1, 2, 0]);
verifyEqual(testCase, emdW(1:2), zeros(1, 2));
verifyEqual(testCase, emdH(1:2), zeros(1, 2));
verifyTrue(testCase, isnan(emdW(3)));
verifyTrue(testCase, isnan(emdH(3)));
end

% ========================================================================
% Overlapping motifs: cross-pair costs are finite
% ========================================================================

function testGreedyExclusionWhenTwoEstimatesShareOneTruth(testCase)
% Khat <= K, with both estimates resembling motif 1: one exact copy and one
% degraded copy.  The exact copy must take motif 1, so the degraded copy has
% to fall back to its cheapest unclaimed motif.  Independent per-estimate
% minima would hand motif 1 out twice.
fx = overlapping_fixture();

[W_dup, H_dup, dupTransport] = ...
    perturb_estimates(fx.W(:, 1, :), fx.H(1, :));
W_hat = cat(2, fx.W(:, 1, :), W_dup);
H_hat = [fx.H(1, :); H_dup];

[emdW, emdH, ids] = helper.similarity_WH_EMD(fx.W, fx.H, W_hat, H_hat);

% baseCrossCost(3,1) < baseCrossCost(2,1), so motif 3 is the fallback.
verifyEqual(testCase, ids, [1, 3]);
verifyEqual(testCase, emdW(1), 0);
verifyEqual(testCase, emdH(1), 0);

% Aligning motif 3 to motif 1 costs baseCrossCost(3,1)/uniformMass.  The
% residual inside the degraded copy can move that by at most its own
% transport cost, by the triangle inequality for the transport distance.
fallback = fx.baseCrossCost(3, 1) / fx.uniformMass;
band = dupTransport + solver_rel_tol() * fallback;
verifyGreaterThanOrEqual(testCase, emdW(2), fallback - band);
verifyLessThanOrEqual(testCase, emdW(2), fallback + band);
verifyTrue(testCase, isfinite(emdH(2)));
end

function testOverlappingMotifsMatchOnTransportCost(testCase)
% Every motif in this fixture recruits the same neurons with the same mass,
% so no pairing is forced to discard mass and every cross-pair EMD is a
% finite transport cost.  The matcher has to compare numbers here instead of
% relying on failed pairs being rejected outright.
fx = overlapping_fixture();
permutation = [3, 1, 2];

[W_hat, H_hat, expectedW, expectedH] = ...
    perturb_estimates(fx.W(:, permutation, :), fx.H(permutation, :));

[emdW, emdH, ids] = helper.similarity_WH_EMD(fx.W, fx.H, W_hat, H_hat);

verifyEqual(testCase, ids, permutation);
verifyEqual(testCase, numel(unique(ids)), numel(ids));
verifyEqual(testCase, emdW, expectedW, 'RelTol', solver_rel_tol());
verifyEqual(testCase, emdH, expectedH, 'RelTol', solver_rel_tol());

% The matched costs have to sit below the cheapest wrong pairing, which is
% what makes the correspondence recoverable in the first place.
cheapestWrongPair = min(nonzeros(fx.baseCrossCost)) / fx.uniformMass;
verifyLessThan(testCase, max(emdW), cheapestWrongPair);
end

% ========================================================================
% Tolerances
% ========================================================================

function tol = solver_rel_tol()
% Relative slack for the TFOCS solve that compute_EMD runs.  This is the only
% tolerance in this file that is not derived from a fixture.
%
% similarity_WH_EMD reports ||M||_1 rather than the full objective, so the
% amplified residual lambdaR*||R||_1 no longer sits on top of the asserted
% values.  What is left is the accuracy of the transport term itself, which at
% the default opts.tol = 1e-6 is good: a 3-bin displacement of a unit-mass
% profile comes back as 3.0003, and the worst relative error measured over the
% loading rows of both fixtures is 0.17%.
tol = 5e-2;
end

function f = residual_fraction()
% Size of the residual between a matched pair, as a fraction of the peak
% amplitude of the estimate.
f = 0.1;
end

function b = residual_distance()
% Number of bins the residual mass is displaced by.
b = 2;
end

% ========================================================================
% Fixtures
% ========================================================================

function fx = disjoint_fixture()
% Each motif recruits a single, private neuron.  Transport moves mass only
% along time within a neuron, so a wrong pairing has no shared support to
% transport along and has to discard both motifs' mass, which with
% lambdaR = 1e3 makes the cross-pair EMD huge -- larger than the cost of
% explaining the estimate from an empty profile, so the rejection rule
% discards it.  Cross pairs are failed matches, which makes the expected ids
% unambiguous.
%
% The profiles sit at lags 3:4 of 7 so a shift of +/-2 bins and a residual
% displaced by 2 bins both stay inside the window.  The loading events sit
% well inside 1:T for the same reason.
N = 3;
K = 3;
L = 7;
T = 40;

W = zeros(N, K, L);
W(1, 1, 3:4) = [1.0, 0.5];
W(2, 2, 3:4) = [0.8, 1.0];
W(3, 3, 3:4) = [1.0, 0.6];

H = zeros(K, T);
H(1, [6, 22]) = [1.0, 0.7];
H(2, [10, 26]) = [0.9, 1.0];
H(3, [14, 32]) = [1.0, 0.8];

fx = pack_fixture(W, H);
end

function fx = overlapping_fixture()
% All three motifs recruit the same four neurons with one unit spike each,
% and all three loading rows sum to the same value.  After the internal
% normalisation every motif therefore carries the same mass on every neuron,
% so no pairing is forced into the unmatched-mass term: every cross-pair EMD
% is a pure transport cost.  The motifs differ only in the order the neurons
% fire, ascending, descending and synchronous.
%
% Neuron n fires at lag offset+n-1 in motif 1, offset+N-n in motif 2 and
% offset+1 in motif 3.  Minimising the total displacement over the lag gives
% the raw pairwise costs below, in units of amplitude times bins:
%
%   motifs 1 and 2:  sum_n |5-2n-l|  = 8 at its best lag
%   motifs 1 and 3:  sum_n |2-n-l|   = 4
%   motifs 2 and 3:  sum_n |n-3-l|   = 4
%
% Those are for unit spikes. similarity_WH_EMD scales each profile to unit
% total mass, so every cross-pair cost divides by the common profile mass.
N = 4;
K = 3;
L = 7;
T = 40;
offset = 2;

W = zeros(N, K, L);
for n = 1:N
    W(n, 1, offset + n - 1) = 1;
    W(n, 2, offset + N - n) = 1;
    W(n, 3, offset + 1) = 1;
end

H = zeros(K, T);
H(1, [6, 22]) = 1;
H(2, [11, 27]) = 1;
H(3, [16, 33]) = 1;

fx = pack_fixture(W, H);

assert(all(fx.massH == fx.massH(1)), ...
    'Overlapping fixture needs equal loading row sums.');
assert(all(fx.massW == fx.massW(1)), ...
    'Overlapping fixture needs equal profile masses.');
fx.uniformMass = fx.massW(1);
fx.baseCrossCost = [0, 8, 4; 8, 0, 4; 4, 4, 0];
end

function fx = pack_fixture(W, H)
fx.W = W;
fx.H = H;
% The two per-motif scales that similarity_WH_EMD normalises away: the profile
% total and the loading row sum.
fx.massH = sum(H, 2)';
fx.massW = sum(W, [1, 3]);
end

% ========================================================================
% Fixture manipulation
% ========================================================================

function [W_hat, H_hat, expectedW, expectedH] = perturb_estimates(W_hat, H_hat)
% Displace residual_fraction of the peak amplitude by residual_distance
% bins, in both the profiles and the loadings.  Mass is moved rather than
% added, so every row mass is preserved, no unmatched mass is needed and the
% resulting EMD is exactly the transport cost.
%
% similarity_WH_EMD scales profile k to unit total mass and loading row k to
% unit sum, so the moved mass is delta/massW in W and dh/massH in H.

delta = residual_fraction() * max(W_hat(:));
dh = residual_fraction() * max(H_hat(:));
distance = residual_distance();

N = size(W_hat, 1);
Khat = size(W_hat, 2);
L = size(W_hat, 3);
T = size(H_hat, 2);
massH = sum(H_hat, 2)';
massW = sum(W_hat, [1, 3]);

expectedW = zeros(1, Khat);
expectedH = zeros(1, Khat);

for jj = 1:Khat
    profile = reshape(W_hat(:, jj, :), N, L);
    [n, l] = last_nonzero(profile);
    assert(l + distance <= L && profile(n, l + distance) == 0, ...
        'Profile residual would be truncated or collide.');
    assert(profile(n, l) >= delta, ...
        'Profile residual would drive an amplitude negative.');
    W_hat(n, jj, l) = W_hat(n, jj, l) - delta;
    W_hat(n, jj, l + distance) = delta;

    [~, t] = last_nonzero(H_hat(jj, :));
    assert(t + distance <= T && H_hat(jj, t + distance) == 0, ...
        'Loading residual would be truncated or collide.');
    assert(H_hat(jj, t) >= dh, ...
        'Loading residual would drive a loading negative.');
    H_hat(jj, t) = H_hat(jj, t) - dh;
    H_hat(jj, t + distance) = dh;

    % Moving mass leaves the profile total untouched, so normalising to unit
    % mass just divides the displaced amount by that total.
    expectedW(jj) = delta * distance / massW(jj);
    expectedH(jj) = dh * distance / massH(jj);
end
end

function [row, col] = last_nonzero(A)
[rows, cols] = find(A ~= 0);
[col, k] = max(cols);
row = rows(k);
end

function [W_scaled, H_scaled] = rescale_factors(W, H, scales)
W_scaled = W;

for k = 1:numel(scales)
    W_scaled(:, k, :) = W_scaled(:, k, :) * scales(k);
end

H_scaled = H ./ scales(:);
end
