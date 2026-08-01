function Wout = shift_profiles(W, s, Lout)
% SHIFT_PROFILES  Move motif profiles along the lag axis, zero filled.
%
%   Wout = helper.shift_profiles(W, s)
%   Wout = helper.shift_profiles(W, s, Lout)
%
% W is N x K x L, or N x L for a single motif, and every profile is moved s
% lag bins later (positive s means later). Mass pushed past either end of the
% window is dropped rather than wrapped around, so this is a shift of the lag
% window and not a circular shift.
%
% Lout sets the width of the output window and defaults to L. Passing the
% width of a set of estimated profiles lets a shifted ground-truth profile be
% compared against them directly, which is how the motif matching searches
% over lags.
%
% A 2-D input returns a 2-D output.
%
% See also helper.shift_loadings

assert(isscalar(s) && isnumeric(s) && isfinite(s) && s == round(s), ...
    'helper:shift_profiles:InvalidShift', ...
    'Shift must be a finite integer scalar.');

szW = size(W);
is2D = numel(szW) == 2;

if is2D
    W = reshape(W, [szW(1), 1, szW(2)]);
end

[N, K, L] = size(W);

if nargin < 3 || isempty(Lout)
    Lout = L;
end

assert(isscalar(Lout) && Lout >= 1 && Lout == round(Lout), ...
    'helper:shift_profiles:InvalidWidth', ...
    'Lout must be a positive integer scalar.');

Wout = zeros(N, K, Lout);

for l = 1:L
    target = l + s;
    if target >= 1 && target <= Lout
        Wout(:, :, target) = W(:, :, l);
    end
end

if is2D
    Wout = reshape(Wout, [N, Lout]);
end
