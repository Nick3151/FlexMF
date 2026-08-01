function Hout = shift_loadings(H, s, Tout)
% SHIFT_LOADINGS  Move temporal loadings along the time axis, zero filled.
%
%   Hout = helper.shift_loadings(H, s)
%   Hout = helper.shift_loadings(H, s, Tout)
%
% H is K x T and every row is moved s time bins later (positive s means
% later). Mass pushed past either end of the recording is dropped rather than
% wrapped around, so this is a shift in time and not a circular shift.
%
% Tout sets the length of the output and defaults to T.
%
% Moving a motif profile s bins later and its loadings s bins earlier leaves
% the reconstruction W*H unchanged, so undoing a profile lag on the loadings
% means calling this with the negated lag.
%
% See also helper.shift_profiles

assert(isscalar(s) && isnumeric(s) && isfinite(s) && s == round(s), ...
    'helper:shift_loadings:InvalidShift', ...
    'Shift must be a finite integer scalar.');

[K, T] = size(H);

if nargin < 3 || isempty(Tout)
    Tout = T;
end

assert(isscalar(Tout) && Tout >= 1 && Tout == round(Tout), ...
    'helper:shift_loadings:InvalidWidth', ...
    'Tout must be a positive integer scalar.');

Hout = zeros(K, Tout);

for t = 1:T
    target = t + s;
    if target >= 1 && target <= Tout
        Hout(:, target) = H(:, t);
    end
end
