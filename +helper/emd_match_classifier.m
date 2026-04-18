function [is_match, ratio] = emd_match_classifier(emd_curve, center_fraction, threshold)
% EMD_MATCH_CLASSIFIER  Distinguish successful from failed sequence matches
% based on the shape of the EMD-vs-time-shift curve.
%
% A successful match produces a clear local minimum near the center of the
% curve (zero-shift region). A failed match (different neurons activated)
% produces a curve with no such dip — the center value is comparable to or
% higher than the flanks.
%
% The classifier computes the "central dip ratio":
%   ratio = min(EMD in central window) / mean(EMD in flanks)
% A ratio well below 1 indicates a deep central dip → successful match.
%
% USAGE:
%   [is_match, ratio] = emd_match_classifier(emd_curve)
%   [is_match, ratio] = emd_match_classifier(emd_curve, center_fraction, threshold)
%
% INPUTS:
%   emd_curve        - 1-D array of EMD values computed at each time shift.
%                      The array is assumed to be ordered from most negative
%                      to most positive shift, with zero shift at the center.
%   center_fraction  - Fraction of the curve length treated as the central
%                      window (default: 0.4, i.e. middle 40%).
%   threshold        - Ratio below which a central dip is declared a match
%                      (default: 0.85). Tune this on labelled examples from
%                      your own data.
%
% OUTPUTS:
%   is_match         - Logical scalar. true = successful match.
%   ratio            - The central dip ratio (lower → stronger match signal).
%
% EXAMPLE:
%   % Simulate a successful match curve
%   shifts = -20:20;
%   emd_success = 480 + 0.45*shifts.^2 - 120*exp(-shifts.^2/30) + randn(size(shifts))*5;
%
%   % Simulate a failed match curve (no central dip)
%   emd_fail    = 600 + 40*exp(-shifts.^2/800) + 0.15*shifts.^2 + randn(size(shifts))*5;
%
%   [m1, r1] = emd_match_classifier(emd_success);
%   [m2, r2] = emd_match_classifier(emd_fail);
%
%   fprintf('Success curve → match=%d, ratio=%.3f\n', m1, r1);
%   fprintf('Fail    curve → match=%d, ratio=%.3f\n', m2, r2);
 
% ── Defaults ─────────────────────────────────────────────────────────────
if nargin < 2 || isempty(center_fraction)
    center_fraction = 0.4;
end
if nargin < 3 || isempty(threshold)
    threshold = 0.85;
end
 
% ── Input validation ──────────────────────────────────────────────────────
emd_curve = emd_curve(:);          % ensure column vector
n = numel(emd_curve);
assert(n >= 5, 'emd_curve must have at least 5 elements.');
assert(center_fraction > 0 && center_fraction < 1, ...
    'center_fraction must be in (0, 1).');
assert(threshold > 0 && threshold <= 1, ...
    'threshold must be in (0, 1].');
 
% ── Central window indices ────────────────────────────────────────────────
mid = ceil(n / 2);                          % 1-based centre index
half_w = max(1, round(n * center_fraction / 2));
 
c_start = max(1,   mid - half_w);
c_end   = min(n,   mid + half_w);
 
% ── Compute ratio ────────────────────────────────────────────────────────
center_vals = emd_curve(c_start : c_end);
flank_vals  = emd_curve([1 : c_start-1, c_end+1 : n]);
 
if isempty(flank_vals)
    warning('emd_match_classifier:noFlanks', ...
        'center_fraction is so large there are no flank samples. Returning ratio=1.');
    ratio    = 1;
    is_match = false;
    return
end
 
center_min  = min(center_vals);
flanks_mean = mean(flank_vals);
 
if flanks_mean == 0
    warning('emd_match_classifier:zeroFlanks', ...
        'Flank mean is zero; cannot compute ratio. Returning ratio=0.');
    ratio    = 0;
    is_match = true;
    return
end
 
ratio    = center_min / flanks_mean;
is_match = ratio < threshold;
 
end