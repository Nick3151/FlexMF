function plot_range = estimate_plot_range(H, n_occurrences)
%ESTIMATE_PLOT_RANGE Select [1, end] covering occurrences of every sequence.
if nargin < 2 || isempty(n_occurrences)
    n_occurrences = 5;
end

K = size(H, 1);
T = size(H, 2);
active_thresh = 1e-3 * max(abs(H(:)));
last_onset = T * ones(1, K);
for k = 1:K
    onsets = find(diff([0, abs(H(k, :)) > active_thresh]) == 1);
    if ~isempty(onsets)
        last_onset(k) = onsets(min(n_occurrences, numel(onsets)));
    end
end
plot_range = [1, min(T, max(last_onset))];
end
