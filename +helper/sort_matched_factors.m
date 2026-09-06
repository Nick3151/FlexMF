function [W_sorted, H_sorted] = sort_matched_factors(W, H, ids)
%SORT_MATCHED_FACTORS  Put matched estimates at their ground-truth indices.
%
%   ids(j) = i means estimated factor j matched ground-truth factor i (0 = none).
%   Matched estimates are placed at position i; unmatched ones fill the rest.
Khat = numel(ids);
factor_order = zeros(1, Khat);
used = false(1, Khat);

matched = find(ids > 0);
for jj = matched
    ii = ids(jj);
    if ii >= 1 && ii <= Khat && factor_order(ii) == 0
        factor_order(ii) = jj;
        used(jj) = true;
    end
end

unmatched = find(~used);
empty = find(factor_order == 0);
n = min(numel(empty), numel(unmatched));
factor_order(empty(1:n)) = unmatched(1:n);

W_sorted = W(:, factor_order, :);
H_sorted = H(factor_order, :);
end
