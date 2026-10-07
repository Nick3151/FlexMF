function s = hoyer_sparsity(H)
% Mean Hoyer sparsity over all K rows of H (K x T).
% Per row: (sqrt(T) - ||h||_1 / ||h||_2) / (sqrt(T) - 1), in [0, 1];
% 1 = a single nonzero entry, 0 = constant row. All-zero rows count as 1,
% so the mean is always over K rows regardless of how many factors survive.

T = size(H, 2);
l1 = sum(abs(H), 2);
l2 = sqrt(sum(H.^2, 2));
s_rows = (sqrt(T) - l1 ./ l2) / (sqrt(T) - 1);
s_rows(l2 == 0) = 1;
s = mean(s_rows);

end
