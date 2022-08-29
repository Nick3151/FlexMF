function y = Beckmann_UOT_obj(K, T, mu, X, mode)
% The linear operator for the objective function of unbalanced EMD
% X = [M^T; R]
% A(X) = [M^T; mu*R]
% A*(y) = [y1; mu*y2]

switch mode
    case 0
        y = {[K,T], [K,T]};
    case 1
        y = zeros(K, T);
        y(1, :) = X(1, :);
        y(2:end, :) = mu*X(2:end, :);
    case 2
        y = zeros(K, T);
        y(1, :) = X(1, :);
        y(2:end, :) = mu*X(2:end, :);
end