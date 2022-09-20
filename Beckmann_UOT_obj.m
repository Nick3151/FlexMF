function y = Beckmann_UOT_obj(N, T, mu, X, mode)
% The linear operator for the objective function of unbalanced EMD
% X = [M; R] size 2N*T
% A(X) = [M; mu*R]
% A*(Y) = [Y1; mu*Y2]

switch mode
    case 0
        y = {[2*N,T], [2*N,T]};
    case 1
        y = zeros(2*N, T);
        y(1:N, :) = X(1:N, :);
        y(N+1:end, :) = mu*X(N+1:end, :);
    case 2
        y = zeros(2*N, T);
        y(1:N, :) = X(1:N, :);
        y(N+1:end, :) = mu*X(N+1:end, :);
end