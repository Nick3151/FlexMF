function y = Beckmann_UOT_constraint(N, T, X, mode)
% The linear operator for the constraint of unbalanced EMD
% X = [M; R] size 2N*T
% A(X) = MD'-R
% A*(y) = [YD;-Y]
% divergence matrix
D = eye(T) - diag(ones(T-1,1),-1);
switch mode
    case 0
        y = {[2*N,T], [N,T]};
    case 1
        M = X(1:N,:);
        R = X(N+1:end,:);
        y = M*D'-R;
    case 2
        y = zeros(2*N, T);
        y(1:N,:) = X*D;
        y(N+1:end,:) = -X;
end