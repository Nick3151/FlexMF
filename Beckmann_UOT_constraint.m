function y = Beckmann_UOT_constraint(K, T, X, mode)
% The linear operator for the constraint of unbalanced EMD
% X = [M^T; R]
% A(X) = [DM';DM';...;DM']-R
% A*(y) = [sum_i(yi*D);y1;y2;...;yk]
% divergence matrix
D = eye(T) - diag(ones(T-1,1),-1);
switch mode
    case 0
        y = {[K,T], [K-1,T]};
    case 1
        M = X(1,:)';
        R = X(2:end,:);
        y = repmat((D*M)', K-1, 1)-R;
    case 2
        y = zeros(K, T);
        y(1, :) = sum(X,1)*D;
        y(2:end, :) = -X;
end