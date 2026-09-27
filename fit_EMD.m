function y = fit_EMD(N, T, M_, mode)
% Divergence operator on M alone (no W,H factorization term)
% f(M_) = M_*D', so that the residual is R = (X1-X2) + f(M_)
% f*(Y) = Y*D

D = eye(T) - diag(ones(T-1,1),-1);
D(T,T) = 0;

switch mode
    case 0
        y = {[N,T], [N,T]};
    case 1
        y = M_*D';
    case 2
        y = M_*D;
end
