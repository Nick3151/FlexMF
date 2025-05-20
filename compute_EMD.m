function [d, M, R, out] = compute_EMD(X1, X2, opts, varargin)
% Unbalanced EMD between two matrices, along the second(T) dimension

p  = inputParser;
addOptional(p, 'mu', 1e2);
parse(p, varargin{:})

[N,T] = size(X1);
assert(isequal(size(X2), [N,T]), 'Dimensions of the two matrices should be the same!')
mu = p.Results.mu;

A = @(Y, mode)Beckmann_UOT_constraint(N, T, Y, mode);
W = @(Y, mode)Beckmann_UOT_obj(N, T, mu, Y, mode);

b = X2-X1;
[Y, out] = solver_sBPDN_W(A,W,b,0,.001,[],[],opts);
M = Y(1:N,:);
R = Y(N+1:2*N,:);
d = norm(M(:),1);