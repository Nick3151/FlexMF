function [d, M, R, out] = compute_EMD(X1, X2, opts, varargin)
% Unbalanced EMD between two matrices, or two 1d sequences, 
% along the time dimension. d is the full UOT objective
% ||M||_1 + lambdaR*||R||_1.
% R = X1 + M*D' - X2 is eliminated from the solve and penalized as
% lambdaR*||R||_1, so the transport constraint holds exactly.

p  = inputParser;
addOptional(p, 'lambdaR', 1e1);
addOptional(p, 'continuationOptions', [])
parse(p, varargin{:})

[N,T] = size(X1);
assert(isequal(size(X2), [N,T]), 'Dimensions of the two matrices should be the same!')
lambdaR = p.Results.lambdaR;
continuationOptions = p.Results.continuationOptions;

if isequal(X1,X2)
    d = 0;
    M = zeros(N,T);
    R = zeros(N,T);
    out = [];
else
    if isempty(continuationOptions)
        continuationOptions = continuation();
    end
    
    op_fit = @(Y, mode)fit_EMD(N, T, Y, mode);
    
    [M, out] = tfocs_SCD(prox_l1, {op_fit, X1-X2}, proj_linf(lambdaR), .1, [], [], opts, continuationOptions);
    R = helper.correct_warp(X1, M) - X2;
    d = norm(M(:),1) + lambdaR*norm(R(:),1);
end