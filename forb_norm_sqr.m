function [f, G] = forb_norm_sqr(alpha, X, t)
if nargin == 3
    error('This function does not support proximity minimization.');
else
    G = alpha*X;
    f = alpha/2*sum(X.^2, 'all');
end