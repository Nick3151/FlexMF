function [f, g] = off_diag_frob_norm_sqr(lambda, Q, X, t)
% 0.5*lambda*norm(Q.*X, 'fro')^2
% Q_ij=1 when i /neq j
% Q_ij=0 when i=j
% gradient = lambda*sum(Q.*X, 'all')

assert(isequal(size(Q), size(X)), 'Dimension does not match!');
if nargin==4
    error('This function does not support proximity minimization.')
elseif nargin==3
    g = lambda*Q.*X;
    f = 0.5*lambda*sum(Q(:).*(X(:).^2));
end

    