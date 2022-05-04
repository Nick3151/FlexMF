function [f, Z] = off_diag_l1_norm(lambda, X, t)
if nargin==3
    Z = sign(X).*max(abs(X)-t/(2*lambda), 0); %apply shrinkage operator
    n = size(X,1);
    Z(1:n+1:end) = diag(X); %keep the diagonal entries the same in proximal operator
elseif nargout==2
    error('This function is not differentiable.')
end
f = (sum(abs(X), 'all')-sum(diag(X)))/(2*lambda);
    