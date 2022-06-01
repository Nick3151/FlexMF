function y = smooth_cross_ortho_H(A, K, H, mode)
% smooth cross orthogonal operator on H
% A = WTXS or [WTXS; I]
[M, T] = size(A);

switch mode
    case 0
        y = {[K,T], [M,K]};
    case 1
        y = A*H';
    case 2
        y = H'*A;
end