function y = smooth_cross_ortho_H(W, X, H, mode)
% smooth cross orthogonal operator on H
[N, K, L] = size(W);
[~, T] = size(X);
assert((size(X,1)==N), 'Dimensions do not match!')
switch mode
    case 0
        y = {[K,T], [K,K]};
    case 1
        WTX = helper.transconv(W, X);
        smoothkernel = ones(1,(2*L)-1);
        WTXS = conv2(WTX, smoothkernel, 'same');
        y = WTXS*H';
    case 2
        WTX = helper.transconv(W, X);
        smoothkernel = ones(1,(2*L)-1);
        WTXS = conv2(WTX, smoothkernel, 'same');
        y = H'* WTXS;
end