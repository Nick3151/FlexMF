function y = smooth_cross_ortho_H(W, X, H, mode)
[K, T] = size(H);
[N, ~, L] = size(W);
assert((size(W,2)==K) && (isequal(size(X),[N,T])), 'Dimensions do not match!')
switch mode
    case 0
        y = {[K,T], [K,K]};
    case 1
        WTX = helper.transconv(W, X);
        smoothkernel = ones(1,(2*L)-1);
        WTXS = conv2(WTX, smoothkernel, 'same');
        y = WTXS*H';
    case 2
end