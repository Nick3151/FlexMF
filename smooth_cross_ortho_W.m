function y = smooth_cross_ortho_W(X, H, L, W, mode)
% smooth cross orthogonal operator on H
[N, T] = size(X);
[K, ~] = size(H);
assert((size(H,2)==T), 'Dimensions do not match!')
switch mode
    case 0
        y = {[N,K,L], [K,K]};
    case 1
        WTX = helper.transconv(W, X);
        smoothkernel = ones(1,(2*L)-1);
        WTXS = conv2(WTX, smoothkernel, 'same');
        y = WTXS*H';
    case 2
        y = zeros([N,K,L]);
        H = [zeros(K,L),H,zeros(K,L)];
        X = [zeros(N,L),X,zeros(N,L)];
        smoothkernel = ones(1,(2*L)-1);
        SHT = conv2(H, smoothkernel, 'same')';
        for l = 1 : L
            y(:,:,l) = circshift(X, [0,-l+1])'*SHT*W;
        end
end