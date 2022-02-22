function loss = total_loss(W, H, X, lambda)
    [N,K,L] = size(W);
    % Total loss function
    smoothkernel = ones(1,(2*L)-1);  % for factor competition
    WTX = helper.transconv(W, X);
    R = conv2(WTX, smoothkernel, 'same')*H';
    R = lambda*(norm(R(:),1)-norm(diag(R),1));
    loss = norm(helper.reconstruct(W, H)-X, 'fro')^2/2 + R;
end