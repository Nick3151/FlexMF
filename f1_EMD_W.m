function y = f1_EMD_W(X, H, L, W_, mode)
% off-diagnal part of the smooth cross orthogonal operator on W_
% W_ = [W_flat M R]
% W_flat: N*(KL)
% f1(W_) = Q.*(transconv(W,X)*S*H')
% f1*(Y) = [flatten(X(<--l)*S*H'*(Y.*Q)'); zeros(N,2T)]

[N, T] = size(X);
[K, ~] = size(H);
assert((size(H,2)==T), 'Dimensions do not match!')

Q = ones(K);
Q(1:K+1:end) = 0;   % off diagonal mask
smoothkernel = ones(1,(2*L)-1);

switch mode
    case 0
        y = {[N,(K*L+2*T)], [K,K]};
    case 1
        W_flat = W_(:,1:K*L);
        W = reshape(W_flat, [N,K,L]);       
        WTX = helper.transconv(W, X);        
        WTXS = conv2(WTX, smoothkernel, 'same');
        y = Q.*(WTXS*H');
    case 2
        y_tmp = zeros([N,K,L]);
        H = [zeros(K,L),H,zeros(K,L)];
        X = [zeros(N,L),X,zeros(N,L)];
        SHT = conv2(H, smoothkernel, 'same')';
        for l = 1 : L
            y_tmp(:,:,l) = circshift(X, [0,-l+1])*SHT*(Q.*W_)';
        end
        y = [reshape(y_tmp, [N,K*L]), zeros(N,2*T)];
end