function y = fit_EMD_W(H, N, L, W_, mode)
% Unbalanced-EMD residual operator on W_
% W_ = [W_flat M]
% f(W_) = div(M)-conv(W,H), so that the residual is R = X + f(W_)
% f*(Y) = [-flatten(Y*H(-->l)'), Y*D]

[K,T] = size(H);
% divergence matrix
D = eye(T) - diag(ones(T-1,1),-1);
D(T,T) = 0;

switch mode
    case 0
        y = {[N,(K*L+T)], [N,T]};
    case 1
        W_flat = W_(:,1:K*L);
        W = reshape(W_flat, [N,K,L]);
        M = W_(:,K*L+(1:T));
        y = M*D'-helper.reconstruct(W,H);
    case 2
        y_tmp = zeros([N,K,L]);
        H_pad = [zeros(K,L),H,zeros(K,L)];
        W_pad = [zeros(N,L),W_,zeros(N,L)];
        for l = 1 : L
            y_tmp(:,:,l) = W_pad * circshift(H_pad, [0,l-1])';
        end
        y = [-reshape(y_tmp, [N,K*L]), W_*D];
end
