function y = W_EMD_W(M0, K, L, W_, mode)
% W of W_
% W_ = [W_flat M]
% f(W_) = W
% f*(Y) = [Y, zeros(N,T)]

[N,T] = size(M0);

switch mode
    case 0
        y = {[N,(K*L+T)], [N,K,L]};
    case 1
        W_flat = W_(:,1:K*L);
        y = reshape(W_flat, [N,K,L]);
    case 2
        y = [reshape(W_, [N,K*L]), zeros(N,T)];
end