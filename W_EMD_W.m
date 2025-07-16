function y = W_EMD_W(M0, K, L, W_, mode)
% W of H_
% W_ = [W_flat M R]
% f(W_) = W
% f*(Y) = [Y, zeros(N,2*T)]

[N,T] = size(M0);

switch mode
    case 0
        y = {[N,(K*L+2*T)], [N,K,L]};
    case 1
        W_flat = W_(:,1:K*L);
        y = reshape(W_flat, [N,K,L]);
    case 2
        y = [reshape(W_, [N,K*L]), zeros(N,2*T)];
end