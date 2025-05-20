function y = f2_EMD_W(M, K, L, W_, mode)
% W_ = [W_flat M R]
% f2(W_) = M
% f2*(Y) = [zeros(N,K*L), Y, zeros(N,T)]

[N,T] = size(M);

switch mode
    case 0
        y = {[N,(K*L+2*T)], [N,T]};
    case 1
        M = W_(:,K*L+(1:T));
        y = M;
    case 2
        y = [zeros(N,K*L), W_, zeros(N,T)];
end