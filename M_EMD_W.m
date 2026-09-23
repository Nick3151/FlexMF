function y = M_EMD_W(M0, K, L, W_, mode)
% Motion field M of W_
% W_ = [W_flat M]
% f(W_) = M
% f*(Y) = [zeros(N,K*L), Y]

[N,T] = size(M0);

switch mode
    case 0
        y = {[N,(K*L+T)], [N,T]};
    case 1
        M = W_(:,K*L+(1:T));
        y = M;
    case 2
        y = [zeros(N,K*L), W_];
end