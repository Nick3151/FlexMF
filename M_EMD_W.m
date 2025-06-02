function y = M_EMD_W(M0, K, L, W_, mode)
% Motion field M of H_
% W_ = [W_flat M R]
% f(W_) = M
% f*(Y) = [zeros(N,K*L), Y, zeros(N,T)]

[N,T] = size(M0);

switch mode
    case 0
        y = {[N,(K*L+2*T)], [N,T]};
    case 1
        M = W_(:,K*L+(1:T));
        y = M;
    case 2
        y = [zeros(N,K*L), W_, zeros(N,T)];
end