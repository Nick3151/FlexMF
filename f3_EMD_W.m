function y = f3_EMD_W(R, K, L, W_, mode)
% W_ = [W_flat M R]
% f3(W_) = R
% f3*(Y) = [zeros(N,(K*L+T)), Y]

[N,T] = size(R);

switch mode
    case 0
        y = {[N,(K*L+2*T)], [N,T]};
    case 1
        R = W_(:,K*L+T+(1:T));
        y = R;
    case 2
        y = [zeros(N,K*L+T), W_];
end