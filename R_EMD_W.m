function y = R_EMD_W(R0, K, L, W_, mode)
% Residual Error R of W_
% W_ = [W_flat M R]
% f(W_) = R
% f*(Y) = [zeros(N,(K*L+T)), Y]

[N,T] = size(R0);

switch mode
    case 0
        y = {[N,(K*L+2*T)], [N,T]};
    case 1
        R = W_(:,K*L+T+(1:T));
        y = R;
    case 2
        y = [zeros(N,K*L+T), W_];
end