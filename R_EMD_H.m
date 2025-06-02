function y = R_EMD_H(R0, K, H_, mode)
% Residual Error R of H_
% H_ = [H M R]'
% f(H_) = R
% f*(Y) = [zeros(K+N,T);Y]

[N,T] = size(R0);

switch mode
    case 0
        y = {[K+2*N,T], [N,T]};
    case 1
        R = H_(K+N+(1:N),:);
        y = R;
    case 2
        y = [zeros(K+N,T);H_];
end