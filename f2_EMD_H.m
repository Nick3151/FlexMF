function y = f2_EMD_H(M0, K, H_, mode)
% H_ = [H M R]'
% f2(H_) = M
% f2*(Y) = [zeros(K,T);Y;zeros(N,T)]

[N,T] = size(M0);

switch mode
    case 0
        y = {[K+2*N,T], [N,T]};
    case 1
        M = H_(K+(1:N),:);
        y = M;
    case 2
        y = [zeros(K,T);H_;zeros(N,T)];
end