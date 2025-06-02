function y = M_EMD_H(M0, K, H_, mode)
% Motion field M of H_
% H_ = [H M R]'
% f(H_) = M
% f*(Y) = [zeros(K,T);Y;zeros(N,T)]

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