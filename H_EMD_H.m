function y = H_EMD_H(M0, K, H_, mode)
% H of H_
% H_ = [H M]'
% f(H_) = H
% f*(Y) = [Y; zeros(N,T)]

[N,T] = size(M0);

switch mode
    case 0
        y = {[K+N,T], [K,T]};
    case 1
        H = H_(1:K,:);
        y = H;
    case 2
        y = [H_; zeros(N,T)];
end