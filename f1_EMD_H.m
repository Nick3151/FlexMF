function y = f1_EMD_H(A, N, H_, mode)
% off-diagnal part of the smooth cross orthogonal operator on H_
% A = WTXS 
% H_ = [H M R]'
% f1*(Y) = [(Y.*Q)'*A; zeros(2*N,T)]

[K, T] = size(A);
Q = ones(K);
Q(1:K+1:end) = 0;   % off diagonal mask

switch mode
    case 0
        y = {[K+2*N,T], [K,K]};
    case 1
        H = H_(1:K,:);
        y = Q.*(A*H');
    case 2
        y = [(Q.*H_)'*A; zeros(2*N,T)];
end