function y = obj_EMD_H(M0, K, mu, H_, mode)
% objective function of H_
% H_ = [H M R]'
% f(H_) = [M mu*R]'
% Y = [Y1; Y2]
% f*(Y) = [zeros(K,T);Y1;mu*Y2]

[N,T] = size(M0);

switch mode
    case 0
        y = {[K+2*N,T], [2*N,T]};
    case 1
        M = H_(K+(1:N),:);
        R = H_(K+N+(1:N),:);
        y = [M;mu*R];
    case 2
        Y1 = H_(1:N,:);
        Y2 = H_(N+(1:N),:);
        y = [zeros(K,T);Y1;mu*Y2];
end