function y = constraint_EMD_H(W, T, H_, mode)
% linear constraint operator on H_
% H_ = [H M R]'
% f(H_) = div(M)-R-conv(W,H)
% f*(Y) = [-transconv(W,Y);Y*D;-Y]

[N,K,L] = size(W);
% divergence matrix
D = eye(T) - diag(ones(T-1,1),-1);

switch mode
    case 0
        y = {[K+2*N,T], [N,T]};
    case 1
        H = H_(1:K,:);
        M = H_(K+(1:N),:);
        R = H_(K+N+(1:N),:);
        y = M*D'-R-helper.reconstruct(W,H);
    case 2
        y = [-helper.transconv(W,H_);H_*D;-H_];
end
