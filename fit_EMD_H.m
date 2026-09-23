function y = fit_EMD_H(W, T, H_, mode)
% Unbalanced-EMD residual operator on H_
% H_ = [H M]'
% f(H_) = div(M)-conv(W,H), so that the residual is R = X + f(H_)
% f*(Y) = [-transconv(W,Y);Y*D]

[N,K,~] = size(W);
% divergence matrix
D = eye(T) - diag(ones(T-1,1),-1);
D(T,T) = 0;

switch mode
    case 0
        y = {[K+N,T], [N,T]};
    case 1
        H = H_(1:K,:);
        M = H_(K+(1:N),:);
        y = M*D'-helper.reconstruct(W,H);
    case 2
        y = [-helper.transconv(W,H_);H_*D];
end
