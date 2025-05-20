function y = f4_EMD_W(H, N, L, W_, mode)
% linear constraint operator on W_
% W_ = [W_flat M R]
% f4(W_) = div(M)-R-conv(W,H)
% f4*(Y) = [-flatten(Y*H(-->l)'), Y*D, -Y]

[K,T] = size(H);
% divergence matrix
D = eye(T) - diag(ones(T-1,1),-1);

switch mode
    case 0
        y = {[N,(K*L+2*T)], [N,T]};
    case 1
        W_flat = W_(:,1:K*L);
        W = reshape(W_flat, [N,K,L]);   
        M = W_(:,K*L+(1:T));
        R = W_(:,K*L+T+(1:T));
        y = M*D'-R-helper.reconstruct(W,H);
    case 2
        y_tmp = zeros([N,K,L]);
        H_pad = [zeros(K,L),H,zeros(K,L)];
        W_pad = [zeros(N,L),W_,zeros(N,L)];
        for l = 1 : L
            y_tmp(:,:,l) = W_pad * circshift(H_pad, [0,l-1])';
        end
        y = [-reshape(y_tmp, [N,K*L]), W_*D, -W_];
end