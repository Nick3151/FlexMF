function y = tensor_conv_W(H, N, L, W, mode)
% convolution operator on W
[K, T] = size(H);

switch mode
    case 0
        y = {[N,K,L], [N,T]};
    case 1
        y = helper.reconstruct(W, H);
    case 2
        y = zeros([N,K,L]);
        H = [zeros(K,L),H,zeros(K,L)];
        W = [zeros(N,L),W,zeros(N,L)];
        for l = 1 : L
            y(:,:,l) = W * circshift(H, [0,l-1])';
        end
end