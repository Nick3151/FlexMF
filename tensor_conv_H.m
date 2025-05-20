function y = tensor_conv_H(W, T, H, mode)
% convolution operator on H
[N, K, L] = size(W);

switch mode
    case 0
        y = {[K,T], [N,T]};
    case 1
        y = helper.reconstruct(W, H);
    case 2
        y = helper.transconv(W, H);
end