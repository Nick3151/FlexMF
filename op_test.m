clear all
close all
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
N = 20;
K = 10;
L = 30;
T = 500;
W = rand(N, K, L);
X = rand(N, T);
H = rand(K, T);
X = [zeros(N,L),X,zeros(N,L)];
H = [zeros(K,L),H,zeros(K,L)];
T = size(X, 2);
smoothkernel = ones(1,(2*L)-1);  % for factor competition
WTX = helper.transconv(W, X);
WTXS = conv2(abs(WTX), smoothkernel, 'same');
A = [WTXS; eye(T)];

op = @(H, mode)smooth_cross_ortho_H(A, K, H, mode);
linop_test(op, 'R2R');

op = @(H, mode)tensor_conv_H(W, T, H, mode);
linop_test(op, 'R2R');

op = @(W, mode)tensor_conv_W(H, N, L, W, mode);
linop_test(op, 'R2R');

op = @(W, mode)smooth_cross_ortho_W(X, H, L, W, mode);
linop_test(op, 'R2R');