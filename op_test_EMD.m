clear all
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
N = 20;
T = 500;
X = rand(N, T);
mu = 100;

% op = @(X, mode)Beckmann_UOT_obj(N, T, mu, X, mode);
% linop_test(op, 'R2R');
% 
% op = @(X, mode)Beckmann_UOT_constraint(N, T, X, mode);
% linop_test(op, 'R2R');

K = 10;
A = rand(K, T);
op = @(H_, mode)f1_EMD_H(A, N, H_, mode);
linop_test(op, 'R2R');

M = rand(N, T);
op = @(H_, mode)f2_EMD_H(M, K, H_, mode);
linop_test(op, 'R2R');

R = rand(N, T);
op = @(H_, mode)f3_EMD_H(R, K, H_, mode);
linop_test(op, 'R2R');

L = 50;
W = rand(N,K,L);
op = @(H_, mode)f4_EMD_H(W, T, H_, mode);
linop_test(op, 'R2R');

H = rand(K,T);
X_pad = [zeros(N,L),X,zeros(N,L)];
H_pad = [zeros(K,L),H,zeros(K,L)];
op = @(W_, mode)f1_EMD_W(X_pad, H_pad, L, W_, mode);
linop_test(op, 'R2R');

op = @(W_, mode)f2_EMD_W(M, K, L, W_, mode);
linop_test(op, 'R2R');

op = @(W_, mode)f3_EMD_W(R, K, L, W_, mode);
linop_test(op, 'R2R');

op = @(W_, mode)f4_EMD_W(H, N, L, W_, mode);
linop_test(op, 'R2R');