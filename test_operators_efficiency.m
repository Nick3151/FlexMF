clear all 
close all
clc

K = 3;
N = 30;
L = 50;
T = 100;


W = randn(N,K,L);
H = randn(K,T);
tic
Xhat = helper.reconstruct(W,H);
toc

tic
X = 0;
for k=1:K
    Wk = squeeze(W(:,k,:));
    Hk = H(k,:);
    X = X+convn(Wk, Hk, 'full');
end
toc