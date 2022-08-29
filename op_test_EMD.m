clear all
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
N = 20;
T = 500;
X = rand(N, T);
mu = 100;

op = @(X, mode)Beckmann_UOT_obj(N, T, mu, X, mode);
linop_test(op, 'R2R');

op = @(X, mode)Beckmann_UOT_constraint(N, T, X, mode);
linop_test(op, 'R2R');