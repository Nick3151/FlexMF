function score = compute_error_balance_score_3d(X,K,L,lambda,lambdaL1W,alpha)
% a = lambda*|||WXS|HT||_1
% b = lambdaL1W*||W||_1
% c = ||X-X_hat||_2
% score = (a^3+b^3+c^3)/3 / ((a+b+c)/3)^3

[W_hat,H_hat,~,~,loadings,power] = FlexMF(X,'K',K, 'L', L, 'maxiter', 50,...
    'lambda', lambda, 'alpha', alpha, 'lambdaL1W', lambdaL1W, 'lambdaL1H', 0, 'neg_prop', 0.2, 'showPlot', 0);
[recon_error, reg_cross, reg_W, reg_H] = helper.get_FlexMF_cost(X,W_hat,H_hat);

epsilon = 1e-6;
recon_error = recon_error + epsilon;
reg_cross = reg_cross*lambda + epsilon;
reg_W = reg_W*lambdaL1W + epsilon;

score = (recon_error^3 + reg_cross^3 + reg_W^3)/3 / ((recon_error+reg_cross+reg_W)/3)^3;
end