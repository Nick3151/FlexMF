function score = compute_error_balance_score(X,K,L,lambda,lambdaL1H,lambdaL1W,alpha)
% a = lambda*|||WXS|HT||_1
% b = lambdaL1W*||W||_1
% c = ||X-X_hat||_2
% score = (a+b+c)/(3*(abc)^(1/3))

[W_hat,H_hat,cost,loadings,power] = FlexMF(X,'K',K, 'L', L, 'maxiter', 50,...
    'lambda', lambda, 'alpha', alpha, 'lambdaL1W', lambdaL1W, 'lambdaL1H', lambdaL1H, 'showPlot', 0);
[recon_error, reg_cross, reg_W, reg_H] = helper.get_FlexMF_cost(X,W_hat,H_hat);
recon_error = recon_error/2;
reg_cross = reg_cross*lambda;
reg_W = reg_W*lambdaL1W;
reg_H = reg_H*lambdaL1H;
score = (recon_error+reg_cross+reg_W)/(3*(recon_error*reg_cross*reg_W)^(1/3));
end