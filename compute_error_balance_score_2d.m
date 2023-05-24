function score = compute_error_balance_score_2d(X,K,L,lambda,alpha)
% a = lambda*|||WXS|HT||_1
% b = ||X-X_hat||_2
% score = (a^2+b^2)/2 / ((a+b)/2)^2

[W_hat,H_hat,~,errors,loadings,power] = FlexMF(X,'K',K, 'L', L, 'maxiter', 50,...
    'lambda', lambda, 'alpha', alpha, 'lambdaL1W', 0, 'lambdaL1H', 0, 'neg_prop', 0, 'showPlot', 0);
[recon_error, reg_cross, reg_W, reg_H] = helper.get_FlexMF_cost(X,W_hat,H_hat);

epsilon = 1e-6;
scale = 1.5;    % Make recon error a bit larger than regularization error
recon_error = recon_error + epsilon;
reg_cross = reg_cross*lambda*scale + epsilon;

score = (recon_error^2 + reg_cross^2)/2 / ((recon_error+reg_cross)/2)^2;
end