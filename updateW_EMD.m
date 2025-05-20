function [W, M, R, out] = updateW_EMD(W0, H, X, M0, R0, params)
% Update W with Earth-mover's distance(EMD) and smooth-orthogonal regularization

[N, K, L] = size(W0);
[~, T] = size(X);
opts_default = tfocs_SCD;
opts = opts_default;
opts.continuation = 1;
opts.tol = 1e-4;
opts.stopCrit = 4;
opts.maxIts = 500;
opts.alg = 'N83';
% opts.debug = true;

%% Initialization
W0_flat = reshape(W0, [N,K*L]);
W0_ = [W0_flat, M0, R0];

%% Linear operators
op_f1 = @(W_, mode)f1_EMD_W(X, H, L, W_, mode);
op_f2 = @(W_, mode)f2_EMD_W(M0, K, L, W_, mode);
op_f3 = @(W_, mode)f3_EMD_W(R0, K, L, W_, mode);
op_f4 = @(W_, mode)f4_EMD_W(H, N, L, W_, mode);

%% Optimize with tfocs
lambda = params.lambda;
lambda_R = params.lambda_R;
mu = 1e-3;

if lambda > 0
    [W_, out] = tfocs_SCD(proj_Rplus_W(K*L), {op_f1, 0; op_f2, 0; op_f3, 0; op_f4, X}, ...
        {proj_linf(lambda), proj_linf(1), proj_linf(lambda_R), proj_Rn}, mu, W0_, [], opts);
else
    [W_, out] = tfocs_SCD(proj_Rplus_W(K*L), {op_f2, 0; op_f3, 0; op_f4, X}, ...
        {proj_linf(1), proj_linf(lambda_R), proj_Rn}, mu, W0_, [], opts);
end

Wflat = W_(:,1:K*L);
W = reshape(Wflat, [N,K,L]);
M = W_(:,K*L+(1:T));
R = W_(:,K*L+T+(1:T));
dW = norm(W(:)-W0(:));
dM = norm(M(:)-M0(:));
dR = norm(R(:)-R0(:));

%% Print intermediate results
if params.verbal
    smoothkernel = ones(1,(2*L)-1);
    WTX = helper.transconv(W, X);
    WTXS = conv2(WTX, smoothkernel, 'same');
    WTXSHT = WTXS*H';
    Q = ones(K);
    Q(1:K+1:end) = 0;
    Xhat = helper.reconstruct(W, H);
    fprintf('reg=%f\n',sum(Q(:).*WTXSHT(:)));
    fprintf('recon=%f\n',sum((X(:)-Xhat(:)).^2)/2);
    fprintf('L1_M=%f\n',sum(abs(M(:))));
    fprintf('L1_R=%f\n',sum(abs(R(:))));
    fprintf('dW=%f\n', dW);
    fprintf('dM=%f\n', dM);
    fprintf('dR=%f\n', dR);
end