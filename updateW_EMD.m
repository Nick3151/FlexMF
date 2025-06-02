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
op_cross_orth_W = @(W_, mode)cross_orth_EMD_W(X, H, L, W_, mode);
op_M = @(W_, mode)M_EMD_W(M0, K, L, W_, mode);
op_R = @(W_, mode)R_EMD_W(R0, K, L, W_, mode);
op_W = @(W_, mode)W_EMD_W(M0, K, L, W_, mode);
op_constraint = @(W_, mode)constraint_EMD_W(H, N, L, W_, mode);

norm_cross_orth2 = linop_normest(op_cross_orth_W).^2;
norm_M2 = linop_normest(op_M).^2;
norm_R2 = linop_normest(op_R).^2;
norm_W2 = linop_normest(op_W).^2;
norm_constraint2 = linop_normest(op_constraint).^2;

proxScale_corss_orth = sqrt(norm_cross_orth2/norm_constraint2);
proxScale_M = sqrt(norm_M2/norm_constraint2);
proxScale_R = sqrt(norm_R2/norm_constraint2);
proxScale_W = sqrt(norm_W2/norm_constraint2);

%% Optimize with tfocs
lambda = params.lambda;
lambda_R = params.lambda_R;
lambdaL1W = params.lambdaL1W;
mu = 1e-3;

affineF = {linop_compose(op_M, 1/proxScale_M), 0; ...
           linop_compose(op_R, 1/proxScale_R), 0; ...
           op_constraint, X};
conjnegF = {proj_linf(proxScale_M), proj_linf(lambda_R*proxScale_R), proj_Rn};

if lambda>0
    affineF(end+1,:) = {linop_compose(op_cross_orth_W, 1/proxScale_corss_orth), 0};
    conjnegF{end+1} = proj_linf(lambda*proxScale_corss_orth);
end

if lambdaL1W>0
    affineF(end+1,:) = {linop_compose(op_W, 1/proxScale_W), 0};
    conjnegF{end+1} = proj_linf(lambdaL1W*proxScale_W);
end

[W_, out] = tfocs_SCD(proj_Rplus_W(K*L), affineF, conjnegF, mu, W0_, [], opts);

Wflat = W_(:,1:K*L);
W = reshape(Wflat, [N,K,L]);
M = W_(:,K*L+(1:T));
R = W_(:,K*L+T+(1:T));
dW = norm(W(:)-W0(:));
dM = norm(M(:)-M0(:));
dR = norm(R(:)-R0(:));

D = eye(T) - diag(ones(T-1,1),-1);
constraint = M*D'-R-helper.reconstruct(W,H)+X;

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
    fprintf('L1_W=%f\n',sum(abs(W(:))));
    fprintf('L1_M=%f\n',sum(abs(M(:))));
    fprintf('L1_R=%f\n',sum(abs(R(:))));
    fprintf('dW=%f\n', dW);
    fprintf('dM=%f\n', dM);
    fprintf('dR=%f\n', dR);
end