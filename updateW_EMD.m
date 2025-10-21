function [W, M, R, out] = updateW_EMD(W0, H, X, M0, R0, params)
% Update W with Earth-mover's distance(EMD) and smooth-orthogonal regularization

[N, K, L] = size(W0);
[~, T] = size(X);
opts_default = tfocs_SCD;
opts = opts_default;
opts.continuation = 1;
opts.tol = 1e-3;
opts.stopCrit = 4;
opts.maxIts = 500;
opts.alg = 'N83';
continue_opts = continuation();
% opts.debug = true;

if ~params.verbal
    opts.printEvery = 0;
    continue_opts.verbose = 0;
end

%% Initialization
W0_flat = reshape(W0, [N,K*L]);
X_pad = [zeros(N,L),X,zeros(N,L)];
H_pad = [zeros(K,L),H,zeros(K,L)];
M0_pad = [zeros(N,L),M0,zeros(N,L)];
R0_pad = [zeros(N,L),R0,zeros(N,L)];
% W0_ = [W0_flat, M0, R0];
W0_ = [W0_flat, M0_pad, R0_pad];

%% Linear operators
op_cross_orth_W = @(W, mode)smooth_cross_ortho_W(X_pad, H_pad, L, W, mode);
% op_cross_orth_W = @(W_, mode)cross_orth_EMD_W(X, H, L, W_, mode);
% op_cross_orth_W = @(W_, mode)cross_orth_EMD_W(X_pad, H_pad, L, W_, mode);
op_M = @(W_, mode)M_EMD_W(M0_pad, K, L, W_, mode);
op_R = @(W_, mode)R_EMD_W(R0_pad, K, L, W_, mode);
op_W = @(W_, mode)W_EMD_W(M0_pad, K, L, W_, mode);
op_constraint = @(W_, mode)constraint_EMD_W(H_pad, N, L, W_, mode);
op_recon = @(W, mode)tensor_conv_W(H_pad, N, L, W, mode);
op_TV = @(W, mode)total_variation_W(N, K, L, W, mode);

norm_cross_orth2 = linop_normest(op_cross_orth_W).^2;
norm_M2 = linop_normest(op_M).^2;
norm_R2 = linop_normest(op_R).^2;
norm_W2 = linop_normest(op_W).^2;
norm_constraint2 = linop_normest(op_constraint).^2;
norm_recon2 = linop_normest(op_recon).^2;

% proxScale_corss_orth = sqrt(norm_cross_orth2/norm_constraint2);
proxScale_corss_orth = sqrt(norm_cross_orth2/norm_recon2);
proxScale_M = sqrt(norm_M2/norm_constraint2);
proxScale_R = sqrt(norm_R2/norm_constraint2);
proxScale_W = sqrt(norm_W2/norm_constraint2);

%% Optimize with tfocs
lambda = params.lambda;
lambda_R = params.lambda_R;
lambda_M = params.lambda_M;
lambdaL1W = params.lambdaL1W;
lambda_TV = params.lambda_TV;
mu = 1e-3;

% divergence matrix
D = eye(T+2*L) - diag(ones(T+2*L-1,1),-1);
% Linear constraint
conjnegF = {proj_Rn};
% affineF = {op_constraint, X_pad};
affineF = {op_recon, -M0_pad*D'+R0_pad-X_pad};

% if lambda_M>0
%     affineF(end+1,:) = {linop_compose(op_M, 1/proxScale_M), 0};
%     conjnegF{end+1} = proj_linf(lambda_M*proxScale_M);
% end
% 
% if lambda_R>0
%     affineF(end+1,:) = {linop_compose(op_R, 1/proxScale_R), 0};
%     conjnegF{end+1} = proj_linf(lambda_R*proxScale_R);
% end

if lambda>0 && proxScale_corss_orth>0
    affineF(end+1,:) = {linop_compose(op_cross_orth_W, 1/proxScale_corss_orth), 0};
%     conjnegF{end+1} = proj_linf(lambda*proxScale_corss_orth);

    Q = ones(K);
    Q(1:K+1:end) = 0;   % off diagonal mask
    conjnegF{end+1} = proj_scale_linf(lambda*proxScale_corss_orth*Q);
end

if lambdaL1W>0 
    affineF(end+1,:) = {op_W, 0};
    conjnegF{end+1} = proj_linf(lambdaL1W);
end

if lambda_TV>0
    affineF(end+1,:) = {linop_compose(op_TV, op_W), 0};
    conjnegF{end+1} = proj_linf(lambda_TV);
end

% [W_, out] = tfocs_SCD(proj_Rplus_W(K*L), affineF, conjnegF, mu, W0_, [], opts, continue_opts);
[W, out] = tfocs_SCD(proj_Rplus, affineF, conjnegF, mu, W0, [], opts, continue_opts);

% Wflat = W_(:,1:K*L);
% W = reshape(Wflat, [N,K,L]);
% M = W_(:,L+K*L+(1:T));
% R = W_(:,3*L+K*L+T+(1:T));
M = M0;
R = R0;
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
%     fprintf('L1_W=%f\n',sum(abs(W(:))));
%     fprintf('L1_M=%f\n',sum(abs(M(:))));
%     fprintf('L1_R=%f\n',sum(abs(R(:))));
    fprintf('L1_W/X=%f\n',norm(W(:),1)/norm(X(:),1));
    fprintf('L1_M/X=%f\n',norm(M(:),1)/norm(X(:),1));
    fprintf('L1_R/X=%f\n',norm(R(:),1)/norm(X(:),1));
    fprintf('Constraint=%f\n', sum(constraint(:).^2/2))
%     fprintf('dW=%f\n', dW);
%     fprintf('dM=%f\n', dM);
%     fprintf('dR=%f\n', dR);
    fprintf('dW/X=%f\n', dW/norm(X(:)));
    fprintf('dM/X=%f\n', dM/norm(X(:)));
    fprintf('dR/X=%f\n', dR/norm(X(:)));
end