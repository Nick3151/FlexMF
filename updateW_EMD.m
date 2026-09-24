function [W, M, R, out] = updateW_EMD(W0, H, X, M0, ~, params)
% Update W with Earth-mover's distance(EMD) and smooth-orthogonal regularization
% The residual R = X + div(M) - conv(W,H) is eliminated from the solve and
% penalized as lambda_R*||R||_1, so the transport constraint holds exactly.

[N, K, L] = size(W0);
[~, T] = size(X);
opts_default = tfocs_SCD;
opts = opts_default;
opts.continuation = 1;
opts.tol = 1e-3;
opts.stopCrit = 4;
opts.maxIts = 200;
opts.alg = 'N83';
continue_opts = continuation();
% opts.debug = true;

if isfield(params, 'muDecrement')
    continue_opts.muDecrement = params.muDecrement;
end

if ~params.verbal
    opts.printEvery = 0;
    continue_opts.verbose = 0;
end

%% Initialization
Xcorr = helper.correct_warp(X,M0);
W0_flat = reshape(W0, [N,K*L]);
X_pad = [zeros(N,L),Xcorr,zeros(N,L)];
H_pad = [zeros(K,L),H,zeros(K,L)];
W0_ = [W0_flat, M0];

%% Linear operators
op_cross_orth_W = @(W_, mode)cross_orth_EMD_W(X_pad, H_pad, L, W_, mode);
op_M = @(W_, mode)M_EMD_W(M0, K, L, W_, mode);
op_W = @(W_, mode)W_EMD_W(M0, K, L, W_, mode);
op_fit = @(W_, mode)fit_EMD_W(H, N, L, W_, mode);
op_TV = @(W, mode)total_variation_W(N, K, L, W, mode);

norm_cross_orth2 = linop_normest(op_cross_orth_W).^2;
norm_M2 = linop_normest(op_M).^2;
norm_W2 = linop_normest(op_W).^2;
norm_TV2 = linop_normest(op_TV).^2;
norm_fit2 = linop_normest(op_fit).^2;

proxScale_cross_orth = sqrt(norm_cross_orth2/norm_fit2);
proxScale_M = sqrt(norm_M2/norm_fit2);
proxScale_W = sqrt(norm_W2/norm_fit2);
proxScale_TV = sqrt(norm_TV2/norm_fit2);

%% Optimize with tfocs
lambda = params.lambda;
lambda_R = params.lambda_R;
lambda_M = params.lambda_M;
lambdaL1W = params.lambdaL1W;
lambda_TV = params.lambda_TV;
if isfield(params, 'mu') && ~isempty(params.mu)
    mu = params.mu;
else
    mu = 1e-1;
end
assert(lambda_R > 0, 'updateW_EMD requires lambda_R > 0.');

% lambda_R*||X + div(M) - conv(W,H)||_1
affineF = {op_fit, X};
conjnegF = {proj_linf(lambda_R)};
dualNames = {'fit'};
dualScales = 1;

if lambda_M>0
    % Homotopy: linearly ramp lambda_M from lambda_M/homotopy to lambda_M
    % over the first homotopy iterations (homotopy=0 disables ramp)
    nHomotopy = params.homotopy;
    if nHomotopy > 0 && params.currentiter <= nHomotopy
        lambda_M_eff = lambda_M * params.currentiter / nHomotopy;
    else
        lambda_M_eff = lambda_M;
    end
    affineF(end+1,:) = {linop_compose(op_M, 1/proxScale_M), 0};
    conjnegF{end+1} = proj_linf(lambda_M_eff*proxScale_M);
    dualNames{end+1} = 'M';
    dualScales(end+1) = proxScale_M;
end

if lambda>0 && proxScale_cross_orth>0
    affineF(end+1,:) = {linop_compose(op_cross_orth_W, 1/proxScale_cross_orth), 0};
    conjnegF{end+1} = proj_linf(lambda*proxScale_cross_orth);
    dualNames{end+1} = 'cross_W';
    dualScales(end+1) = proxScale_cross_orth;
end

if lambdaL1W>0 
    affineF(end+1,:) = {linop_compose(op_W, 1/proxScale_W), 0};
    conjnegF{end+1} = proj_linf(lambdaL1W*proxScale_W);
    dualNames{end+1} = 'L1_W';
    dualScales(end+1) = proxScale_W;
end

if lambda_TV>0
    affineF(end+1,:) = {linop_compose(op_TV, op_W, 1/(proxScale_TV*proxScale_W)), 0};
    conjnegF{end+1} = proj_linf(lambda_TV*proxScale_TV*proxScale_W);
    dualNames{end+1} = 'TV_W';
    dualScales(end+1) = proxScale_TV*proxScale_W;
end

if isfield(params, 'dual')
    dual0 = params.dual;
else
    dual0 = struct();
end
z0 = helper.dual_warm_start(dual0, dualNames, dualScales, affineF);

[W_, out] = tfocs_SCD(proj_Rplus_W(K*L), affineF, conjnegF, mu, W0_, z0, opts, continue_opts);
out.dual_unscaled = helper.dual_unscale(dual0, out.dual, dualNames, dualScales);

Wflat = W_(:,1:K*L);
W = reshape(Wflat, [N,K,L]);
M = W_(:,K*L+(1:T));
R = helper.correct_warp(X,M) - helper.reconstruct(W,H);
dW = norm(W(:)-W0(:));
dM = norm(M(:)-M0(:));

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
    fprintf('L1_W/X=%f\n',norm(W(:),1)/norm(X(:),1));
    fprintf('L1_M/X=%f\n',norm(M(:),1)/norm(X(:),1));
    fprintf('L1_R/X=%f\n',norm(R(:),1)/norm(X(:),1));
    fprintf('dW/X=%f\n', dW/norm(X(:)));
    fprintf('dM/X=%f\n', dM/norm(X(:)));
end
