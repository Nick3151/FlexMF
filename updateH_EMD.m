function [H, M, R, out] = updateH_EMD(W, H0, X, M0, ~, params)
% Update H with Earth-mover's distance(EMD) and smooth-orthogonal regularization
% The residual R = X + div(M) - conv(W,H) is eliminated from the solve and
% penalized as lambda_R*||R||_1, so the transport constraint holds exactly.

[N, K, L] = size(W);
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
H0_ = [H0; M0];

%% Linear operators
Xcorr = helper.correct_warp(X,M0);    % Correct data with warping/jitering
smoothkernel = ones(1,(2*L)-1);  % for factor competition
WTX = helper.transconv(W, Xcorr);
WTXS = conv2(abs(WTX), smoothkernel, 'same');
A = WTXS;
op_cross_orth_H = @(H_, mode)cross_orth_EMD_H(A, N, H_, mode);
op_M = @(H_, mode)M_EMD_H(M0, K, H_, mode);
op_H = @(H_, mode)H_EMD_H(M0, K, H_, mode);
op_fit = @(H_, mode)fit_EMD_H(W, T, H_, mode);

norm_cross_orth2 = linop_normest(op_cross_orth_H).^2;
norm_M2 = linop_normest(op_M).^2;
norm_H2 = linop_normest(op_H).^2;
norm_fit2 = linop_normest(op_fit).^2;

proxScale_cross_orth = sqrt(norm_cross_orth2/norm_fit2);
proxScale_M = sqrt(norm_M2/norm_fit2);
proxScale_H = sqrt(norm_H2/norm_fit2);
%% Optimize with tfocs
lambda = params.lambda;
lambda_R = params.lambda_R;
lambda_M = params.lambda_M;
lambdaL1H = params.lambdaL1H;
Reweight = params.Reweight;
if isfield(params, 'mu') && ~isempty(params.mu)
    mu = params.mu;
else
    mu = 1e-1;
end
assert(lambda_R > 0, 'updateH_EMD requires lambda_R > 0.');

% lambda_R*||X + div(M) - conv(W,H)||_1
affineF = {op_fit, X};
conjnegF = {proj_linf(lambda_R)};

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
end

if lambda>0 && proxScale_cross_orth>0
    affineF(end+1,:) = {linop_compose(op_cross_orth_H, 1/proxScale_cross_orth), 0};
    conjnegF{end+1} = proj_linf(lambda*proxScale_cross_orth);
end

if lambdaL1H>0
    affineF(end+1,:) = {linop_compose(op_H, 1/proxScale_H), 0};
    % IRL1 from iter 2 onward (after FlexMF row-normalizes H); iter 1 uses uniform L1
    if Reweight && params.currentiter > 1
        epsilon = 1e-2;
        conjnegF{end+1} = proj_abs_box(lambdaL1H./(abs(H0)+epsilon)*proxScale_H);
    else
        conjnegF{end+1} = proj_linf(lambdaL1H*proxScale_H);
    end
end

[H_, out] = tfocs_SCD(proj_Rplus_H(K), affineF, conjnegF, mu, H0_, [], opts, continue_opts);

H = H_(1:K,:);
M = H_(K+(1:N),:);
R = helper.correct_warp(X,M) - helper.reconstruct(W,H);
dH = norm(H(:)-H0(:));
dM = norm(M(:)-M0(:));

%% Print intermediate results
if params.verbal
    AH = A*H';
    Q = ones(K);
    Q(1:K+1:end) = 0;
    Xhat = helper.reconstruct(W, H);
    fprintf('reg=%f\n',sum(Q(:).*AH(:)));
    fprintf('recon=%f\n',sum((X(:)-Xhat(:)).^2)/2);
    fprintf('L1_H/X=%f\n',norm(H(:),1)/norm(X(:),1));
    fprintf('L1_M/X=%f\n',norm(M(:),1)/norm(X(:),1));
    fprintf('L1_R/X=%f\n',norm(R(:),1)/norm(X(:),1));
    fprintf('dH/X=%f\n', dH/norm(X(:)));
    fprintf('dM/X=%f\n', dM/norm(X(:)));
end
