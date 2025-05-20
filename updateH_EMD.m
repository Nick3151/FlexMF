function [H, M, R, out] = updateH_EMD(W, H0, X, M0, R0, params)
% Update H with Earth-mover's distance(EMD) and smooth-orthogonal regularization

[N, K, L] = size(W);
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
H0_ = [H0; M0; R0];

%% Linear operators
smoothkernel = ones(1,(2*L)-1);  % for factor competition
WTX = helper.transconv(W, X);
WTXS = conv2(abs(WTX), smoothkernel, 'same');
A = WTXS;
op_f1 = @(H_, mode)f1_EMD_H(A, N, H_, mode);
op_f2 = @(H_, mode)f2_EMD_H(M0, K, H_, mode);
op_f3 = @(H_, mode)f3_EMD_H(R0, K, H_, mode);
op_f4 = @(H_, mode)f4_EMD_H(W, T, H_, mode);


normf12 = linop_normest(op_f1).^2;
normf22 = linop_normest(op_f2).^2;
normf32 = linop_normest(op_f3).^2;
normf42 = linop_normest(op_f4).^2;

proxScale1 = sqrt(normf12/normf42);
proxScale2 = sqrt(normf22/normf42);
proxScale3 = sqrt(normf32/normf42);
%% Optimize with tfocs
lambda = params.lambda;
lambda_R = params.lambda_R;
mu = 1e-3;

if lambda>0
    [H_, out] = tfocs_SCD(proj_Rplus_H(K), {linop_compose(op_f1, 1/proxScale1), 0; ...
        linop_compose(op_f2, 1/proxScale2), 0; linop_compose(op_f3, 1/proxScale3), 0; op_f4, X}, ...
        {proj_linf(lambda*proxScale1), proj_linf(proxScale2), proj_linf(lambda_R*proxScale3), proj_Rn}, mu, H0_, [], opts);
else
    [H_, out] = tfocs_SCD(proj_Rplus_H(K), {linop_compose(op_f2, 1/proxScale2), 0; ...
        linop_compose(op_f3, 1/proxScale3), 0; op_f4, X}, ...
        {proj_linf(proxScale2), proj_linf(lambda_R*proxScale3), proj_Rn}, mu, H0_, [], opts);
end

H = H_(1:K,:);
M = H_(K+(1:N),:);
R = H_(K+N+(1:N),:);
dH = norm(H(:)-H0(:));
dM = norm(M(:)-M0(:));
dR = norm(R(:)-R0(:));

D = eye(T) - diag(ones(T-1,1),-1);
constraint = M*D'-R-helper.reconstruct(W,H)+X;

%% Print intermediate results
if params.verbal
    AH = A*H';
    Q = ones(K);
    Q(1:K+1:end) = 0;
    Xhat = helper.reconstruct(W, H);
    fprintf('reg=%f\n',sum(Q(:).*AH(:)));
    fprintf('recon=%f\n',sum((X(:)-Xhat(:)).^2)/2);
    fprintf('L1_M=%f\n',sum(abs(M(:))));
    fprintf('L1_R=%f\n',sum(abs(R(:))));
    fprintf('Constraint=%f\n', sum(constraint(:).^2/2))
    fprintf('dH=%f\n', dH);
    fprintf('dM=%f\n', dM);
    fprintf('dR=%f\n', dR);
end