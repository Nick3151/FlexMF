function [H, M, R, out] = updateH_EMD(W, H0, X, M0, R0, params)
% Update H with Earth-mover's distance(EMD) and smooth-orthogonal regularization

[N, K, L] = size(W);
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
H0_ = [H0; M0; R0];

%% Linear operators
smoothkernel = ones(1,(2*L)-1);  % for factor competition
WTX = helper.transconv(W, X);
WTXS = conv2(abs(WTX), smoothkernel, 'same');
A = WTXS;
op_cross_orth_H = @(H_, mode)cross_orth_EMD_H(A, N, H_, mode);
op_M = @(H_, mode)M_EMD_H(M0, K, H_, mode);
op_R = @(H_, mode)R_EMD_H(R0, K, H_, mode);
op_H = @(H_, mode)H_EMD_H(M0, K, H_, mode);
lambda_R = params.lambda_R;
op_obj = @(H_, mode)obj_EMD_H(M0, K, lambda_R, H_, mode);
op_constraint = @(H_, mode)constraint_EMD_H(W, T, H_, mode);


norm_cross_orth2 = linop_normest(op_cross_orth_H).^2;
norm_M2 = linop_normest(op_M).^2;
norm_R2 = linop_normest(op_R).^2;
norm_H2 = linop_normest(op_H).^2;
norm_constraint2 = linop_normest(op_constraint).^2;

proxScale_corss_orth = sqrt(norm_cross_orth2/norm_constraint2);
proxScale_M = sqrt(norm_M2/norm_constraint2);
proxScale_R = sqrt(norm_R2/norm_constraint2);
proxScale_H = sqrt(norm_H2/norm_constraint2);
%% Optimize with tfocs
lambda = params.lambda;
lambda_R = params.lambda_R;
lambda_M = params.lambda_M;
lambdaL1H = params.lambdaL1H;
Reweight = params.Reweight;
mu = 1e-3;

% affineF = {linop_compose(op_M, 1/proxScale_M), 0; ...
%            linop_compose(op_R, 1/proxScale_R), 0; ...
%            op_constraint, X};
% conjnegF = {proj_linf(proxScale_M), proj_linf(lambda_R*proxScale_R), proj_Rn};
% 
% if lambda>0
%     affineF(end+1,:) = {linop_compose(op_cross_orth_H, 1/proxScale_corss_orth), 0};
%     conjnegF{end+1} = proj_linf(lambda*proxScale_corss_orth);
% end
% 
% if lambdaL1H>0
%     affineF(end+1,:) = {linop_compose(op_H, 1/proxScale_H), 0};
%     conjnegF{end+1} = proj_linf(lambdaL1H*proxScale_H);
% end
% 
% [H_, out] = tfocs_SCD(proj_Rplus_H(K), affineF, conjnegF, mu, H0_, [], opts);

conjnegF = {proj_Rn, proj_linf(lambda_R), proj_linf(lambda_M)};
affineF = {op_constraint, X; ...
           op_R, 0; ...
           op_M, 0 };

if lambda>0
    affineF(end+1,:) = {op_cross_orth_H, 0};
    conjnegF{end+1} = proj_linf(lambda);
end

smooth_win = 10;
[K,T] = size(H0);
H0_smooth = zeros(K,T);
for k=1:K
    H0_smooth(k,:) = filtfilt(ones(1,smooth_win)/smooth_win, 1, H0(k,:));
end

if lambdaL1H>0
    affineF(end+1,:) = {op_H, 0};
    if Reweight && (params.currentiter > 10)
        epsilon = 1e-2;
        conjnegF{end+1} = proj_abs_box(lambdaL1H./(abs(H0_smooth)+epsilon));
    else
        conjnegF{end+1} = proj_linf(lambdaL1H);
    end
end

[H_, out] = tfocs_SCD(proj_Rplus_H(K), affineF, conjnegF, mu, H0_, [], opts, continue_opts);

% [H_, out] = solver_sBPDN_W(op_constraint,op_obj,-X,0,mu,[],[],opts);


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
    fprintf('L1_H=%f\n',sum(abs(H(:))));
    fprintf('L1_M=%f\n',sum(abs(M(:))));
    fprintf('L1_R=%f\n',sum(abs(R(:))));
    fprintf('Constraint=%f\n', sum(constraint(:).^2/2))
    fprintf('dH=%f\n', dH);
    fprintf('dM=%f\n', dM);
    fprintf('dR=%f\n', dR);
end