function H = updateH(W, H0, X, params)
[N, K, L] = size(W);
[~, T] = size(X);
opts_default = tfocs;
opts = opts_default;
opts.maxIts = 500;
opts.tol = 1e-6;
opts.restart = 50;
if ~params.showPlot 
    opts_default.printEvery = 0;
    opts.printEvery = 0;
end
% ops.alg = 'N83';

if params.lambda > 0
    % Use smooth cross orthogonal panelty regularization
    % Apply unconstrained split Bregman method to solve 
    smoothkernel = ones(1,(2*L)-1);  % for factor competition
    WTX = helper.transconv(W, X);
    WTXS = conv2(abs(WTX), smoothkernel, 'same');
    if params.lambdaL1H > 0
        B = zeros(K+T,K);
        D = zeros(K+T,K);
        Q = ones(K+T,K);
        for k = 1:K
            Q(k,k) = 0;   % off diagonal mask
        end
        A = [WTXS; eye(T)];
    else
        B = zeros(K);
        D = zeros(K);
        Q = ones(K);
        Q(1:K+1:end) = 0;   % off diagonal mask
        A = WTXS;
    end
    
    tol_H = 1e-3;
    max_iter = 20;
    for i=1:max_iter
        % Step 1: Update H
        op_recon = @(H, mode)tensor_conv_H(W, T, H, mode);
%         op_recon_error = tfunc_scale(smooth_quad, 1, op_recon, -X);
        op_reg = @(H, mode)smooth_cross_ortho_H(A, K, H, mode);
%         op_reg_error = tfunc_scale(smooth_quad(params.lambda), 1, op_reg, B-D);
%         op_smooth = tfunc_sum(op_recon_error, op_reg_error);
%         op_smooth = @(varargin)off_diag_frob_norm_sqr(params.lambda, Q, varargin{:});
        smoothF = {smooth_quad, smooth_quad(params.alpha)};
        affineF = {op_recon, -X; op_reg, B-D};
        H = tfocs(smoothF, affineF, proj_Rplus, H0, opts);
        
        % Plot to show progress
        if params.showPlot 
            Xhat = helper.reconstruct(W, H);
            SimpleWHPlot(W, H, Xhat); 
            drawnow
        end
        
        dH = sqrt(mean((H(:)-H0(:)).^2));
        if dH < tol_H
            fprintf('Step size tolerance of H reached\n')
            break
        end
        
        % Step 2: Update D
        AH = A*H';
        D = tfocs(smooth_quad(params.alpha), {1, -AH-B}, prox_l1(params.lambda*Q), D, opts_default);
        
        % Step 3: Update B
        B = B + AH - D;
        
        if params.showPlot 
            fprintf('dH=%f\n',dH);
            fprintf('reg=%f\n',sum(Q(:).*AH(:)));
            fprintf('D=%f\n',sum(Q(:).*D(:)));
            fprintf('B=%f\n',sum(Q(:).*B(:)));
        end
        
        H0 = H;
    end
%     op_recon = @(H, mode)tensor_conv_H(W, T, H, mode);
%     op_reg = @(H, mode)smooth_cross_ortho_H(W, X, H, mode);
%     epsilon = 100;
%     opts.nonneg = true;
%     H = solver_sBPDN_W(op_recon, op_reg, X, epsilon, 0, H0, [], opts);

elseif params.lambdaL1H > 0
    % Solve LASSO
    op_recon = @(H, mode)tensor_conv_H(W, T, H, mode);
    H = tfocs(smooth_quad, {op_recon,-X}, prox_l1pos(params.lambdaL1H), H0, opts_default);
else
    % No regularization
    op_recon = @(H, mode)tensor_conv_H(W, T, H, mode);
    H = tfocs(smooth_quad, {op_recon,-X}, proj_Rplus, H0, opts_default);
end