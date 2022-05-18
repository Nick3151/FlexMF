function H = updateH(W, H0, X, params)
[N, K, L] = size(W);
[~, T] = size(X);
opts = tfocs;
opts.maxIts = 500;
if params.lambda > 0
    % Use smooth cross orthogonal panelty regularization
    % Apply unconstrained split Bregman method to solve 
    smoothkernel = ones(1,(2*L)-1);  % for factor competition
    
    B = zeros(K);
    D = zeros(K);
    Q = ones(K);
    Q(1:K+1:end) = 0;   % off diagonal mask
    
    tol = 1e-3;
    max_iter = 10;
    for i=1:max_iter
        % Step 1: Update H
        op_recon = @(H, mode)tensor_conv_H(W, T, H, mode);
        op_recon_error = tfunc_scale(smooth_quad, 1, op_recon, -X);
        op_reg = @(H, mode)smooth_cross_ortho_H(W, X, H, mode);
        op_reg_error = tfunc_scale(smooth_quad(params.lambda), 1, op_reg, B-D);
        op_smooth = tfunc_sum(op_recon_error, op_reg_error);
        H = tfocs(op_smooth, [], proj_Rplus, H0, opts);
        
        if sqrt(sum((H-H0).^2, 'all')) < tol
            break
        end
        
        % Step 2: Update D
        WTX = helper.transconv(W, X);
        WTXS = conv2(WTX, smoothkernel, 'same');
        WTXSHT = WTXS*H';
        D = tfocs(smooth_quad(params.lambda), {1, -WTXSHT-B}, prox_l1(Q), D, opts);
        
        % Step 3: Update B
        B = B + WTXSHT - D;
        H0 = H;
    end
    
elseif params.lambdaL1H > 0
    % Solve LASSO
    op_recon = @(H, mode)tensor_conv_H(W, T, H, mode);
    H = tfocs(smooth_quad, {op_recon,-X}, prox_l1pos(params.lambdaL1H), H0, opts);
else
    % No regularization
    op_recon = @(H, mode)tensor_conv_H(W, T, H, mode);
    H = tfocs(smooth_quad, {op_recon,-X}, [], H0, opts);
end