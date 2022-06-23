function W = updateW(W0, H, X, params)
[N, K, L] = size(W0);
[~, T] = size(X);
opts = tfocs;
opts.maxIts = 500;
opts.tol = 1e-6;
if ~params.showPlot 
    opts.printEvery = 0;
end

% if params.lambda > 0
%     % Use smooth cross orthogonal panelty regularization
%     % Apply unconstrained split Bregman method to solve 
%     smoothkernel = ones(1,(2*L)-1);  % for factor competition
%     
%     B = zeros(K);
%     D = zeros(K);
%     Q = ones(K);
%     Q(1:K+1:end) = 0;   % off diagonal mask
%     
%     tol = 1e-3;
%     max_iter = 10;
%     for i=1:max_iter
%         % Step 1: Update W
%         op_recon = @(W, mode)tensor_conv_W(H, N, L, W, mode);
% %         op_recon_error = tfunc_scale(smooth_quad, 1, op_recon, -X);
%         op_reg = @(W, mode)smooth_cross_ortho_W(X, H, L, W, mode);
% %         op_reg_error = tfunc_scale(smooth_quad(params.lambda), 1, op_reg, B-D);
% %         op_smooth = tfunc_sum(op_recon_error, op_reg_error);
%         smoothF = {smooth_quad, smooth_quad(params.lambda)};
%         affineF = {op_recon, -X; op_reg, B-D};
%         W = tfocs(smoothF, affineF, [], W0, opts);
%         
%         % Plot to show progress
%         if params.showPlot 
%             Xhat = helper.reconstruct(W, H);
%             SimpleWHPlot(W, H, Xhat); 
%             drawnow
%         end
%         
%         if sqrt(sum((W-W0).^2, 'all')) < tol
%             fprintf('Step size tolerance of W reached\n')
%             break
%         end
%         
%         % Step 2: Update D
%         WTX = helper.transconv(W, X);
%         WTXS = conv2(WTX, smoothkernel, 'same');
%         WTXSHT = WTXS*H';
%         D = tfocs(smooth_quad(params.lambda), {1, -WTXSHT-B}, prox_l1(Q), D);
%         
%         % Step 3: Update B
%         B = B + WTXSHT - D;
%         W0 = W;
%     end
    
if params.lambdaL1W > 0
    % Solve LASSO
    op_recon = @(W, mode)tensor_conv_W(H, N, L, W, mode);
    W = tfocs(smooth_quad, {op_recon,-X}, prox_l1(params.lambdaL1W), W0, opts);
else
    % No regularization
    op_recon = @(W, mode)tensor_conv_W(H, N, L, W, mode);
    W = tfocs(smooth_quad, {op_recon,-X}, proj_Rplus, W0, opts);
end