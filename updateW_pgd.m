function W = updateW_pgd(W0, H, X, params, E, C)
[N, K, L] = size(W0);
[~, T] = size(X);

X_pad = [zeros(N,L), X, zeros(N,L)];
H_pad = [zeros(K,L), H, zeros(K,L)];

if params.neg_prop==0 && params.lambda > 0

    smoothkernel = ones(1,(2*L)-1);
    B  = zeros(K);
    D  = zeros(K);
    Q  = ones(K);
    Q(1:K+1:end) = 0;   % off-diagonal mask

    % TV auxiliary variables — init if not passed in
    if nargin < 5 || isempty(E)
        E = zeros(N, K, L);
        C = zeros(N, K, L);
    end

    tol_W    = 1e-2;
    max_iter = 10;
    W        = W0;

    for i = 1:max_iter

        %% --- Step 1: Update W via projected gradient descent ---
        for inner = 1:params.inner_iter

            % (a) Gradient of reconstruction term: d/dW 1/2||W⊛H - X||^2
            Xhat      = helper.reconstruct(W, H_pad);               % N x (T+2L)
            R_recon   = Xhat - X_pad;                               % N x (T+2L)
            grad_data = tensor_conv_grad_W(H_pad, N, L, R_recon);   % N x K x L

            % (b) Gradient of cross-orth coupling: d/dW alpha/2||(W⊛^T X)SH^T - D + B||^2
            WTX       = helper.transconv(W, X_pad);                 % K x K
            WTXS      = conv2(WTX, smoothkernel, 'same');           % K x K
            WTXSHT    = WTXS * H_pad';                              % K x K
            R_co      = WTXSHT - D + B;                             % K x K
            grad_co   = transconv_grad_W(X_pad, H_pad, ...
                            smoothkernel, R_co, N, K, L);           % N x K x L

            % (c) Gradient of TV coupling: d/dW beta/2||∇_l W - E + C||^2
            R_tv      = grad_l(W) - E + C;                         % N x K x L
            grad_tv   = grad_l_adj(R_tv);                          % N x K x L

            % (d) Total gradient
            grad_W = grad_data + params.alpha_W * grad_co + params.beta * grad_tv;

            % (e) Gradient step + project onto W >= 0
            W = max(W - params.lr * grad_W, 0);
        end

        % Plot to show progress
        Xhat = helper.reconstruct(W, H_pad);
        if params.showPlot
            SimpleWHPlot_patch(W, H_pad, 'Data', Xhat);
            drawnow
        end

        % Check W convergence
        dW = sqrt(sum((W - W0).^2, 'all') / K);
        if dW < tol_W
            if params.verbal
                fprintf('Step size tolerance of W reached\n')
            end
            break
        end

        %% --- Step 2: Update D (unchanged from original) ---
        WTX    = helper.transconv(W, X_pad);
        WTXS   = conv2(WTX, smoothkernel, 'same');
        WTXSHT = WTXS * H_pad';

        % soft threshold off-diagonal entries
        D = sign(WTXSHT + B) .* max(abs(WTXSHT + B) - params.lambda * Q / params.alpha_W, 0);

        %% --- Step 3: Update B (unchanged from original) ---
        B = B + WTXSHT - D;

        %% --- Step 4: Update E via soft threshold ---
        E = shrink(grad_l(W) + C, params.lambda_TV / params.beta);

        %% --- Step 5: Update C ---
        C = C + grad_l(W) - E;

        if params.verbal
            fprintf('dW          = %f\n', dW);
            fprintf('reg         = %.4f\n', norm(Q(:).* WTXSHT(:), 1));
            fprintf('recon       = %.4f\n', sum((X_pad(:) - Xhat(:)).^2) / 2);
            fprintf('threshold   = %.4f\n', params.lambda / params.alpha_W);
            fprintf('norm WTXSHT = %.4f\n', norm(WTXSHT, 'fro'));
            fprintf('norm D      = %.4f\n', norm(Q(:).*D(:), 1));
            fprintf('norm B      = %.4f\n', norm(Q(:).*B(:), 1));
            fprintf('\n')
        end

        W0 = W;
    end

elseif params.lambdaL1W > 0
    % LASSO — unchanged
    op_recon = @(W, mode)tensor_conv_W(H_pad, N, L, W, mode);
    opts = get_tfocs_opts(params);
    W = tfocs(smooth_quad, {op_recon, -X_pad}, prox_l1(params.lambdaL1W), W0, opts);

else
    % No regularization — unchanged
    op_recon = @(W, mode)tensor_conv_W(H_pad, N, L, W, mode);
    opts = get_tfocs_opts(params);
    W = tfocs(smooth_quad, {op_recon, -X_pad}, proj_Rn, W0, opts);
end

end % function


%% =========================================================
%  Helper functions
%% =========================================================

function dW = tensor_conv_grad_W(H_pad, N, L, R)
% Gradient of 1/2||W⊛H - X||^2 w.r.t. W
% = correlation of residual R with H_pad, sliced per lag
% R: N x (T+2L),  H_pad: K x (T+2L)
% Output: N x K x L
[K, T_pad] = size(H_pad);
dW = zeros(N, K, L);
for l = 1:L
    % For lag l, W(:,k,l) contributes to R via a shift of H by l-1
    H_shift = H_pad(:, l:l+T_pad-L);          % K x (T_pad-L+1) — align with R
    dW(:,:,l) = R(:, 1:size(H_shift,2)) * H_shift';
end
end


function grad = transconv_grad_W(X_pad, H_pad, smoothkernel, R_co, N, K, L)
% Gradient of alpha/2||(W⊛^T X)SH^T - D + B||^2 w.r.t. W
% R_co = (W⊛^T X)SH^T - D + B,  size K x K
% This is the adjoint of the transconv operator applied to R_co * H * S^T
% Output: N x K x L
[~, T_pad] = size(H_pad);
RS   = conv2(R_co, smoothkernel, 'same');      % K x K  (undo the S smoothing)
RSH  = RS * H_pad;                             % K x (T_pad)
grad = zeros(N, K, L);
for l = 1:L
    % adjoint of transconv at lag l: correlate X_pad with RSH shifted by l
    X_shift = X_pad(:, l:l+T_pad-L);
    grad(:,:,l) = X_shift * RSH(:, 1:size(X_shift,2))';
end
end


function dW = grad_l(W)
% Forward finite difference along L (dim 3)
% W: N x K x L  ->  dW: N x K x L  (zero-padded at boundary)
dW = diff(W, 1, 3);
dW = cat(3, dW, zeros(size(W,1), size(W,2), 1));
end


function g = grad_l_adj(V)
% Adjoint of grad_l (negative backward difference along dim 3)
% V: N x K x L  ->  g: N x K x L
g = -diff(cat(3, zeros(size(V,1), size(V,2), 1), V), 1, 3);
end


function out = shrink(x, tau)
% Element-wise soft threshold
out = sign(x) .* max(abs(x) - tau, 0);
end


function opts = get_tfocs_opts(params)
opts = tfocs;
opts.maxIts    = 1000;
opts.tol       = 1e-4;
opts.restart   = 50;
opts.alg       = 'N83';
if ~params.showPlot
    opts.printEvery = 0;
end
end