function [W, H, errors,grads,loadings,power] = FlexMF(X, varargin)
%
% USAGE: 
%
% [W, H, cost, loadings, power] = seqNMF(X, ...    % X is the data matrix
%       'K', 10, 'L', 20, 'lambda', .1, ...        % Other inputs optional
%       'W_init', W_init, 'H_init', H_init, ...
%       'showPlot', 1, 'maxiter', 20, 'tolerance', -Inf, 'shift', 1, ... 
%       'lambdaL1W', 0, 'lambdaL1H', 0, ...
%       'lambdaOrthoH', 0, 'lambdaOrthoW', 0, 'M', M)
%
% ------------------------------------------------------------------------
% DESCRIPTION:
%
%   Factorizes the NxT data matrix X into K factors 
%   Factor exemplars are returned in the NxKxL tensor W
%   Factor timecourses are returned in the KxT matrix H
%
%                                    ----------    
%                                L  /         /|
%                                  /         / |
%        ----------------         /---------/  |          ----------------
%        |              |         |         |  |          |              |
%      N |      X       |   =   N |    W    |  /   (*)  K |      H       |           
%        |              |         |         | /           |              |
%        ----------------         /----------/            ----------------
%               T                      K                         T
% See paper: 
%   XXXXXXXXXXXXXXXXX
%
% ------------------------------------------------------------------------
%
% INPUTS:
%
% Name              Default                             Description
%  X                                                    Data matrix (NxT) to factorize
% 'K'               10                                  Number of factors
% 'L'               100                                 Length (timebins) of each factor exemplar
% 'lambda'          .001                                Regularization parameter
% 'W_init'          max(X(:))*rand(N,K,L)               Initial W
% 'H_init'          max(X(:))*rand(K,T)./(sqrt(T/3))    Initial H (rows have norm ~1 if max(data) is 1)
% 'showPlot'        1                                   Plot every iteration? no=0
% 'maxiter'         100                                 Maximum # iterations to run
% 'tolerance'       -Inf                                Stop if improved less than this;  Set to -Inf to always run maxiter
% 'shift'           1                                   Shift factors to center; Helps avoid local minima
% 'lambdaL1W'       0                                   L1 sparsity parameter; Increase to make W's more sparse
% 'lambdaL1H'       0                                   L1 sparsity parameter; Increase to make H's more sparse
% 'W_fixed'         0                                   Fix W during the fitting proceedure   
% 'SortFactors'     1                                   Sort factors by loadings
% 'lambdaOrthoH'    0                                   ||HSH^T||_1,i~=j; Encourages events-based factorizations
% 'lambdaOrthoW'    0                                   ||Wflat^TWflat||_1,i~=j; ; Encourages parts-based factorizations
% 'useWupdate'      1                                   Wupdate for cross orthogonality often doesn't change results much, and can be slow, so option to remove  
% 'M'               ones(N,T)                           Masking matrix if excluding a random test set from the fit
% ------------------------------------------------------------------------
% OUTPUTS:
%
% W                         NxKxL tensor containing factor exemplars
% H                         KxT matrix containing factor timecourses
% errors                      1x(#Iterations+1) vector containing 
%                               reconstruction error at each iteration. 
%                               cost(1) is error before 1st iteration.
% loadings                  1xK vector containing loading of each factor 
%                               (Fraction power in data explained by each factor)
% power                     Fraction power in data explained 
%                               by whole reconstruction
%
%                           Note, if doing fit with masked (held-out) data,
%                               the cost and power do not include masked
%                               (M==0) test set elements
% ------------------------------------------------------------------------
% CREDITS:
%   Emily Mackevicius and Andrew Bahle, 2/1/2018
%
%   Original CNMF algorithm: Paris Smaragdis 2004
%   (https://link.springer.com/chapter/10.1007/978-3-540-30110-3_63)
%   Adapted from NMF toolbox by Colin Vaz 2015 (http://sail.usc.edu)
%
%   Please cite our paper: 
%       https://www.biorxiv.org/content/early/2018/03/02/273128
%% parse function inputs

% Parse inputs
[X,N,T,K,L,params] = parse_seqNMF_params(X, varargin);

%% initialize
W = params.W_init;
H = params.H_init;

Xhat = helper.reconstruct(W, H); 
mask = find(params.M == 0); % find masked (held-out) indices 
X(mask) = Xhat(mask); % replace data at masked elements with reconstruction, so masked datapoints do not effect fit

smoothkernel = ones(1,(2*L)-1);  % for factor competition
smallnum = max(X(:))*1e-6; 
lasttime = 0;

% Calculate reconstruction error 
errors = zeros(params.maxiter+1, 1);
errors(1) = norm(X-Xhat, 'fro');

W0 = zeros(N,K,L);
grads_H_recon_all = [];
grads_H_WXH_all = [];
grads_H_HH_all = [];
grads_H_L1H_all = [];
grads_W_recon_all = [];
grads_W_WXH_all = [];
grads_W_WW_all = [];
grads_W_L1W_all = [];
etaH_all = [];
etaW_all = [];
for iter = 1 : params.maxiter
    tic
    fprintf('Iter %d\n', iter);
    % Compute terms for standard CNMF H update 
    WTX = helper.transconv(W, X);
    WTXhat = helper.transconv(W, Xhat);
         
    % Compute regularization terms for H update
    if params.lambda>0
        WTXS = conv2(WTX, smoothkernel, 'same');
        WTXSHT = WTXS*H';
        dRdH = params.lambda.*(~eye(K)).*sign(WTXSHT)*WTXS;  
    else 
        dRdH = 0; 
    end
    if params.lambdaOrthoH>0
        HS = conv2(H, smoothkernel, 'same');
        HSHT = HS*H';
        dHHdH = params.lambdaOrthoH*(~eye(K)).*sign(HSHT)*HS;
    else
        dHHdH = 0;
    end
    if params.lambdaL1H>0
        dL1HdH = params.lambdaL1H*sign(H);
    else
        dL1HdH = 0;
    end
    dLdH = WTXhat - WTX + dRdH + dL1HdH + dHHdH; % include L1 sparsity, if specified
    grad_H_recon = mean(abs(WTXhat - WTX),'all');
    grad_H_WXH = mean(abs(dRdH),'all');
    grad_H_HH = mean(abs(dHHdH),'all');
    grad_H_L1H = mean(abs(dL1HdH),'all');
    fprintf('Reconstruction gradient over H: %f\n', grad_H_recon);
    fprintf('Regularization WXH gradient over H: %f\n', grad_H_WXH);
    fprintf('Regularization HH gradient over H: %f\n', grad_H_HH);
    fprintf('Regularization L1H gradient over H: %f\n', grad_H_L1H);
    grads_H_recon_all = [grads_H_recon_all grad_H_recon];
    grads_H_WXH_all = [grads_H_WXH_all grad_H_WXH];
    grads_H_HH_all = [grads_H_HH_all grad_H_HH];
    grads_H_L1H_all = [grads_H_L1H_all grad_H_L1H];
    
    % Update H using line search
%     disp('line search on H')
%     tic
    alpha = 0.4;
    beta = 0.6;
    etaH = 1;
    while helper.total_loss(W,H-etaH*dLdH,X,params) > helper.total_loss(W,H,X,params)-alpha*norm(dLdH,'fro')^2*etaH
        etaH = beta*etaH;
    end
%     toc
    H = H - etaH*dLdH;
    H(H<0) = 0;
    fprintf('Line search step on H: %f\n', etaH);
    etaH_all = [etaH_all etaH];
    
    if ~params.W_fixed
    % Update each W_l separately
        grad_W_recon = [];
        grad_W_WXH = [];
        grad_W_WW = [];
        grad_W_L1W = [];
        etaW_tmp = [];
        
        Xhat = helper.reconstruct(W, H); 
        mask = find(params.M == 0); % find masked (held-out) indices 
        X(mask) = Xhat(mask); % replace data at masked elements with reconstruction, so masked datapoints do not effect fit
        
        WTXSHT = WTXS*H';
        if params.lambdaOrthoW>0
            Wflat = sum(W,3);
        end
        if params.lambda>0 && params.useWupdate
            XS = conv2(X, smoothkernel, 'same'); 
        end
        if params.lambdaL1W>0
            dL1WdW = params.lambdaL1W*sign(W);
        else
            dL1WdW = 0;
        end
        for l = 1 : L % could parallelize to speed up for long L
            % Compute terms for standard CNMF W update
            H_shifted = circshift(H,[0,l-1]);
            XHT = X * H_shifted';
            XhatHT = Xhat * H_shifted';

            % Compute regularization terms for W update
            if params.lambda>0 && params.useWupdate % Often get similar results with just H update, so option to skip W update
                dRdW = params.lambda.*XS*(H_shifted')*((~eye(K)).*sign(WTXSHT)); 
            else
                dRdW = 0;
            end
            if params.lambdaOrthoW>0
                dWWdW = params.lambdaOrthoW*Wflat*(~eye(K)).*sign(Wflat'*Wflat);
            else
                dWWdW = 0;
            end
            
            dLdW = XhatHT - XHT + dRdW + squeeze(dL1WdW(:,:,l)) + dWWdW; % include L1 and Worthogonality sparsity, if specified
            grad_W_recon = [grad_W_recon mean(abs(XhatHT - XHT), 'all')];
            grad_W_WXH = [grad_W_WXH mean(abs(dRdW), 'all')];
            grad_W_WW = [grad_W_WW mean(abs(dWWdW), 'all')];
            grad_W_L1W = [grad_W_L1W mean(abs(squeeze(dL1WdW(:,:,l))), 'all')];
            
            % Update W using line search
%             disp('line search on W')
%             tic
            alpha = 0.4;
            beta = 0.6;
            etaW = 1;
            while helper.total_loss(W-etaW*dLdW,H,X,params) > helper.total_loss(W,H,X,params)-alpha*norm(dLdW,'fro')^2*etaW
                etaW = beta*etaW;
            end
%             toc
            W(:, :, l) = W(:, :, l) - etaW*dLdW;
            etaW_tmp = [etaW_tmp etaW];
        end
    end
    fprintf('Reconstruction gradient over W: %f\n', mean(grad_W_recon));
    fprintf('Regularization WXH gradient over W: %f\n', mean(grad_W_WXH));
    fprintf('Regularization WW gradient over W: %f\n', mean(grad_W_WW));
    fprintf('Regularization L1W gradient over W: %f\n', mean(grad_W_L1W));
    fprintf('Line search step on W: %f\n', mean(etaW_tmp));
    grads_W_recon_all = [grads_W_recon_all mean(grad_W_recon)];
    grads_W_WXH_all = [grads_W_WXH_all mean(grad_W_WXH)];
    grads_W_WW_all = [grads_W_WW_all mean(grad_W_WW)];
    grads_W_L1W_all = [grads_W_L1W_all mean(grad_W_L1W)];
    etaW_all = [etaW_all mean(etaW_tmp)];
    toc
    
    % Calculate cost for this iteration
    Xhat = helper.reconstruct(W, H);    
    mask = find(params.M == 0); % find masked (held-out) indices 
    X(mask) = Xhat(mask); % replace data at masked elements with reconstruction, so masked datapoints do not effect fit
    errors(iter+1) = norm(X-Xhat, 'fro');
    
    % Shift to center factors
    if params.shift
        [W, H] = helper.shiftFactors(W, H);  
%         W = W+smallnum; % add small number to shifted W's, since multiplicative update cannot effect 0's
    end
    
    % Renormalize so rows of H have constant energy
    norms = sqrt(sum(H.^2, 2))';
    H = diag(1 ./ (norms+eps)) * H;
    for l = 1 : L
        W(:, :, l) = W(:, :, l) * diag(norms);
    end 
    
    % Stopping criteria... Stop if reach maxiter or if change in cost function is less than the tolerance
    if (iter == params.maxiter) || (norm(W(:)-W0(:)) < params.tolerance)
        errors = errors(1 : iter+1);  % trim vector
        lasttime = 1; 
        if iter>1
            params.lambda = 0; % Do one final CNMF iteration (no regularization, just prioritize reconstruction)
        end
    end
    
    % Plot to show progress
    if params.showPlot 
        SimpleWHPlot(W, H, Xhat,0); 
        title(sprintf('iteration #%i',iter));
        drawnow
    end
    
    W0 = W;
    
    if lasttime
        break
    end
end
grads = struct('grads_H_recon_all', grads_H_recon_all, ...
               'grads_H_WXH_all', grads_H_WXH_all, ...
               'grads_H_HH_all', grads_H_HH_all, ...
               'grads_H_L1H_all', grads_H_L1H_all, ...
               'etaH_all', etaH_all, ...
               'grads_W_recon_all', grads_W_recon_all, ...
               'grads_W_WXH_all', grads_W_WXH_all, ...
               'grads_W_WW_all', grads_W_WW_all, ...
               'grads_W_L1W_all', grads_W_L1W_all, ...
               'etaW_all', etaW_all);
   
% Undo zeropadding by truncating X, Xhat and H
X = X(:,L+1:end-L);
Xhat = Xhat(:,L+1:end-L);
H = H(:,L+1:end-L);

% Compute explained power of whole reconstruction and each factor
power = (sum(X(:).^2)-sum((X(:)-Xhat(:)).^2))/sum(X(:).^2);  % fraction power explained by whole reconstruction
[loadings,ind] = sort(helper.computeLoadingPercentPower(X,W,H),'descend'); % fraction power explained by each factor

% sort factors by loading power
if params.SortFactors
    W = W(:,ind,:);
    H = H(ind,:);
end

    function [X,N,T,K,L,params] = parse_seqNMF_params(X, inputs)
        % parse inputs, set unspecified parameters to the defaults
        
        % Get data dimensions
        [N, T] = size(X);

        p = inputParser; % 
        %USAGE: addOptional(p,'parametername',defaultvalue);
        addOptional(p,'K',10);
        addOptional(p,'L',100);
        addOptional(p,'lambda',.001);
        addOptional(p,'showPlot',1);
        addOptional(p,'maxiter',100);
        addOptional(p,'tolerance',-Inf);
        addOptional(p,'shift',1);
        addOptional(p,'lambdaL1W',0);
        addOptional(p,'lambdaL1H',0);
        addOptional(p,'W_fixed',0);
        addOptional(p,'W_init', nan); % depends on K--initialize post parse
        addOptional(p,'H_init', nan); % depends on K--initialize post parse
        addOptional(p,'SortFactors', 1); % sort factors by loading?
        addOptional(p,'lambdaOrthoW',0); % for this regularization: ||Wflat^TWflat||_1,i~=j
        addOptional(p,'lambdaOrthoH',0); % for this regularization: ||HSH^T||_1,i~=j
        addOptional(p,'useWupdate',1); % W update for cross orthogonality often doesn't change results much, and can be slow, so option to skip it 
        addOptional(p,'M',nan); % Masking matrix: default is ones; set elements to zero to hold out as masked test set
        addOptional(p, 'neg_prop', 0.2); % proportion of negative indices
        parse(p,inputs{:});
        L = p.Results.L; 
        K = p.Results.K; 
        params = p.Results; 
        
        % zeropad data by L
        X = [zeros(N,L),X,zeros(N,L)];
        [N, T] = size(X);

        % initialize W_init and H_init, if not provided
        if isnan(params.W_init)
            neg_indices = (rand(N, K, L) < params.neg_prop);
            params.W_init = max(X(:))*rand(N, K, L);
            params.W_init(neg_indices) = -params.W_init(neg_indices);
        end
        if isnan(params.H_init)
            params.H_init = max(X(:))*rand(K,T)./(sqrt(T/3)); % normalize so frobenius norm of each row ~ 1
        else
            params.H_init = [zeros(K,L),params.H_init,zeros(K,L)];
        end
        if isnan(params.M)
            params.M = ones(N,T);
        else
            params.M = [ones(N,L),params.M,ones(N,L)];
        end
    end
end