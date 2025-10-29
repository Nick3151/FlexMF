function [W, H, cost, errors, loadings, power, M, R] = FlexMF(X, varargin)
%
% USAGE: 
%
% [W, H, errors, loadings, power] = FlexMF(X, ...    % X is the data matrix
%       'K', 10, 'L', 20, 'lambda', .1, ...        % Other inputs optional
%       'W_init', W_init, 'H_init', H_init, ...
%       'showPlot', 1, 'maxiter', 20, 'tolerance', -Inf, 'shift', 1, ... 
%       'lambdaL1W', 0, 'lambdaL1H', 0, ...
%       'M', M)
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
% 'alpha_H'         1e-3                                Regularization parameter for H
% 'alpha_W'         1e-6                                Regularization parameter for W
% 'W_init'          max(X(:))*rand(N,K,L)               Initial W
% 'H_init'          max(X(:))*rand(K,T)./(sqrt(T/3))    Initial H (rows have norm ~1 if max(data) is 1)
% 'showPlot'        1                                   Plot every iteration? no=0
% 'maxiter'         100                                 Maximum # iterations to run
% 'tolerance'       1e-3                                Stop if improved less than this;  Set to -Inf to always run maxiter
% 'alg'             'N83'                               Algorithm
% 'shift'           1                                   Shift factors to center; Helps avoid local minima
% 'lambdaL1W'       0                                   L1 sparsity parameter; Increase to make W's more sparse
% 'lambdaL1H'       0                                   L1 sparsity parameter; Increase to make H's more sparse
% 'Reweight'        0                                   Whether to use reweighted L1 minimization
% 'W_fixed'         0                                   Fix W during the fitting proceedure   
% 'SortFactors'     1                                   Sort factors by loadings
% 'useWupdate'      1                                   Wupdate for cross orthogonality often doesn't change results much, and can be slow, so option to remove  
% 'M'               ones(N,T)                           Masking matrix if excluding a random test set from the fit
% 'neg_prop'        0.2                                 Proportion of negative indices
% 'EMD              0                                   Optimize EMD instead of reconstruction error
% 'lambda_R'        1                                   Penalty coefficient on residual term for unbalanced EMD
% 'lambda_M'        1e-4                                Penalty coefficient on motion field for unbalanced EMD
% 'lambda_TV'       0                                   TV norm of W parmater; Increase to make W more smooth along the time dimension
% 'verbal'          1                                   Print intermediate output?
% ------------------------------------------------------------------------
% OUTPUTS:
%
% W                         NxKxL tensor containing factor exemplars
% H                         KxT matrix containing factor timecourses
% cost                      (#Iterations+1)x1 vector containing 
%                               reconstruction error(or EMD) at each iteration. 
%                               cost(1) is error before 1st iteration.
% errors                    (#Iterations+1)x4 matrix containing
%                               reconstrcution and regularization  
%                               errors for each iteration
% loadings                  1xK vector containing loading of each factor 
%                               (Fraction power in data explained by each factor)
% power                     Fraction power in data explained 
%                               by whole reconstruction
%
%                           Note, if doing fit with masked (held-out) data,
%                               the cost and power do not include masked
%                               (M==0) test set elements
% M,R                       M and R for EMD optimization
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
W_pre = params.W_init;
H_pre = params.H_init;
W = params.W_init;

Xhat = helper.reconstruct(W_pre, H_pre); 
mask = find(params.M == 0); % find masked (held-out) indices 
X(mask) = Xhat(mask); % replace data at masked elements with reconstruction, so masked datapoints do not effect fit

lasttime = 0;
cost = zeros(params.maxiter+1, 1);
errors = zeros(params.maxiter+1, 4);
if params.EMD
    opts = tfocs_SCD;
    opts.continuation = 1;
    opts.tol = 1e-4;
    opts.stopCrit = 4;
    opts.maxIts = 500;
    opts.printEvery = 0;
    opts.alg = 'N83';
    continue_opts = continuation();
    continue_opts.verbose = 0;
    M_pre = zeros(N,T);
    R_pre = zeros(N,T);
    X_ = [X M_pre R_pre];
    Xhat_ = [Xhat M_pre R_pre];
%     cost(1) = compute_EMD(X,Xhat,opts, 'continuationOptions', continue_opts);
    cost(1) = norm(X_(:)-Xhat_(:));
    L1_Ms = zeros(params.maxiter+1, 2); % L1 norms of M after updating W/H
    L1_Rs = zeros(params.maxiter+1, 2); % L1 norms of R after updating W/H
    L1_Ws = zeros(params.maxiter+1, 2); % L1 norms of W after updating W/H
    L1_Hs = zeros(params.maxiter+1, 2); % L1 norms of H after updating W/H
else
%     cost(1) = sqrt(mean((X(:)-Xhat(:)).^2));
    cost(1) = norm(X(:)-Xhat(:));
end
[recon_err, reg_cross, reg_W, reg_H] = helper.get_FlexMF_cost(X,W_pre,H_pre);
errors(1,:) = [recon_err, reg_cross, reg_W, reg_H];

for iter = 1 : params.maxiter
    params.currentiter = iter;
    if params.verbal
        fprintf('Iter %d\n', iter);
    end
    if iter > 5
        dcost = cost(iter) - mean(cost((iter-5):(iter-1)));
        if params.EMD
            X_ = [X M_pre R_pre];
            dcost_norm = dcost/norm(X_(:));
        else
            dcost_norm = dcost/norm(X(:));
        end
        if params.verbal
            fprintf('dcost=%f\n', dcost);
            fprintf('dcost/X=%f\n', dcost_norm);
        end
    end
    % Stopping criteria... Stop if reach maxiter or if change in cost function is less than the tolerance
    % if (iter == params.maxiter) || ((iter>5) && (dW < params.tolerance))
    if (iter == params.maxiter) || ((iter>5) && (abs(dcost_norm) < params.tolerance))
        cost = cost(1 : iter+1);  % trim vector
        errors = errors(1 : iter+1, :);
        if params.EMD
            L1_Ms = L1_Ms(1: iter+1, :);
            L1_Rs = L1_Rs(1: iter+1, :);
            L1_Hs = L1_Hs(1: iter+1, :);
            L1_Ws = L1_Ws(1: iter+1, :);
        end
        lasttime = 1; 
%         if iter>1
%             params.lambda = 0; % Do one final CNMF iteration (no regularization, just prioritize reconstruction)
%         end
    end
    
    if ~params.W_fixed
    %     tic
        if params.verbal
            fprintf('Updating W\n');
        end
        W0 = W_pre;
        if params.EMD
            M0 = M_pre;
            R0 = R_pre;
%             M0 = zeros(N,T);
%             R0 = zeros(N,T);
            [W, M, R, out] = updateW_EMD(W0, H_pre, X, M0, R0, params);
            L1_Ms(iter+1,1) = norm(M(:),1)/norm(X(:),1);
            L1_Rs(iter+1,1) = norm(R(:),1)/norm(X(:),1);
            L1_Ws(iter+1,1) = norm(W(:),1)/norm(X(:),1);
            L1_Hs(iter+1,1) = norm(H_pre(:),1)/norm(X(:),1);
        else
            W = updateW(W0, H_pre, X, params); 
        end
    else
        W = W_pre;
        if params.EMD
            M = M_pre;
            R = R_pre;
        end
    end
    
    if params.verbal
        fprintf('Updating H\n');
    end

    H0 = H_pre;
    if params.EMD
        M0 = M;
        R0 = R;
%     M0 = zeros(N,T);
%     R0 = zeros(N,T);
        [H, M, R, out] = updateH_EMD(W, H0, X, M0, R0, params);
        L1_Ms(iter+1,2) = norm(M(:),1)/norm(X(:),1);
        L1_Rs(iter+1,2) = norm(R(:),1)/norm(X(:),1);
        L1_Ws(iter+1,2) = norm(W(:),1)/norm(X(:),1);
        L1_Hs(iter+1,2) = norm(H(:),1)/norm(X(:),1);
    else
        H = updateH(W, H0, X, params);   
    end

%     toc
    if ~params.W_fixed
        % Shift to center factors
        if params.shift
            [W, H] = helper.shiftFactors(W, H);  
        end
        
        % Renormalize so rows of H have constant energy
        norms = sqrt(sum(H.^2, 2))';
        H = diag(1 ./ (norms+eps)) * H;
        for l = 1 : L
            W(:, :, l) = W(:, :, l) * diag(norms);
        end 
    end
    
    % Calculate cost for this iteration
    Xhat = helper.reconstruct(W, H);    
    mask = find(params.M == 0); % find masked (held-out) indices 
    X(mask) = Xhat(mask); % replace data at masked elements with reconstruction, so masked datapoints do not effect fit
    if params.EMD
%         cost(iter+1) = compute_EMD(X,Xhat,opts, 'continuationOptions', continue_opts);
        X_ = [X M_pre R_pre];
        Xhat_ = [Xhat M R];
        cost(iter+1) = norm(X_(:)-Xhat_(:));
    else
%         cost(iter+1) = sqrt(mean((X(:)-Xhat(:)).^2));
        cost(iter+1) = norm(X(:)-Xhat(:));
    end
    [recon_err, reg_cross, reg_W, reg_H] = helper.get_FlexMF_cost(X,W,H);
    errors(iter+1,:) = [recon_err, reg_cross, reg_W, reg_H];
%     dW = sqrt(mean((W(:)-W_pre(:)).^2));
%     fprintf('dW=%f\n', dW);
    
    % Plot to show progress
    if params.showPlot 
        Xhat = helper.reconstruct(W, H);
        SimpleWHPlot_patch(W, H, 'Data', Xhat); 
        title(sprintf('iteration #%i',iter));
        drawnow
    end
    
    W_pre = W;
    H_pre = H;
    if params.EMD
        M_pre = M;
        R_pre = R;
    end
    if lasttime
        break
    end
end
   
% Undo zeropadding by truncating X, Xhat and H
X = X(:,L+1:end-L);
Xhat = Xhat(:,L+1:end-L);
H = H(:,L+1:end-L);
if params.EMD
    M = M(:,L+1:end-L);
    R = R(:,L+1:end-L);
end

% Compute explained power of whole reconstruction and each factor
power = (sum(X(:).^2)-sum((X(:)-Xhat(:)).^2))/sum(X(:).^2);  % fraction power explained by whole reconstruction
[loadings,ind] = sort(helper.computeLoadingPercentPower(X,W,H),'descend'); % fraction power explained by each factor

% sort factors by loading power
if params.SortFactors
    W = W(:,ind,:);
    H = H(ind,:);
end

if params.verbal && params.EMD
    figure;
    plot(1:.5:iter+1.5, reshape(L1_Ms', 1, []), 1:.5:iter+1.5, reshape(L1_Rs', 1, []))
    hold on
    plot(1:.5:iter+1.5, reshape(L1_Ws', 1, []), 1:.5:iter+1.5, reshape(L1_Hs', 1, []))
    legend('||M||_1', '||R||_1', '||W||_1', '||H||_1')
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
        addOptional(p,'alpha_H',1e-3);
        addOptional(p,'alpha_W',1e-6);
        addOptional(p,'showPlot',1);
        addOptional(p,'verbal',1);
        addOptional(p,'maxiter',100);
        addOptional(p,'tolerance',1e-3);
        addOptional(p,'alg','N83');
        addOptional(p,'shift',1);
        addOptional(p,'lambdaL1W',0);
        addOptional(p,'lambdaL1H',0);
        addOptional(p,'Reweight',0);
        addOptional(p,'W_fixed',0);
        addOptional(p,'W_init', nan); % depends on K--initialize post parse
        addOptional(p,'H_init', nan); % depends on K--initialize post parse
        addOptional(p,'SortFactors', 1); % sort factors by loading?
        addOptional(p,'useWupdate',1); % W update for cross orthogonality often doesn't change results much, and can be slow, so option to skip it 
        addOptional(p,'M',nan); % Masking matrix: default is ones; set elements to zero to hold out as masked test set
        addOptional(p, 'neg_prop', 0.2); % proportion of negative indices
        addOptional(p, 'EMD', 0);  % Optimize EMD instead of reconstruction error
        addOptional(p, 'lambda_R', 1); % Penalty coefficient on residual term for unbalanced EMD
        addOptional(p, 'lambda_M', 1e-4); % Penalty coefficient on motion field for unbalanced EMD
        addOptional(p, 'lambda_TV', 0); % TV norm of W along the time dimension
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