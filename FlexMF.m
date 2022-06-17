function [W, H, errors, loadings, power] = FlexMF(X, varargin)
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
% 'useWupdate'      1                                   Wupdate for cross orthogonality often doesn't change results much, and can be slow, so option to remove  
% 'M'               ones(N,T)                           Masking matrix if excluding a random test set from the fit
% 'neg_prop'        0.2                                 Proportion of negative indices
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
W_pre = params.W_init;
H_pre = params.H_init;
W = params.W_init;
lasttime = 0;
errors = zeros(params.maxiter, 1);

for iter = 1 : params.maxiter
    tic
    fprintf('Iter %d\n', iter);
    
    fprintf('Updating W\n');
%     W0 = max(W_pre(:))*rand(N, K, L);
    W0 = W_pre;
%     H_tmp = H_pre + 0.05*max(H_pre(:))*rand(K,T); 
    W = updateW(W0, H_pre, X, params);   
    
    fprintf('Updating H\n');
%     H0 = max(H_pre(:))*rand(K,T); 
    H0 = H_pre;
    W_tmp = W + 0.05*max(W_pre(:))*rand(N, K, L);
    H = updateH(W, H0, X, params);   

    toc
    
    % Calculate cost for this iteration
    Xhat = helper.reconstruct(W, H);    
    mask = find(params.M == 0); % find masked (held-out) indices 
    X(mask) = Xhat(mask); % replace data at masked elements with reconstruction, so masked datapoints do not effect fit
    errors(iter) = norm(X-Xhat, 'fro');
    
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
    
    % Stopping criteria... Stop if reach maxiter or if change in cost function is less than the tolerance
    if (iter == params.maxiter) || (norm(W(:)-W_pre(:)) < params.tolerance)
        errors = errors(1 : iter);  % trim vector
        lasttime = 1; 
%         if iter>1
%             params.lambda = 0; % Do one final CNMF iteration (no regularization, just prioritize reconstruction)
%         end
    end
    
    % Plot to show progress
    if params.showPlot 
        Xhat = helper.reconstruct(W, H);
        SimpleWHPlot(W, H, Xhat); 
        title(sprintf('iteration #%i',iter));
        drawnow
    end
    
    W_pre = W;
    H_pre = H;
    
    if lasttime
        break
    end
end
   
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
        addOptional(p,'lambda',.01);
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