function [emds_W, emds_H, ids] = similarity_WH_EMD(W, H, W_hat, H_hat)
% Measure the similarity between the estimated factors W_hat and real
% factors W by matching each factor
% Use EMD metrics

szW = size(W);
if length(szW) == 2
    K = 1;
    N = szW(1);
    L = szW(2);
    W = reshape(W, [N,1,L]);
elseif length(szW) == 3
    N = szW(1);
    K = szW(2);
    L = szW(3);
end

szWhat = size(W_hat);
assert(N==szWhat(1), 'W and W_hat should have the same N!')
if length(szWhat) == 2
    Khat = 1;
    W_hat = reshape(W_hat, [N,1,szWhat(2)]);
elseif length(szW) == 3
    Khat = szWhat(2);
end
Lhat = size(W_hat, 3);

%% Normalize so that each row in H and H_hat has unit sum
norms = sum(H, 2)';
H = diag(1 ./ (norms+eps)) * H;
for l = 1 : L
    W(:, :, l) = W(:, :, l) * diag(norms);
end 

norms_hat = sum(H_hat, 2)';
H_hat = diag(1 ./ (norms_hat+eps)) * H_hat;
for l = 1 : Lhat
    W_hat(:, :, l) = W_hat(:, :, l) * diag(norms_hat);
end 

%% EMD options
opts = tfocs_SCD;
opts.continuation = 1;
opts.tol = 1e-4;
opts.stopCrit = 4;
opts.maxIts = 500;
opts.printEvery = 0;
opts.alg = 'N83';
continue_opts = continuation();
continue_opts.verbose = 0;

norms_W = sum(W,[1,3]);
norms_What = sum(W_hat,[1,3]);
S = nan(K,Khat);
shift = nan(K,Khat);

for ii = 1:K
    for jj = 1:Khat
%         fprintf('K=%d, K_hat=%d\n', ii, jj);
%         tic
        wk = squeeze(W(:,ii,:));
        wk_hat = squeeze(W_hat(:,jj,:));
        wpad = cat(2,zeros(N,Lhat),wk,zeros(N,Lhat));
        Stmp = nan(1,2*Lhat+1);
        
        % Ignore zero factors
        if (norms_W(ii)<1e-3*max(norms_W)) || (norms_What(jj)<1e-3*max(norms_What))
            continue
        end
        for l=-Lhat:Lhat            
            wtmp = circshift(wpad, [0,l]);
            wtmp = wtmp(:,(Lhat+1):(end-L));
            Stmp(l+Lhat+1) = compute_EMD(wtmp, wk_hat, opts, 'continuationOptions', continue_opts);
        end
        [S(ii,jj), idx] = min(Stmp);
        shift(ii,jj) = idx-Lhat-1;
%         toc
    end
end
% S(isnan(S)) = 0;
% S(S<0) = eps;
%%
temp = S;
emds_W = nan(1,Khat);
emds_H = nan(1,Khat);
ids = zeros(1, Khat);

% Matching each non-zero factor to all the factors of another reconstruction
for ii = 1:min(K,Khat)
    if ~any(temp(:))
        break
    end
    [r,c]= find(temp == min(temp(:)));
    emds_W(c(1)) = temp(r(1), c(1));
    Hk = H(r(1),:);
    Hk_hat = H_hat(c(1),:);
    Hpad = cat(2, zeros(1,Lhat),Hk,zeros(1,Lhat));
    % Shift H to opposite direction
    Htmp = circshift(Hpad,-shift(r(1), c(1)));
    Htmp = Htmp((Lhat+1):(end-Lhat));
    emds_H(c(1)) = compute_EMD(Htmp, Hk_hat, opts, 'continuationOptions', continue_opts);
    ids(c(1)) = r(1);

    temp(r(1),:) = nan;
    temp(:,c(1)) = nan;

end

emds_H(isnan(emds_W)) = [];
ids(isnan(emds_W)) = [];
emds_W(isnan(emds_W)) = [];
end