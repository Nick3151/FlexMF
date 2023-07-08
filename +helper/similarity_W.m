function [coeffs, ids] = similarity_W(W, W_hat)
% Measure the similarity between the estimated factors W_hat and real
% factors W by matching each factor

% Shift factors to center
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
if length(szWhat) == 2
    Khat = 1;
    assert(N==szWhat(1) && L==szWhat(2), 'Dimensions of W and W_hat do not match!')
    W_hat = reshape(W_hat, [N,1,L]);
elseif length(szW) == 3
    Khat = szWhat(2);
    assert(N==szWhat(1) && L==szWhat(3), 'Dimensions of W and W_hat do not match!')
end

%%
for ii = 1:K
    for jj = 1:Khat
        wk = squeeze(W(:,ii,:));
        wk_hat = squeeze(W_hat(:,jj,:));
        for l=-L:L
            wpad = cat(2,zeros(N,L),wk,zeros(N,L));
            wtmp = circshift(wpad, [0,l]);
            wtmp = wtmp(:,(L+1):(end-L));
            Stmp(l+L+1) = (wtmp(:)'*wk_hat(:))/((sqrt(wtmp(:)'*wtmp(:))*sqrt(wk_hat(:)'*wk_hat(:)))+eps);
        end
        S(ii,jj) = max(Stmp);
    end
end
% S(isnan(S)) = 0;
% S(S<0) = eps;
%%
temp = S;
coeffs = zeros(1,Khat);
ids = zeros(1, Khat);

% Matching each non-zero factor to all the factors of another reconstruction
for ii = 1:min(K,Khat)
    if ~any(temp(:))
        break
    end
    [r,c]= find(temp == max(temp(:)));
    coeffs(c(1)) = temp(r(1), c(1));
    ids(c(1)) = r(1);

    temp(r(1),:) = 0;
    temp(:,c(1)) = 0;

end

coeffs(~coeffs) = [];
ids(~ids) = [];
end