function [coeffs_W, coeffs_H, ids] = similarity_WH(W, H, W_hat, H_hat)
% Measure the similarity between the estimated factors W_hat and real
% factors W by matching each factor

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

%%
for ii = 1:K
    for jj = 1:Khat
        wk = squeeze(W(:,ii,:));
        wk_hat = squeeze(W_hat(:,jj,:));
        for l=-Lhat:Lhat
            wpad = cat(2,zeros(N,Lhat),wk,zeros(N,Lhat));
            wtmp = circshift(wpad, [0,l]);
            wtmp = wtmp(:,(Lhat+1):(end-L));
            Stmp(l+Lhat+1) = (wtmp(:)'*wk_hat(:))/((sqrt(wtmp(:)'*wtmp(:))*sqrt(wk_hat(:)'*wk_hat(:)))+eps);
        end
        [S(ii,jj), idx] = max(Stmp);
        shift(ii,jj) = idx-Lhat-1;
    end
end
% S(isnan(S)) = 0;
% S(S<0) = eps;
%%
temp = S;
coeffs_W = zeros(1,Khat);
coeffs_H = zeros(1,Khat);
ids = zeros(1, Khat);

% Matching each non-zero factor to all the factors of another reconstruction
for ii = 1:min(K,Khat)
    if ~any(temp(:))
        break
    end
    [r,c]= find(temp == max(temp(:)));
    coeffs_W(c(1)) = temp(r(1), c(1));
    Hk = H(r(1),:);
    Hk_hat = H_hat(c(1),:);
    Hpad = cat(2, zeros(1,Lhat),Hk,zeros(1,Lhat));
    % Shift H to opposite direction
    Htmp = circshift(Hpad,-shift(r(1), c(1)));
    Htmp = Htmp((Lhat+1):(end-Lhat));
    coeffs_H(c(1)) = (Htmp(:)'*Hk_hat(:))/((sqrt(Htmp(:)'*Htmp(:))*sqrt(Hk_hat(:)'*Hk_hat(:)))+eps);
    ids(c(1)) = r(1);

    temp(r(1),:) = 0;
    temp(:,c(1)) = 0;

end

coeffs_H(~coeffs_W) = [];
ids(~coeffs_W) = [];
coeffs_W(~coeffs_W) = [];
end