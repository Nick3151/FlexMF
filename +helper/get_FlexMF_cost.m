function [recon_err, reg_cross, reg_W, reg_H] = get_FlexMF_cost(X,W,H)

[N,K,L] = size(W);
[~,T] = size(H);
Xhat = helper.reconstruct(W, H);
smoothkernel = ones(1,(2*L)-1);  
WTX = zeros(K, T);

for l = 1 : L
    X_shifted = [X(:, l:T) zeros(N, l-1)];
    WTX = WTX + W(:, :, l)' * X_shifted;
end   

% Compute regularization terms for H update
WTXSHT = (conv2(abs(WTX), smoothkernel, 'same')*H'); 
WTXSHT = WTXSHT.*~eye(K);
recon_err = sum((X(:)-Xhat(:)).^2)/2;
%sqrt(sum((X(:)-Xhat(:)).^2)) + lambda.*norm(WTXSHT,1);
reg_cross = norm(WTXSHT(:),1);
reg_W = norm(W(:),1);
reg_H = norm(H(:),1);

end