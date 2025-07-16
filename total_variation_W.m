function y = total_variation_W(N, K, L, W, mode)
% Total variation on W long time dimension
% V_k = W_k*D'
% V*_k = Y_k*D

% divergence matrix
D = eye(L) - diag(ones(L-1,1),-1);

switch mode
    case 0
        y = {[N,K,L], [N,K,L]};
    case 1
        assert(isequal(size(W), [N,K,L]), 'Incorrect dimension of W')
        y = zeros(N,K,L);
        for k=1:K
            Wk = squeeze(W(:,k,:));
            y(:,k,:) = Wk*D';
        end
    case 2
        y = zeros(N,K,L);
        for k=1:K
            Wk = squeeze(W(:,k,:));
            y(:,k,:) = Wk*D;
        end
end
