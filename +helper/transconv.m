function WTX = transconv(W,X)
    % Transpose convolution
    
    % get size of W and H
    [N,K,L] = size(W);
    [~,T] = size(X);
    
    % zeropad by L
    X = [zeros(N,L),X,zeros(N,L)];
    T = T+2*L;
    WTX = zeros(K, T);
    
    for l = 1 : L
        X_shifted = circshift(X,[0,-l+1]); 
        WTX = WTX + W(:, :, l)' * X_shifted;
    end 
    
    % undo zer0padding
    WTX = WTX(:,(L+1):(end-L)); 
end