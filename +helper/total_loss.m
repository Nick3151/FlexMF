function loss = total_loss(W, H, X, params)
    [N,K,L] = size(W);
    % Total loss function
    if params.lambda>0
        smoothkernel = ones(1,(2*L)-1);  % for factor competition
        WTX = helper.transconv(W, X);
        RWXH = conv2(WTX, smoothkernel, 'same')*H';
        RWXH = params.lambda*(norm(RWXH(:),1)-norm(diag(RWXH),1));
    else
        RWXH = 0;
    end
    
    if params.lambdaOrthoH>0
        HS = conv2(H, smoothkernel, 'same');
        HSHT = HS*H';
        RHH = params.lambdaOrthoH/2*(norm(HSHT(:),1)-norm(diag(HSHT),1));
    else
        RHH = 0;
    end
    
    if params.lambdaOrthoW>0
        Wflat = sum(W,3);
        WTW = Wflat'*Wflat;
        RWW = params.lambdaOrthoW/2*(norm(WTW(:),1)-norm(diag(WTW),1));
    else
        RWW = 0;
    end
    
    if params.lambdaL1H>0
        RH = params.lambdaL1H*norm(H(:),1);
    else
        RH = 0;
    end
    
    if params.lambdaL1W>0
        RW = params.lambdaL1W*norm(W(:),1);
    else
        RW = 0;
    end
    loss = norm(helper.reconstruct(W, H)-X, 'fro')^2/2 + RWXH + RH + RW + RHH + RWW;
end