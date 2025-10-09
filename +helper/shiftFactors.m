function [W,H] = shiftFactors(W,H)
    % shift factors by center of mass
    
    % get size of W and H
    [K,T] = size(H);
    szW = size(W);
    if length(szW) == 2
        assert(K==1, 'dimensions of W and H do not match!')
        N = szW(1);
        L = szW(2);
        W = reshape(W, [N,1,L]);
    elseif length(szW) == 3
        assert(K==szW(2), 'dimensions of W and H do not match!')
        N = szW(1);
        L = szW(3);
    end
    
    if L>1 % if L=1, no room to shift
    
        center = max(floor(L/2),1); % middle bin

        % pad with zeros, for use with circshift. data is zeropadded within seqNMF, so don't need to pad H
        Wpad = cat(3,zeros(N,K,L),W,zeros(N,K,L));

        for k = 1:K
            % compute center of mass
            temp = sum(abs(squeeze(W(:,k,:))),1);
            cmass = max(floor(sum(temp.*(1:length(temp)))/sum(temp)),1);          
            %Wpad(:,k,:) = circshift(squeeze(Wpad(:,k,:)),center-cmass,2);
            %Changing for compatibility with 2016, this is the same
            Wpad(:,k,:) = circshift(squeeze(Wpad(:,k,:)),[0,center-cmass]); 
            %H(k,:) = circshift(H(k,:),cmass-center); 
            H(k,:) = circshift(H(k,:),[0,cmass-center]); 
        end


        % undo zero pad
        W = Wpad(:,:,(L+1):(end-L));
        
        if K==1 && length(szW) == 2
            W = reshape(W,[N,L]);
        end

    end