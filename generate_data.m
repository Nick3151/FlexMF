function [data,W,H,V_hat] = generate_data(T,Nneurons,Dt,NeuronNoise,SeqNoiseTime,SeqNoiseNeuron,stretch,bin,neg,seed)
%rng(2001)
if seed == 0
    rng shuffle
else
    rng(seed)
end

additional_neurons = 0;

%% Parameters
% V = data.sequences(10000,randi(10,2,1)+5,randi(4,2,1),.01,rand(2),ones(2,1)*0.95);
% Nneurons = [5,15,10,8]; % the number of neurons in each sequence
% Dt = [2,1,3,3]; % the number of time steps between each neuron in the sequence
% Pseq = []; % the probability of the sequence occuring
% NeuronNoise = 0.01; % the noise in a neurons firing rate
% SeqNoiseTime = [0.2,0.2,0.1,0.1]; % the noise in the sequence aka jitter (p of each neuron jittered 1 dt)
% SeqNoiseNeuron = [0.95,0.95,0.95,0.95]; % the probability that a neuron participates in a given seq
% T = 1000;
% Share = []; % the propotion of the chain that is shared in other sequences
%% Calculate useful things
N = sum(Nneurons)+additional_neurons; % Total number of neurons
K = length(Nneurons); % The number of sequences
lseq = Dt.*Nneurons; % the length of each sequences

j = 1;
neurons = {};
for ii = 1:length(Nneurons)
    neurons{ii} = j:j+Nneurons(ii)-1;
    j = j+Nneurons(ii);
end

%% MAKE H's
xx = zeros(1,T);
H = zeros(K,T);
%randomly distribute seq starting points preventing the same sequence from initiation during itself
% for ii = 1:length(lseq)
%     pos(ii,:) = cumsum(randi(lseq(ii)+50,50,1)+lseq(ii)); %sp sets how often the seq happen
%     temp = pos(ii,:);
%     H(ii,temp(temp<T))= 1;
% %     H(ii,logical(xx)) = 0;
% %     for jj = 1:length(pos(ii,:))
% %         xx(pos(ii,jj):pos(ii,jj)+lseq(ii)) = 1;
% %     end
%       
% end
nn = K*1000; % make smaller
if stretch > 0
    stretches = randi([-stretch, stretch],nn,1);
else
    stretches = zeros(nn,1);
end

if seed == 0
    rng shuffle
else
    rng(seed)
end

% time courses with a minimum gap.
temp = [];
for j = 1:nn
    % need to fix this to work for stretches = 0
    temp = [temp;randi(100,1,1)+max(lseq)];
%     temp = [temp;randi(100,1,1)+(max(lseq)/max(Dt)*(max(Dt)+stretches(ii)))];
end
temp = cumsum(temp);

%indx = randi(nseq,nn,1);
indx = ones(1000,1)';
for ii = 2:K
    indx = [indx,ii*ones(nn/K,1)'];
end
indx = indx(randperm(nn));
for ii = 1:K
    H(ii,temp((indx == ii))) = 1;
    Hs{ii} = stretches((indx == ii));
end
H = H(:,1:T);

%% Make Data using noise parameters in the reconstruction
%leng = max(Dt)+ max(lseq) + (stretch)*max(lseq); 
leng = (max(lseq)/max(Dt)*(max(Dt)+stretch));
L = leng+150;
W = zeros(N,K,L);
% H(:,T-(2*(max(lseq)/max(Dt)*(max(Dt)+(stretch)))):T) = 0;
H(:,T-(2*(max(lseq))):T) = 0;
[~,T] = size(H);
V_hat = zeros(N,T);

%Dont forget these things!
%NeuronNoise = [0.1,0.05,0.2,0.13]; % the noise in a neurons firing rate
%SeqNoiseTime = [0.1,0.2,0.1,0.1]; % the noise in the sequence aka jitter (p of each neuron jittered 1 dt)
%SeqNoiseNeuron = [0.9,0.9,0.9,0.9]; % the probability that a neuron participates in a given seq
for ii = 1:K % go through each factor
    ind = find(H(ii,:));
    
    % neg: proportion of negative indices
    neg_indices = (rand(1,length(neurons{ii})) < neg);
    
    Dt_temp = Dt(ii);
    Wi = zeros(N,L);          
    temp = ones(1,length(neurons{ii}));
    temp(neg_indices) = -temp(neg_indices);
    temp2 = zeros(length(neurons{ii}),Dt_temp*Nneurons(ii));
    temp2(:,1:Dt_temp:Dt_temp*Nneurons(ii)) = diag(temp);    
    Wi(neurons{ii},50:49+size(temp2,2)) = temp2;  
    W(:,ii,:) = Wi;
    
    for jj = 1:sum(H(ii,:)) % go through each iteration of the sequence
        tempH = zeros(1,size(H,2));
        tempH(ind(jj)) = 1;
                   
        if stretch > 0 % change the dt for each instance
            Dt_temp = Dt(ii)+Hs{ii}(jj);%+(randi(stretch))
            %*(-1+(2*(rand(1)>0.5))); If you want compression as well   
            tempW = zeros(N,L);          
            temp2 = zeros(length(neurons{ii}),Dt_temp*Nneurons(ii));
            temp2(:,1:Dt_temp:Dt_temp*Nneurons(ii)) = diag(temp);    
            tempW(neurons{ii},50:49+size(temp2,2)) = temp2;  
        else 
            tempW = Wi;      
        end

        % neurons are jittered with some lambda
        %shifts = poissrnd(SeqNoiseTime(ii),N,1).*(1-2*(rand(N,1)>0.5));
        shifts = round(normrnd(0,SeqNoiseTime(ii),N,1));
        %shifts(abs(shifts) >5) = 0; %stop poiss from getting to big
        for idx = 1:N
            tempW(idx,:) = circshift(tempW(idx,:),shifts(idx)*1);
        end
        
        % neurons participate with some p
        tempW(rand(N,1)>SeqNoiseNeuron(ii),:) = 0;
    
        [tempW, tempH] = helper.shiftFactors(tempW, tempH);
        % Compensate for factor shift due to stretching
        shift = circshift(1:length(tempH),floor((Dt(ii)*Nneurons(ii)- Dt_temp*Nneurons(ii))/2));
        tempH = tempH(:,shift);
        newData = conv2(tempH,tempW);
        %V_hat = V_hat + newData(:,ceil(size(tempW,2)/2):end - floor(size(tempW,2)/2));
        V_hat = V_hat + newData(:,1:(end - size(tempW,2)+1));

    end
end


%% could add indepent noise later
V_hat = V_hat + (rand(size(V_hat))<NeuronNoise);
%V_hat  = V_hat./V_hat;
V_hat(isnan(V_hat)) = 0;
%%
%data = V_hat;
if ~bin
    filtbio = [zeros(1,10*10) exp(-(1:10*10)/10)]; 
    data = conv2(V_hat,filtbio,'same');
else
    data = V_hat;
end

[W, H] = helper.shiftFactors(W, H);

if ~bin
    for ii = 1:size(W,2)
        W(:,ii,:) = conv2(squeeze(W(:,ii,:)),filtbio,'same');
    end
end
rng shuffle
end
