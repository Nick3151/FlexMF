function [data,W,H,X_hat] = generate_data(T,Nneurons,Dt,varargin)
% Parameters
% Nneurons = [5,15,10,8]; % the number of neurons in each sequence
% Dt = [2,1,3,3]; % the number of time steps between each neuron in the sequence
% Pseq = []; % the probability of the sequence occuring
% Magnitudes = ones(number_of_seqences, 1); % the activation magnitudes of each motif
% noise = 0.001; % probability of added noise in each bin
% jitter = zeros(number_of_seqences,1); % Jitter time std
% participation = 1.*ones(number_of_seqences,1); % Participation probability = 100%
% gap = 100;  % the maximum gap between sequences
% warp = 0; % the maximum warping time
% overlap_t = 0; if temporal overlap between different motifs is allowed
% overlap_n = 0; The proportion of shared neurons between different motifs
% len_burst = 10; The time of continuous firing
% dynamic = 1; if transient dynamics is being simulated
% neg = 0; % Proportion of negative indices in W
% seed = 0; % Random seed
% additional_neurons = 0; % Number of additional neurons

%% Parse inputs
K = length(Nneurons); % The number of motifs
validArrayK = @(x) isnumeric(x) && length(x)==K;

p  = inputParser;
addOptional(p, 'Magnitudes', ones(K, 1), validArrayK)
addOptional(p, 'noise', .001, @isscalar)
addOptional(p, 'jitter', zeros(K, 1), validArrayK)
addOptional(p, 'participation', ones(K, 1), validArrayK)
addOptional(p, 'gap', 100, @isscalar)
addOptional(p, 'warp', 0, @isscalar)
addOptional(p, 'overlap_t', 0, @isscalar)
addOptional(p, 'overlap_n', 0, @isscalar)
addOptional(p, 'len_burst', 10, @isscalar)
addOptional(p, 'dynamic', 1, @isscalar)
addOptional(p, 'neg', 0, @isscalar)
addOptional(p, 'seed', 0, @isscalar)
addOptional(p, 'additional_neurons', 0, @isscalar)
parse(p, varargin{:})

Magnitudes = p.Results.Magnitudes;
noise = p.Results.noise;
jitter = p.Results.jitter;
participation = p.Results.participation;
warp = p.Results.warp;
gap = p.Results.gap;
overlap_t = p.Results.overlap_t;
overlap_n = p.Results.overlap_n;
len_burst = p.Results.len_burst;
dynamic = p.Results.dynamic;
neg = p.Results.neg;
seed = p.Results.seed;
additional_neurons = p.Results.additional_neurons;

if seed == 0
    rng shuffle
else
    rng(seed)
end

%% Calculate useful things
lseq = Dt.*Nneurons; % the length of each sequences

j = 1;
neurons = {};
assert((overlap_n<1) && (overlap_n>=0), 'Overlap_n should be between [0,1)!')
num_overlap = floor(overlap_n*Nneurons);
for k = 1:length(Nneurons)
    neurons{k} = j:j+Nneurons(k)-1;
    j = j+Nneurons(k)-num_overlap(k);
end

N = j+num_overlap(K)-1+additional_neurons; % Total number of neurons

%% MAKE H's
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
if warp > 0
    stretches = randi([-warp, warp],nn,1);
else
    stretches = zeros(nn,1);
end

if seed == 0
    rng shuffle
else
    rng(seed+1)
end

% time courses with a minimum gap.
temp = [];
for j = 1:nn
    if overlap_t
        temp = [temp;randi(gap)];
    else
        temp = [temp;randi(gap)+max(lseq)];
    end
end
temp = cumsum(temp);

%indx = randi(nseq,nn,1);
indx = ones(1000,1)';
for k = 2:K
    indx = [indx,k*ones(1000,1)'];
end
indx = indx(randperm(nn));
for k = 1:K
    H(k,temp((indx == k))) = Magnitudes(k);
    Hs{k} = stretches((indx == k));
end
H = H(:,1:T);

%% Make Data using noise parameters in the reconstruction
%leng = max(Dt)+ max(lseq) + (stretch)*max(lseq); 
leng = (max(lseq)/max(Dt)*(max(Dt)+warp));
L = max(lseq)+20;
W = zeros(N,K,L);
% H(:,T-(2*(max(lseq)/max(Dt)*(max(Dt)+(stretch)))):T) = 0;
H(:,(T-leng):T) = 0;
[~,T] = size(H);
X_hat = zeros(N,T);

%Dont forget these things!
%NeuronNoise = [0.1,0.05,0.2,0.13]; % the noise in a neurons firing rate
%SeqNoiseTime = [0.1,0.2,0.1,0.1]; % the noise in the sequence aka jitter (p of each neuron jittered 1 dt)
%SeqNoiseNeuron = [0.9,0.9,0.9,0.9]; % the probability that a neuron participates in a given seq
for k = 1:K % go through each factor
    ind = find(H(k,:));
    
    % neg: proportion of negative indices
    neg_indices = (rand(1,length(neurons{k})) < neg);
    
    Dt_temp = Dt(k);
    Wi = zeros(N,L);          
    temp = ones(1,length(neurons{k}));
    temp(neg_indices) = -temp(neg_indices);
    l = Dt_temp*Nneurons(k);
    temp2 = zeros(length(neurons{k}),l);
    temp2(:,1:Dt_temp:Dt_temp*Nneurons(k)) = diag(temp);    
    Wi(neurons{k},(1:l)+2*max(jitter)) = temp2;  
    W(:,k,:) = Wi;
    
    for jj = 1:sum(H(k,:)) % go through each iteration of the sequence
        tempH = zeros(1,size(H,2));
        tempH(ind(jj)) = 1;
                   
        if warp > 0 % change the dt for each instance
            Dt_temp = Dt(k)+Hs{k}(jj);%+(randi(stretch))
            %*(-1+(2*(rand(1)>0.5))); If you want compression as well  
            l = Dt_temp*Nneurons(k);  
            tempW = zeros(N,leng);      
            
            temp2 = zeros(length(neurons{k}),l);
            temp2(:,1:Dt_temp:l) = diag(temp);    
            tempW(neurons{k},(1:l)+2*max(jitter)) = temp2;  
        else 
            tempW = Wi;      
        end

        % neurons are jittered with some lambda
        %shifts = poissrnd(SeqNoiseTime(ii),N,1).*(1-2*(rand(N,1)>0.5));
        shifts = round(normrnd(0,jitter(k),N,1));
        %shifts(abs(shifts) >5) = 0; %stop poiss from getting to big
        for idx = 1:N
            tempW(idx,:) = circshift(tempW(idx,:),shifts(idx)*1);
        end
        
        % neurons participate with some p
        tempW(rand(N,1)>participation(k),:) = 0;
    
        [tempW, tempH] = helper.shiftFactors(tempW, tempH);
        % Compensate for factor shift due to stretching
        shift = circshift(1:length(tempH),floor((Dt(k)*Nneurons(k)- Dt_temp*Nneurons(k))/2));
        tempH = tempH(:,shift);
        newData = conv2(tempH,tempW);
        %V_hat = V_hat + newData(:,ceil(size(tempW,2)/2):end - floor(size(tempW,2)/2));
        X_hat = X_hat + newData(:,1:(end - size(tempW,2)+1));

    end
end

%% Add indepent noise
X_noise = (rand(size(X_hat))<noise);

%% Continuous firing
if isempty(len_burst)
    len_burst = 1;
end

X_hat_tmp = conv2(X_hat, ones(1,len_burst));
X_hat = X_hat_tmp(:,1:T);

for k = 1:K
    W_tmp = conv2(squeeze(W(:,k,:)), ones(1,len_burst));
    W(:,k,:) = W_tmp(:,1:L);
end

X = X_hat + (~X_hat).*X_noise;
X(isnan(X)) = 0;

%%
%data = V_hat;
if dynamic
    filtbio = exp(-(0:8)/4);
    data_tmp = conv2(X,filtbio);
    data = data_tmp(:,1:T);

    for k = 1:size(W,2)
        W_tmp= conv2(squeeze(W(:,k,:)),filtbio);
        W(:,k,:) = W_tmp(:,1:L);
    end
else
    data = X;
end

end
