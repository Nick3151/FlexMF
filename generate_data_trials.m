function [data,W,H,X_hat] = generate_data_trials(Trials, Length, Nmotifs, Nneurons, Magnitudes, Dt, noise, jitter, participation, warp, len_spikes, dynamic, neg, seed)
% Generate data based on trials
if seed == 0
    rng shuffle
else
    rng(seed)
end

additional_neurons = 0;

% Parameters:
% Trials = 30; % total number of trials
% Length = 150; % length of each trial
% number_of_seqences = 3;
% Nmotifs = 6*ones(number_of_seqences, 1); % the number of occurences of each motif
% Nneurons = 10*ones(number_of_seqences, 1); % the number of neurons in each motif
% Magnitudes = ones(number_of_seqences, 1); % the activation magnitudes of each motif
% Dt = 3.*ones(number_of_seqences,1); % gap between each member of the motif
% noise = 0.001; % probability of added noise in each bin
% jitter = zeros(number_of_seqences,1); % Jitter time std
% participation = 1.*ones(number_of_seqences,1); % Participation probability = 100%
% warp = 0; % the maximum warping time
% len_spikes = 20; The time of continuous firing
% dynamic = 1; if transient dynamics is being simulated
% neg = 0; % Proportion of negative indices in W

%% Calculate useful things
N = sum(Nneurons)+additional_neurons; % Total number of neurons
K = length(Nmotifs); % The number of motifs
lmotif = Dt.*Nneurons+1; % the length of each motif
if Dt>0
    lseq_warp = (max(lmotif)/max(Dt)*(max(Dt)+warp)); % warping motif length
else
    lseq_warp = 1;
end
assert(lseq_warp<=Length, 'Sequence length must be less than or equal to the length of each trial!')

j = 1;
neurons = cell(K,1);
for k = 1:K
    neurons{k} = j:j+Nneurons(k)-1;
    j = j+Nneurons(k);
end

%% MAKE H's, sample warping params
H = zeros(K,Trials);
Hs = cell(K,1);

motif_ind = cell(K,1); % Index of trials in which motif occurs
for k = 1:K
    motif_ind{k} = randperm(Trials, Nmotifs(k));
    if warp > 0
        Hs{k} = randi([-warp, warp],Nmotifs(k),1);
    else
        Hs{k} = zeros(Nmotifs(k),1);
    end
    H(k, motif_ind{k}) = ones(1,Nmotifs(k))*Magnitudes(k);
end


%% Make W and X_hat
if seed == 0
    rng shuffle
else
    rng(seed+1)
end

W = zeros(N,K,Length);
X_hat = zeros(N,Length,Trials);
% H_hat = zeros(K,T);
for k = 1:K % go through each factor

    % neg: proportion of negative indices
    neg_indices = (rand(1,Nneurons(k)) < neg);
    
    Dt_temp = Dt(k);
    Wk = zeros(N,Length);          
    temp = ones(1,Nneurons(k));
    temp(neg_indices) = -temp(neg_indices);
    if Dt_temp > 0
        l = Dt_temp*Nneurons(k);
        temp2 = zeros(length(neurons{k}),l);
        temp2(:,1:Dt_temp:l) = diag(temp);
        Wk(neurons{k},(1:l)+2*max(jitter)) = temp2;  
    else
        Wk(neurons{k},1+2*max(jitter)) = temp';
    end
    
    W(:,k,:) = Wk;

    for j = 1:Nmotifs(k) % go through each iteration of the sequence
                   
        if warp > 0 % change the dt for each instance
            tempW = zeros(N,Length); 
            Dt_temp = Dt(k)+Hs{k}(j); 
            if Dt_temp > 0
                l = Dt_temp*Nneurons(k);         
                temp2 = zeros(length(neurons{k}),l);
                temp2(:,1:Dt_temp:l) = diag(temp);    
                tempW(neurons{k},(1:l)+2*max(jitter)) = temp2;  
            else
                tempW(neurons{k},1+2*max(jitter)) = temp';
            end
        else 
            tempW = Wk;      
        end

        % neurons are jittered with some time bins
        %shifts = poissrnd(SeqNoiseTime(ii),N,1).*(1-2*(rand(N,1)>0.5));
        shifts = round(normrnd(0,jitter(k),N,1));
        %shifts(abs(shifts) >5) = 0; %stop poiss from getting to big
        for idx = 1:N
            tempW(idx,:) = circshift(tempW(idx,:),shifts(idx)*1);
        end
        
        % neurons participate with some p
        tempW(rand(N,1)>participation(k),:) = 0;
        X_hat(:,:,motif_ind{k}(j)) = X_hat(:,:,motif_ind{k}(j)) + H(k,motif_ind{k}(j))*tempW;
    end

end

%% Add indepent noise
X_noise = (rand(size(X_hat))<noise);

%% Continuous firing
if isempty(len_spikes)
    len_spikes = 1;
end

for t = 1:Trials
    X_hat_tmp = conv2(squeeze(X_hat(:,:,t)), ones(1,len_spikes));
    X_hat(:,:,t) = X_hat_tmp(:,1:Length);
end

for k = 1:K
    W_tmp = conv2(squeeze(W(:,k,:)), ones(1,len_spikes));
    W(:,k,:) = W_tmp(:,1:Length);
end

X = X_hat + (~X_hat).*X_noise;
X(isnan(X)) = 0;

%% Include calcimu dynamic
if dynamic
    filtbio = exp(-(0:8)/4); 
    data = zeros(N,Length,Trials);
    for t = 1:Trials
        data_tmp = conv2(squeeze(X(:,:,t)),filtbio);
        data(:,:,t) = data_tmp(:,1:Length);
    end
    for k = 1:K
        W_tmp= conv2(squeeze(W(:,k,:)),filtbio);
        W(:,k,:) = W_tmp(:,1:Length);
    end
else
    data = X;
end

rng shuffle

end