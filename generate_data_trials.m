function [data,W,H,X_hat] = generate_data_trials(Trials, Length, Nmotifs, Nneurons, Magnitudes, Dt, noise, jitter, participation, warp, bin, neg, seed)
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
% bin = 0; % Binary data or not
% neg = 0; % Proportion of negative indices in W

%% Calculate useful things
N = sum(Nneurons)+additional_neurons; % Total number of neurons
K = length(Nmotifs); % The number of motifs
lmotif = Dt.*Nneurons; % the length of each motif

lseq_warp = (max(lmotif)/max(Dt)*(max(Dt)+warp)); % warping motif length
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
    l = Dt_temp*Nneurons(k);
    temp2 = zeros(length(neurons{k}),l);
    temp2(:,1:Dt_temp:l) = diag(temp);    
    Wk(neurons{k},ceil((Length-l)/2):ceil((Length-l)/2)-1+l) = temp2;  
    W(:,k,:) = Wk;

    for j = 1:Nmotifs(k) % go through each iteration of the sequence
                   
        if warp > 0 % change the dt for each instance
            Dt_temp = Dt(k)+Hs{k}(j); 
            l = Dt_temp*Nneurons(k);
            tempW = zeros(N,Length);          
            temp2 = zeros(length(neurons{k}),l);
            temp2(:,1:Dt_temp:l) = diag(temp);    
            tempW(neurons{k},ceil((Length-l)/2):ceil((Length-l)/2)-1+l) = temp2;  
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

%% could add indepent noise later
X_noise = (rand(size(X_hat))<noise);
X_hat = X_hat + (~X_hat).*X_noise;
X_hat(isnan(X_hat)) = 0;
%%
%data = V_hat;
if ~bin
    filtbio = exp(-(1:30)/10); 
    data = zeros(N,Length,Trials);
    for t = 1:Trials
        data(:,:,t) = conv2(squeeze(X_hat(:,:,t)),filtbio,'same');
    end
else
    data = X_hat;
end

if ~bin
    for k = 1:K
        W(:,k,:) = conv2(squeeze(W(:,k,:)),filtbio,'same');
    end
end
rng shuffle

end