function [data,W,H,X_hat,motif_ind] = generate_data_trials(Trials, Length, Nmotifs, Nneurons, Dt, varargin)
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
% overlap_t = 0; if temporal overlap between different motifs is allowed
% overlap_n = 0; The proportion of shared neurons between different motifs
% len_burst = 10; The time of continuous firing
% dynamic = 1; if transient dynamics is being simulated
% neg = 0; % Proportion of negative indiuntitledces in Wuntitled
% seed = 0; % Random seed
% additional_neurons = 0; % Number of additional neurons
% spontaneous_rate = 0; % Spontaneous activities of additional neurons

%% Parse inputs
K = length(Nmotifs); % The number of motifs
if size(Nmotifs,1) ~= 1
    Nmotifs = Nmotifs';
end
validArrayK = @(x) isnumeric(x) && length(x)==K;

p  = inputParser;
addOptional(p, 'Magnitudes', ones(K, 1), validArrayK)
addOptional(p, 'noise', .001, @isscalar)
addOptional(p, 'jitter', zeros(K, 1), validArrayK)
addOptional(p, 'participation', ones(K, 1), validArrayK)
addOptional(p, 'warp', 0, @isscalar)
addOptional(p, 'overlap_t', 0, @isscalar)
addOptional(p, 'overlap_n', 0, @isscalar)
addOptional(p, 'len_burst', 10, @isscalar)
addOptional(p, 'dynamic', 1, @isscalar)
addOptional(p, 'neg', 0, @isscalar)
addOptional(p, 'seed', 0, @isscalar)
addOptional(p, 'additional_neurons', 0, @isscalar)
addOptional(p, 'spontaneous_rate', 0, @isscalar)
parse(p, varargin{:})

Magnitudes = p.Results.Magnitudes;
noise = p.Results.noise;
jitter = p.Results.jitter;
participation = p.Results.participation;
warp = p.Results.warp;
overlap_t = p.Results.overlap_t;
overlap_n = p.Results.overlap_n;
len_burst = p.Results.len_burst;
dynamic = p.Results.dynamic;
neg = p.Results.neg;
seed = p.Results.seed;
additional_neurons = p.Results.additional_neurons;
spontaneous_rate = p.Results.spontaneous_rate;

% Set random seed
if seed == 0
    rng shuffle
else
    rng(seed)
end

%% Calculate useful things
lmotif = Dt.*Nneurons+1; % the length of each motif
if Dt>0
    lseq_warp = (max(lmotif)/max(Dt)*(max(Dt)+warp)); % warping motif length
else
    lseq_warp = 1;
end
assert(lseq_warp<=Length, 'Sequence length must be less than or equal to the length of each trial!')

j = 1;
neurons = cell(K,1);
assert((overlap_n<1) && (overlap_n>=0), 'Overlap_n should be between [0,1)!')
num_overlap = floor(overlap_n*Nneurons);
for k = 1:K
    neurons{k} = j:j+Nneurons(k)-1;
    j = j+Nneurons(k)-num_overlap(k);
end
N = j+num_overlap(K)-1+additional_neurons; % Total number of neurons

%% MAKE H's, sample warping params
H = zeros(K,Trials);
Hs = cell(K,1);

motif_ind = cell(K,1); % Index of trials in which motif occurs

if overlap_t
    for k = 1:K
        motif_ind{k} = randperm(Trials, Nmotifs(k));
        if warp > 0
            Hs{k} = randi([-warp, warp],Nmotifs(k),1);
        else
            Hs{k} = zeros(Nmotifs(k),1);
        end
        H(k, motif_ind{k}) = ones(1,Nmotifs(k))*Magnitudes(k);
    end
else
    assert(sum(Nmotifs)<=Trials, 'Do not have enough trials!')
    motif_inds_all = randperm(Trials, sum(Nmotifs));
    inds_tmp = [0,cumsum(Nmotifs)];
    for  k=1:K
        motif_ind{k} = motif_inds_all(inds_tmp(k)+1:inds_tmp(k+1));
        if warp > 0
            Hs{k} = randi([-warp, warp],Nmotifs(k),1);
        else
            Hs{k} = zeros(Nmotifs(k),1);
        end
        H(k, motif_ind{k}) = ones(1,Nmotifs(k))*Magnitudes(k);
    end
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
for k = 1:K % go through each motif

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

%% Add independent noise
X_noise = (rand(size(X_hat))<noise);

%% Spontaneous activities of additional neurons
if additional_neurons
    X_spont = (rand(additional_neurons,Length,Trials)<spontaneous_rate);
    X_hat(N-additional_neurons+1:N,:,:) = X_spont;
end

%% Burst
if isempty(len_burst)
    len_burst = 1;
end

for t = 1:Trials
    X_hat_tmp = conv2(squeeze(X_hat(:,:,t)), ones(1,len_burst));
    X_hat(:,:,t) = X_hat_tmp(:,1:Length);
    if additional_neurons
        X_spont_tmp = conv2(squeeze(X_spont(:,:,t)), ones(1,len_burst));
        X_spont(:,:,t) = X_spont_tmp(:,1:Length);
    end
end

for k = 1:K
    W_tmp = conv2(squeeze(W(:,k,:)), ones(1,len_burst));
    W(:,k,:) = W_tmp(:,1:Length);
end

X = X_hat;
if additional_neurons
    X(N-additional_neurons+1:N,:,:) = X_spont;
end
X = X + (~X).*X_noise;
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