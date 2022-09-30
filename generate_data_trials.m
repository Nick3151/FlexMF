function [data,W,H,V_hat] = generate_data_trials(trials, frames, Nsequences, Nneurons, Dt, noise, jitter, participation, warp, bin, neg, seed)
% Generate data based on trials
if seed == 0
    rng shuffle
else
    rng(seed)
end

additional_neurons = 0;

% Parameters:
% trials = 30; % total number of trials
% frames = 150; % length of each trial
% number_of_seqences = 3;
% Nsequences = 6*ones(number_of_seqences, 1); % the number of occurences of sequences
% Nneurons = 10*ones(number_of_seqences, 1); % the number of neurons in each sequence
% Dt = 3.*ones(number_of_seqences,1); % gap between each member of the sequence
% noise = 0.001; % probability of added noise in each bin
% jitter = zeros(number_of_seqences,1); % Jitter time std
% participation = 1.*ones(number_of_seqences,1); % Participation probability = 100%
% warp = 0; % the maximum warping time
% bin = 0; % Binary data or not
% neg = 0; % Proportion of negative indices in W

assert(sum(Nsequences)<=trials, 'Total number of sequences must be less than trial numbers!')

%% Calculate useful things
N = sum(Nneurons)+additional_neurons; % Total number of neurons
K = length(Nsequences); % The number of sequences
lseq = Dt.*Nneurons; % the length of each sequences
lpad = 10; % zero padding before and after each sequence allows temporal warping and jittering
lseq_warp = (max(lseq)/max(Dt)*(max(Dt)+warp)); % warping seqeuence length
L = max(lseq)+2*lpad;   % Factor length
T = trials*frames; % Total data length
assert(max(lseq)<=frames, 'Sequence length must be less than or equal to the length of each trial!')

j = 1;
neurons = cell(K,1);
for k = 1:K
    neurons{k} = j:j+Nneurons(k)-1;
    j = j+Nneurons(k);
end

%% MAKE H's
H = zeros(K,T);
Hs = cell(K,1);
j = 1;
ind_perm = randperm(trials);
seq_ind = cell(K,1); % Index of trials in which sequence occurs
for k = 1:K
    seq_ind{k} = ind_perm(j:j+Nsequences(k)-1);
    j = j+Nsequences(k);
    if warp > 0
        Hs{k} = randi([-warp, warp],Nsequences(k),1);
    else
        Hs{k} = zeros(Nsequences(k),1);
    end
    for i=1:Nsequences(k)
        ind_trial = seq_ind{k}(i);
        H(k, (ind_trial-1)*frames+randi(frames-L)) = 1; % The time when sequence occurs within the time span of the trial
    end
end


%% Make W and V_hat
if seed == 0
    rng shuffle
else
    rng(seed+1)
end

W = zeros(N,K,L);
V_hat = zeros(N,T);
% H_hat = zeros(K,T);
for k = 1:K % go through each factor
    ind = find(H(k,:));

    % neg: proportion of negative indices
    neg_indices = (rand(1,Nneurons(k)) < neg);
    
    Dt_temp = Dt(k);
    Wk = zeros(N,L);          
    temp = ones(1,Nneurons(k));
    temp(neg_indices) = -temp(neg_indices);
    temp2 = zeros(length(neurons{k}),Dt_temp*Nneurons(k));
    temp2(:,1:Dt_temp:Dt_temp*Nneurons(k)) = diag(temp);    
    Wk(neurons{k},lpad:lpad-1+size(temp2,2)) = temp2;  
    W(:,k,:) = Wk;

    for j = 1:Nsequences(k) % go through each iteration of the sequence
        tempH = zeros(1,T);
        tempH(ind(j)) = 1;
                   
        if warp > 0 % change the dt for each instance
            Dt_temp = Dt(k)+Hs{k}(j);%+(randi(stretch))
            %*(-1+(2*(rand(1)>0.5))); If you want compression as well  
            L_temp = lseq_warp+2*lpad;
            tempW = zeros(N,L_temp);          
            temp2 = zeros(length(neurons{k}),Dt_temp*Nneurons(k));
            temp2(:,1:Dt_temp:Dt_temp*Nneurons(k)) = diag(temp);    
            tempW(neurons{k},lpad:lpad-1+size(temp2,2)) = temp2;  
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
    
%     % Compensate for factor shift due to warping
%     shift = circshift(1:length(tempH),floor((Dt(k)*Nneurons(k)- Dt_temp*Nneurons(k))/2));
%     tempH = tempH(:,shift);
%     H_hat(k,:) = H_hat(k,:) + tempH;
        newData = conv2(tempH,tempW);
        %V_hat = V_hat + newData(:,ceil(size(tempW,2)/2):end - floor(size(tempW,2)/2));
        V_hat = V_hat + newData(:,1:(end - size(tempW,2)+1));

    end

end

%% could add indepent noise later
V_hat = V_hat + (rand(size(V_hat))<noise);
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
    for k = 1:size(W,2)
        W(:,k,:) = conv2(squeeze(W(:,k,:)),filtbio,'same');
    end
end
rng shuffle

end