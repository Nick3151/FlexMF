function [true_pos_rates, false_pos_rates, ids] = match_motif_neurons_overlap(neurons_active_true, neurons_active_motifs, neurons_spont)
% Compute true positive rates and false positive rates for active motif neurons

K = length(neurons_active_true);
Khat = length(neurons_active_motifs);

true_pos_rates = zeros(Khat,1);
false_pos_rates = zeros(Khat,1);

% Match detected motifs to ground truth
for ii=1:K
    for jj=1:Khat
        neurons_motif_tmp = neurons_active_motifs{ii};
        neurons_true_tmp = neurons_active_true{jj};
        S(ii,jj) = length(intersect(neurons_motif_tmp, neurons_true_tmp))/length(neurons_true_tmp);
    end
end

temp = S;
ids = zeros(1, Khat);

for ii = 1:min(K,Khat)
    if ~any(temp(:))
        break
    end
    [r,c]= find(temp == max(temp(:)));
    true_pos_rates(c(1)) = temp(r(1), c(1));
    neurons_motif_tmp = neurons_active_motifs{c(1)};
    false_pos_rates(c(1)) = length(intersect(neurons_motif_tmp, neurons_spont))/length(neurons_spont);
    ids(c(1)) = r(1);

    temp(r(1),:) = 0;
    temp(:,c(1)) = 0;
end

% true_pos_rates(~ids) = [];
% false_pos_rates(~ids) = [];
% ids(~ids) = [];