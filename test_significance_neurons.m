function [h, p, active_times] = test_significance_neurons(Wk, X_tmp, Xall, varargin)
p = inputParser;
addOptional(p, 'nNull', 100);
addOptional(p, 'seed', 0);
addOptional(p, 'alpha', .05);
% Multiple-comparison correction across neurons: 'bonferroni', 'bh' (Benjamini-Hochberg) or 'none'
addOptional(p, 'correction', 'bonferroni', @(x) any(strcmpi(x, {'bonferroni', 'bh', 'none'})));

parse(p, varargin{:});
nNull = p.Results.nNull;
seed = p.Results.seed;
alpha = p.Results.alpha;
correction = lower(p.Results.correction);
rng(seed)

assert(size(Xall,1)==size(Wk,1) && size(X_tmp,1)==size(Wk,1), 'Neurons do not match!')
[N,T] = size(Xall);
[~,L] = size(Wk);
R = size(X_tmp,2);  % number of motif occurrences

overlap = zeros(N,1,R);
overlap_null = zeros(N, nNull, R);

for r=1:R
    overlap(:,1,r) = sum(Wk.*squeeze(X_tmp(:,r,:)),2);
end

for n=1:nNull
    times_motifs_null = randperm(T-L,R);
    for r=1:R
        Xnull = Xall(:,times_motifs_null(r)+(1:L));
        overlap_null(:,n,r) = sum(Wk.*Xnull,2);
    end
end

%% Rank-sum test for each neuron
p = zeros(N,1);
for ni=1:N
    p(ni) = ranksum(reshape(overlap(ni,:,:), 1,[]), reshape(overlap_null(ni,:,:), 1,[]),...
        "tail","right");
end
switch correction
    case 'bonferroni'
        h = double(p<=alpha/N);
    case 'bh'
        % Reject the k smallest p values, with k the largest rank such that p_(k) <= k*alpha/N
        p_sorted = sort(p);
        k_max = find(p_sorted<=(1:N)'*alpha/N, 1, 'last');
        h = zeros(N,1);
        if ~isempty(k_max)
            h = double(p<=p_sorted(k_max));
        end
    case 'none'
        h = double(p<=alpha);
end
active_times = sum(overlap>0,[2,3]);

