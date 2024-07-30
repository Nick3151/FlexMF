function [h, p, active_times] = test_significance_neurons(Wk, X_tmp, Xall, varargin)
p = inputParser;
addOptional(p, 'nNull', 100);
addOptional(p, 'seed', 0);

parse(p, varargin{:});
nNull = p.Results.nNull;
seed = p.Results.seed;
rng(seed)

assert(size(Xall,1)==size(Wk,1) && size(X_tmp,1)==size(Wk,1), 'Neurons do not meatch!')
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
h = zeros(N,1);
for ni=1:N
    [p(ni), h(ni)] = ranksum(reshape(overlap(ni,:,:), 1,[]), reshape(overlap_null(ni,:,:), 1,[]),...
        "tail","right", "alpha", .05/N);
end
active_times = sum(overlap>0,[2,3]);

