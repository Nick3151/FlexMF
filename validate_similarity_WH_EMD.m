% VALIDATE_SIMILARITY_WH_EMD
%
% End-to-end check of the rewritten helper.similarity_WH_EMD on the two
% settings used by analyze_emd_cost_vs_lag, where the old version failed:
%
%   A  seqNMF with Khat=5 > K=3 on generate_data output. The old version let a
%      near-empty junk factor claim a real motif on the first greedy step,
%      because the cheapest point of every wrong pair's cost curve was the lag
%      that shifted the ground truth out of the window, and that cost scaled
%      with the estimate's own mass.
%   B  one exact copy of motif 1, one shifted duplicate of it, and one noise
%      factor. Motifs 2 and 3 have no counterpart, so a matcher with no
%      rejection rule is forced into two spurious pairs.
%
% Labels come from helper.similarity_WH, a shift-tolerant cosine similarity
% that involves no EMD, so the comparison is not circular.

clear

root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(fullfile(root, 'Utils'))
addpath(genpath(fullfile(root, 'FlexMF')))

K = 3;
T = 2000;
[X, W, H, ~] = generate_data(T, 5 * ones(K, 1), 3 * ones(K, 1), 'noise', 0, ...
    'participation', 0.8 * ones(K, 1), 'seed', 1, 'len_burst', 1, 'dynamic', 0);
[N, ~, L] = size(W);

Khat = 5;
frob = norm(X(:));
X = X / frob * Khat;
W = W / frob * Khat;

rng(0);
[WhatA, HhatA] = seqNMF(X, 'K', Khat, 'L', L, 'lambda', 0.005, ...
    'maxiter', 100, 'showPlot', 0);

WhatB = zeros(N, 3, L);
HhatB = zeros(3, T);
WhatB(:, 1, :) = W(:, 1, :);
HhatB(1, :) = H(1, :);
WhatB(:, 2, :) = helper.shift_profiles(reshape(W(:, 1, :), N, L), 1, L);
HhatB(2, :) = helper.shift_loadings(H(1, :), -1);
rng(1);
noiseProfile = rand(N, L) .* (rand(N, L) > 0.7);
WhatB(:, 3, :) = noiseProfile * mean(sum(W, [1 3])) / max(sum(noiseProfile(:)), eps);
HhatB(3, :) = double(rand(1, T) > 0.98);

settings = {'A: seqNMF, Khat=5 > K=3', WhatA, HhatA; ...
            'B: copy, duplicate, noise', WhatB, HhatB};

for s = 1:size(settings, 1)
    name = settings{s, 1};
    Wh = settings{s, 2};
    Hh = settings{s, 3};

    [~, ~, cosIds, cosDetails] = helper.similarity_WH(W, H ./ sum(H, 2), Wh, Hh);

    t0 = tic;
    [emdW, emdH, ids, det] = helper.similarity_WH_EMD(W, H, Wh, Hh);
    elapsed = toc(t0);

    fprintf('\n================ %s ================\n', name);
    fprintf('lag range per ground-truth motif:\n');
    for ii = 1:K
        if isempty(det.lagRange{ii})
            fprintf('  motif %d: skipped (zero factor)\n', ii);
        else
            fprintf('  motif %d: %d lags, %d:%d\n', ii, ...
                numel(det.lagRange{ii}), det.lagRange{ii}(1), det.lagRange{ii}(end));
        end
    end

    fprintf('\ncost of explaining each estimate from an empty profile:\n  %s\n', ...
        num2str(det.refCost, '%12.2f'));

    fprintf('\n%4s %4s %8s %12s %12s %9s %8s %9s\n', ...
        'i', 'j', 'cosine', 'match cost', 'cost/ref', 'best lag', 'reject', 'transport');
    for ii = 1:K
        for jj = 1:size(Wh, 2)
            fprintf('%4d %4d %8.3f %12.4f %12.4f %9s %8d %9s\n', ...
                ii, jj, cosDetails.S_pair(ii, jj), det.cost(ii, jj), ...
                det.cost(ii, jj) / det.refCost(jj), ...
                num2str(det.bestLag(ii, jj)), det.rejected(ii, jj), ...
                num2str(det.transport(ii, jj), '%.4f'));
        end
    end

    fprintf('\nmatching\n');
    fprintf('  ids from EMD    : %s\n', mat2str(ids));
    fprintf('  ids from cosine : %s\n', mat2str(cosIds));
    fprintf('  emds_W          : %s\n', num2str(emdW, '%10.4f'));
    fprintf('  emds_H          : %s\n', num2str(emdH, '%10.4f'));
    fprintf('  elapsed         : %.1f min\n', elapsed / 60);

    % Asked per ground-truth motif rather than per estimate: when a motif has
    % several convincing estimates, one-to-one matching is supposed to give it
    % to exactly one of them, so counting the others as disagreements would
    % penalise the behaviour being tested.
    agree = 0;
    checked = 0;
    for ii = 1:K
        candidates = find(cosDetails.S_pair(ii, :) > 0.8);
        if isempty(candidates)
            continue
        end
        checked = checked + 1;
        agree = agree + any(ids(candidates) == ii);
    end
    fprintf('  motifs given to one of their cosine counterparts: %d/%d\n', ...
        agree, checked);
end
