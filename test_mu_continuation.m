function results = test_mu_continuation(varargin)
%TEST_MU_CONTINUATION  Sweep TFOCS mu and muDecrement for FlexMF EMD updates.
%
%   results = test_mu_continuation
%   results = test_mu_continuation('quick', true)   % default; smaller T
%   results = test_mu_continuation('quick', false)  % T=2000, denser grid
%
% For each (mu, muDecrement) pair:
%   1) Probe one updateH_EMD + updateW_EMD from a shared SeqNMF warm-start
%      (constraint residual, niter, status / iteration-limit flag)
%   2) Run a short FlexMF fit (no Reweight / no TV) and record constraint,
%      L1(M)/||X||_1, mean EMD(W), and wall time
%
% Scores candidates (lower better) and prints the optimal (mu, muDecrement).

%% -------- Options --------
ip = inputParser;
addParameter(ip, 'quick', true, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'seed', 1, @isnumeric);
addParameter(ip, 'do_save', true, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'outDir', 'Simulation_Results', @(s) ischar(s) || isstring(s));
parse(ip, varargin{:});
opt = ip.Results;
opt.quick = logical(opt.quick);
opt.do_save = logical(opt.do_save);
opt.outDir = char(opt.outDir);
if ~exist(opt.outDir, 'dir'), mkdir(opt.outDir); end

%% -------- Paths --------
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
root = fileparts(thisDir);
addpath(fullfile(root, 'TFOCS'))
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))
if exist(fullfile(root, 'seqNMF-master'), 'dir')
    rmpath(genpath(fullfile(root, 'seqNMF-master')));
end
addpath(genpath(thisDir));

%% -------- Data / FlexMF settings (match dynamic demo; no Reweight/TV) --------
K = 3;
Khat = 3;
Nneurons = 10*ones(K,1);
Dt = 3.*ones(K,1);
noise = 0.001;
warp = 2;
len_burst = 5;
Lhat_min = 50;
lambda = 1e-2;
lambda_M = 0.1;
lambda_R = 1;
tolerance = 1e-3;
do_normalize = true;  % scale data Frobenius norm to Khat

if opt.quick
    T = 800;
    maxiter = 5;
    mus = [0.01, 0.1, 1, 10];
    muDecrements = [1, 0.5, 0.2];
    fprintf('Mode: QUICK (T=%d, maxiter=%d)\n', T, maxiter);
else
    T = 2000;
    maxiter = 10;
    mus = [0.01, 0.1, 1, 10, 100];
    muDecrements = [1, 0.5, 0.2, 0.1];
    fprintf('Mode: FULL (T=%d, maxiter=%d)\n', T, maxiter);
end

%% -------- Generate data + SeqNMF warm-start (once) --------
fprintf('Generating warp+noise+burst+dynamics data...\n');
[X, Wtrue, Htrue, ~] = generate_data(T, Nneurons, Dt, ...
    'noise', noise, 'warp', warp, 'len_burst', len_burst, 'dynamic', 1, ...
    'seed', opt.seed);
[N, T] = size(X);
L = size(Wtrue, 3);
Lhat = max(L, Lhat_min);
fprintf('  size(X)=[%d %d], L=%d, Lhat=%d\n', N, T, L, Lhat);

if do_normalize
    frob_norm = norm(X(:));
    X = X / frob_norm * Khat;
    Wtrue = Wtrue / frob_norm * Khat;
    fprintf('  Normalized ||X||_F -> %g\n', Khat);
end

fprintf('SeqNMF warm-start (lambda=%g)...\n', lambda);
tSeq = tic;
[W0, H0] = seqNMF(X, 'K', K, 'L', Lhat, 'lambda', lambda, ...
    'maxiter', 50, 'showPlot', 0);
fprintf('  SeqNMF done in %.1f s\n', toc(tSeq));

%% -------- Sweep --------
nMu = numel(mus);
nDec = numel(muDecrements);
constraint_rel = nan(nMu, nDec);
constraint_H = nan(nMu, nDec);
constraint_W = nan(nMu, nDec);
L1M_rel = nan(nMu, nDec);
cost_final = nan(nMu, nDec);
time_s = nan(nMu, nDec);
niter_H = nan(nMu, nDec);
niter_W = nan(nMu, nDec);
hit_limit_H = false(nMu, nDec);
hit_limit_W = false(nMu, nDec);
status_H = cell(nMu, nDec);
status_W = cell(nMu, nDec);
emds_W_mean = nan(nMu, nDec);

normX1 = norm(X(:), 1);
flexBase = {'K', K, 'L', Lhat, 'EMD', 1, ...
    'lambda', lambda, 'lambda_M', lambda_M, 'lambda_R', lambda_R, ...
    'lambdaL1H', 0, 'lambda_TV', 0, 'Reweight', 0, ...
    'maxiter', maxiter, 'tolerance', tolerance, 'neg_prop', 0, ...
    'W_init', W0, 'H_init', H0, 'showPlot', 0, 'verbal', 0};

fprintf('\nSweeping %d x %d = %d settings...\n', nMu, nDec, nMu*nDec);
for i = 1:nMu
    for j = 1:nDec
        mu = mus(i);
        muDec = muDecrements(j);
        fprintf('\n=== mu=%g, muDecrement=%g ===\n', mu, muDec);

        % Params for a single H/W probe (same fields update*_EMD expects)
        params0 = struct( ...
            'lambda', lambda, 'lambda_M', lambda_M, 'lambda_R', lambda_R, ...
            'lambdaL1H', 0, 'lambdaL1W', 0, 'lambda_TV', 0, ...
            'Reweight', 0, 'homotopy', 10, 'currentiter', 1, ...
            'mu', mu, 'muDecrement', muDec, 'verbal', 0, 'alg', 'N83');

        M00 = zeros(N, T);
        R00 = zeros(N, T);

        [Hp, Mp, Rp, outH] = updateH_EMD(W0, H0, X, M00, R00, params0);
        Xcorr = helper.correct_warp(X, Mp);
        resid = Xcorr - Rp - helper.reconstruct(W0, Hp);
        constraint_H(i,j) = norm(resid(:), 1) / normX1;
        niter_H(i,j) = getfield_default(outH, 'niter', NaN);
        status_H{i,j} = getfield_default(outH, 'status', '');
        hit_limit_H(i,j) = contains_str(status_H{i,j}, 'iteration limit');

        [Wp, Mp2, Rp2, outW] = updateW_EMD(W0, Hp, X, Mp, Rp, params0);
        Xcorr = helper.correct_warp(X, Mp2);
        resid = Xcorr - Rp2 - helper.reconstruct(Wp, Hp);
        constraint_W(i,j) = norm(resid(:), 1) / normX1;
        niter_W(i,j) = getfield_default(outW, 'niter', NaN);
        status_W{i,j} = getfield_default(outW, 'status', '');
        hit_limit_W(i,j) = contains_str(status_W{i,j}, 'iteration limit');
        fprintf('  probe H: constr_rel=%.3g  niter=%g  [%s]\n', ...
            constraint_H(i,j), niter_H(i,j), status_H{i,j});
        fprintf('  probe W: constr_rel=%.3g  niter=%g  [%s]\n', ...
            constraint_W(i,j), niter_W(i,j), status_W{i,j});

        t0 = tic;
        [What, Hhat, cost, ~, ~, ~, M, R] = FlexMF(X, flexBase{:}, ...
            'mu', mu, 'muDecrement', muDec);
        time_s(i,j) = toc(t0);

        Xcorr = helper.correct_warp(X, M);
        resid = Xcorr - R - helper.reconstruct(What, Hhat);
        constraint_rel(i,j) = norm(resid(:), 1) / normX1;
        L1M_rel(i,j) = norm(M(:), 1) / normX1;
        cost_final(i,j) = cost(end);
        [eW, ~, ids] = helper.similarity_WH_EMD(Wtrue, Htrue, What, Hhat);
        emds_W_mean(i,j) = mean(eW, 'omitnan');
        fprintf('  FlexMF: constr_rel=%.3g  L1M/X=%.3g  EMD_W=%.3g  n_det=%d  time=%.1fs\n', ...
            constraint_rel(i,j), L1M_rel(i,j), emds_W_mean(i,j), sum(ids > 0), time_s(i,j));
    end
end

%% -------- Score and pick optimum --------
hit = double(hit_limit_H | hit_limit_W);
score = log10(constraint_rel + 1e-12) ...
    + 5 * hit ...
    + 0.5 * log10(L1M_rel + 1e-12) ...
    + 0.05 * log10(time_s + 1e-12);

[bestScore, idx] = min(score(:));
[iBest, jBest] = ind2sub(size(score), idx);
mu_star = mus(iBest);
muDec_star = muDecrements(jBest);

fprintf('\n========== RESULTS ==========\n');
fprintf('constraint_rel:\n'); disp(constraint_rel);
fprintf('hit iteration limit (H|W):\n'); disp(hit_limit_H | hit_limit_W);
fprintf('L1_M / ||X||_1:\n'); disp(L1M_rel);
fprintf('time (s):\n'); disp(time_s);
fprintf('score (lower better):\n'); disp(score);
fprintf('OPTIMAL: mu = %g, muDecrement = %g  (score=%.3g)\n', ...
    mu_star, muDec_star, bestScore);
fprintf('  constraint_rel=%.4g  L1M/X=%.4g  EMD_W=%.4g  time=%.1fs\n', ...
    constraint_rel(iBest,jBest), L1M_rel(iBest,jBest), ...
    emds_W_mean(iBest,jBest), time_s(iBest,jBest));
fprintf('  probe H status: %s\n', status_H{iBest,jBest});
fprintf('  probe W status: %s\n', status_W{iBest,jBest});
fprintf('==============================\n');
fprintf('Suggested FlexMF call:\n');
fprintf('  FlexMF(..., ''mu'', %g, ''muDecrement'', %g)\n', mu_star, muDec_star);

%% -------- Refit + visualize WH / MR under optimal params --------
fprintf('\nRefitting FlexMF with optimal mu=%g, muDecrement=%g ...\n', mu_star, muDec_star);
[What_star, Hhat_star, cost_star, ~, ~, ~, M_star, R_star] = FlexMF(X, flexBase{:}, ...
    'mu', mu_star, 'muDecrement', muDec_star, 'showPlot', 0, 'verbal', 0);
Xcorr_star = helper.correct_warp(X, M_star);
resid_star = Xcorr_star - R_star - helper.reconstruct(What_star, Hhat_star);
constr_star = norm(resid_star(:), 1) / normX1;
fprintf('  optimal refit: constraint_rel=%.4g  L1M/X=%.4g  cost=%.4g\n', ...
    constr_star, norm(M_star(:), 1) / normX1, cost_star(end));

figure; SimpleWHPlot_patch(What_star, Hhat_star, 'plotAll', 1);
title(sprintf('FlexMF WH (mu=%g, muDec=%g)', mu_star, muDec_star))
set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8])
if opt.do_save
    save2pdf(fullfile(opt.outDir, sprintf('test_mu_opt_WH_mu=%g_muDec=%g.pdf', mu_star, muDec_star)), gcf);
end

figure; SimpleWHPlot_patch(What_star, Hhat_star, 'Data', X, 'plotAll', 1, 'compare', true);
title(sprintf('FlexMF WH + data (mu=%g, muDec=%g)', mu_star, muDec_star))
set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8])
if opt.do_save
    save2pdf(fullfile(opt.outDir, sprintf('test_mu_opt_WH_data_mu=%g_muDec=%g.pdf', mu_star, muDec_star)), gcf);
end

figure; plot_MR(M_star, R_star, sprintf('mu=%g, muDec=%g', mu_star, muDec_star));
if opt.do_save
    save2pdf(fullfile(opt.outDir, sprintf('test_mu_opt_MR_mu=%g_muDec=%g.pdf', mu_star, muDec_star)), gcf);
end

%% -------- Plots --------
figure;
imagesc(log10(constraint_rel + 1e-12)); colorbar
set(gca, 'XTick', 1:nDec, 'XTickLabel', muDecrements, ...
    'YTick', 1:nMu, 'YTickLabel', mus);
xlabel('muDecrement'); ylabel('mu');
title('log_{10}(constraint\_rel)'); hold on
plot(jBest, iBest, 'w*', 'MarkerSize', 14, 'LineWidth', 1.5);
if opt.do_save
    save2pdf(fullfile(opt.outDir, 'test_mu_constraint_rel.pdf'), gcf);
end

figure;
imagesc(double(hit_limit_H | hit_limit_W)); colorbar
set(gca, 'XTick', 1:nDec, 'XTickLabel', muDecrements, ...
    'YTick', 1:nMu, 'YTickLabel', mus, 'CLim', [0 1]);
xlabel('muDecrement'); ylabel('mu');
title('Inner solve hit iteration limit (H or W)');
if opt.do_save
    save2pdf(fullfile(opt.outDir, 'test_mu_hit_limit.pdf'), gcf);
end

figure;
imagesc(score); colorbar
set(gca, 'XTick', 1:nDec, 'XTickLabel', muDecrements, ...
    'YTick', 1:nMu, 'YTickLabel', mus);
xlabel('muDecrement'); ylabel('mu');
title(sprintf('Score (best: mu=%g, muDec=%g)', mu_star, muDec_star));
hold on; plot(jBest, iBest, 'w*', 'MarkerSize', 14, 'LineWidth', 1.5);
if opt.do_save
    save2pdf(fullfile(opt.outDir, 'test_mu_score.pdf'), gcf);
end

%% -------- Pack / save --------
results = struct();
results.mus = mus;
results.muDecrements = muDecrements;
results.constraint_rel = constraint_rel;
results.constraint_H = constraint_H;
results.constraint_W = constraint_W;
results.L1M_rel = L1M_rel;
results.cost_final = cost_final;
results.time_s = time_s;
results.niter_H = niter_H;
results.niter_W = niter_W;
results.hit_limit_H = hit_limit_H;
results.hit_limit_W = hit_limit_W;
results.status_H = status_H;
results.status_W = status_W;
results.emds_W_mean = emds_W_mean;
results.score = score;
results.mu_star = mu_star;
results.muDecrement_star = muDec_star;
results.bestScore = bestScore;
results.What_star = What_star;
results.Hhat_star = Hhat_star;
results.M_star = M_star;
results.R_star = R_star;
results.constraint_rel_star = constr_star;
results.quick = opt.quick;
results.T = T;
results.maxiter = maxiter;
results.lambda = lambda;
results.lambda_M = lambda_M;
results.lambda_R = lambda_R;
results.X = X;
results.Wtrue = Wtrue;
results.Htrue = Htrue;
results.W0 = W0;
results.H0 = H0;

if opt.do_save
    matFile = fullfile(opt.outDir, sprintf('test_mu_continuation_T=%d.mat', T));
    save(matFile, '-struct', 'results');
    fprintf('Saved %s\n', matFile);
end

end

%% -------- Local helpers --------
function v = getfield_default(s, name, default)
if isstruct(s) && isfield(s, name)
    v = s.(name);
else
    v = default;
end
end

function tf = contains_str(s, pattern)
s = lower(char(s));
tf = ~isempty(strfind(s, lower(pattern))); %#ok<STREMP>
end
