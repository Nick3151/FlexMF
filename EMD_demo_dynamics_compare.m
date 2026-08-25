%% Demo: SeqNMF vs FlexMF (EMD) under warping/jittering + additive noise
% Choose temporal distortion (jitter or warp), always with noise, plus
% optional burst length and calcium dynamics. SeqNMF warm-starts FlexMF.
% Compares only:
%   1) SeqNMF
%   2) FlexMF (no Reweight, no TV)
% Validates the EMD constraint after the FlexMF fit, and plots WH / MR.
clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))
rmpath(genpath(fullfile(root, 'seqNMF-master')));
addpath(genpath(fullfile(root, 'FlexMF')));

%% -------- User settings --------
temporal = 'warp';      % 'jitter' | 'warp'
use_burst = true;       % len_burst > 1
use_dynamics = false;    % calcium / transient filter
do_choose_lambda = true;  % true: sweep lambdas for SeqNMF before fitting
do_normalize = true;      % scale data Frobenius norm to Khat
do_save = true;
constraintTol = 0.05;   % flag if ||constraint||_1 / ||X||_1 exceeds this

K = 3;
Khat = 5;
T = 2000;
Nneurons = 10*ones(K,1);
Dt = 3.*ones(K,1);
noise = 0.005;
jitter = 5*ones(K,1);
warp = 8;
len_burst = 3;          % used when use_burst
seed = 1;
maxiter = 50;
Lhat_min = 50;          % longer Lhat helps when dynamics/burst smear motifs

% FlexMF warm-starts from SeqNMF; each method has its own lambda
lambda_seqNMF = 1e-2;
lambda_FlexMF = 5e-2;
lambda_M = 0.05;
lambda_R = 1;
tolerance = 1e-3;
% TFOCS SCD (from test_mu_continuation: best on warp+noise+burst+dynamics)
mu = .1;
muDecrement = 1;

%% -------- Build data config --------
assert(ismember(temporal, {'jitter', 'warp'}), ...
    'temporal must be ''jitter'' or ''warp''.');

gen_args = {'noise', noise, 'seed', seed};
if strcmp(temporal, 'jitter')
    gen_args = [gen_args, {'jitter', jitter}];
else
    gen_args = [gen_args, {'warp', warp}];
end
if use_burst
    gen_args = [gen_args, {'len_burst', len_burst}];
else
    gen_args = [gen_args, {'len_burst', 1}];
end
if use_dynamics
    gen_args = [gen_args, {'dynamic', 1}];
else
    gen_args = [gen_args, {'dynamic', 0}];
end

data_tag = temporal;
data_tag = [data_tag '+noise'];
if use_burst, data_tag = [data_tag '+burst']; end
if use_dynamics, data_tag = [data_tag '+dynamics']; end
data_file = sanitize_name(data_tag);

%% -------- Generate data --------
[X, W, H, ~] = generate_data(T, Nneurons, Dt, gen_args{:});
L = size(W, 3);
Lhat = max(L, Lhat_min);
fprintf('Data: %s | size(X)=[%d %d] | L=%d | Lhat=%d\n', ...
    data_tag, size(X,1), size(X,2), L, Lhat);

%% -------- Normalize --------
if do_normalize
    frob_norm = norm(X(:));
    X = X / frob_norm * Khat;
    W = W / frob_norm * Khat;
    fprintf('Normalized ||X||_F -> %g\n', Khat);
end

plotAll = 1;
figure; SimpleWHPlot_patch(W, H, 'Data', X, 'plotAll', plotAll);
title(sprintf('Generated data (%s)', data_tag), 'FontSize', 16)
set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8])
if do_save
    save2pdf(sprintf('EMD_simulated_data_%s.pdf', data_file))
end

%% -------- Optional SeqNMF lambda sweep --------
if do_choose_lambda
    nLambdas = 20;
    lambdas = sort(logspace(0, -4, nLambdas), 'ascend');
    regularization = zeros(1, nLambdas);
    cost = zeros(1, nLambdas);
    for li = 1:nLambdas
        [What_tmp, Hhat_tmp] = seqNMF(X, 'K', Khat, 'L', Lhat, ...
            'lambda', lambdas(li), 'maxiter', 100, 'showPlot', 0);
        [cost(li), regularization(li), ~] = helper.get_seqNMF_cost(X, What_tmp, Hhat_tmp);
        fprintf('Testing lambda %d/%d (lambda=%g)\n', li, nLambdas, lambdas(li));
    end

    windowSize = 3;
    b = (1/windowSize)*ones(1, windowSize);
    a = 1;
    Rs = filtfilt(b, a, regularization);
    minRs = prctile(regularization, 10); maxRs = prctile(regularization, 90);
    Rs = (Rs - minRs) / (maxRs - minRs);
    Rsc = (regularization - minRs) / (maxRs - minRs);
    Cs = filtfilt(b, a, cost);
    minCs = prctile(cost, 10); maxCs = prctile(cost, 90);
    Cs = (Cs - minCs) / (maxCs - minCs);
    Csc = (cost - minCs) / (maxCs - minCs);

    figure; hold on
    plot(lambdas, Rs, 'b')
    plot(lambdas, Cs, 'r')
    scatter(lambdas, Rsc, 'b', 'markerfacecolor', 'flat');
    scatter(lambdas, Csc, 'r', 'markerfacecolor', 'flat');
    xlabel('Lambda'); ylabel('Cost (au)')
    set(legend('Correlation cost', 'Reconstruction cost'), 'Box', 'on')
    set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
    set(gca, 'color', 'none', 'tickdir', 'out', 'ticklength', [0.025, 0.025])
    title(sprintf('SeqNMF lambda sweep (%s)', data_tag))
    if do_save
        save2pdf(sprintf('Simulate_%s_choose_lambda_SeqNMF', data_file))
    end
end

%% -------- Method definitions --------
% Only SeqNMF and baseline FlexMF (no Reweight, no TV); FlexMF warm-starts from SeqNMF
method_names = {'SeqNMF', 'FlexMF (no Reweight, no TV)'};
nMethod = numel(method_names);

What = cell(nMethod, 1);
Hhat = cell(nMethod, 1);
Mhat = cell(nMethod, 1);
Rhat = cell(nMethod, 1);
times = nan(nMethod, 1);
constraint_rel = nan(nMethod, 1);  % SeqNMF left NaN (no M,R)

%% -------- 1) SeqNMF --------
fprintf('\n=== %s (lambda=%g) ===\n', method_names{1}, lambda_seqNMF);
figure;
set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8])
t0 = tic;
[What{1}, Hhat{1}, ~, ~, ~, ~] = seqNMF(X, 'K', Khat, 'L', Lhat, ...
    'lambda', lambda_seqNMF, 'maxiter', maxiter, 'showPlot', 1);
times(1) = toc(t0);
fprintf('  time: %.2f s\n', times(1));

plot_WH(What{1}, Hhat{1}, X, method_names{1}, plotAll);
if do_save
    save2pdf(sprintf('Simulated_%s_SeqNMF_WH.pdf', data_file), gcf)
end

%% -------- 2) FlexMF (no Reweight, no TV; SeqNMF warm-start) --------
fprintf('\n=== %s ===\n', method_names{2});
fprintf('  lambda=%g, lambda_M=%g, lambda_R=%g, mu=%g, muDecrement=%g\n', ...
    lambda_FlexMF, lambda_M, lambda_R, mu, muDecrement);

figure;
set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8])
t0 = tic;
[What{2}, Hhat{2}, ~, ~, ~, ~, Mhat{2}, Rhat{2}] = FlexMF(X, ...
    'K', Khat, 'L', Lhat, 'EMD', 1, ...
    'lambda', lambda_FlexMF, 'lambda_M', lambda_M, 'lambda_R', lambda_R, ...
    'lambdaL1H', 0, 'lambda_TV', 0, 'Reweight', 0, ...
    'mu', mu, 'muDecrement', muDecrement, ...
    'maxiter', maxiter, 'tolerance', tolerance, 'neg_prop', 0, ...
    'W_init', What{1}, 'H_init', Hhat{1});
times(2) = toc(t0);
fprintf('  time: %.2f s\n', times(2));

% Constraint: X_corr - R ≈ W(*)H
Xcorr = helper.correct_warp(X, Mhat{2});
constraint = Xcorr - Rhat{2} - helper.reconstruct(What{2}, Hhat{2});
constraint_rel(2) = norm(constraint(:), 1) / norm(X(:), 1);
fprintf('  constraint ||.||_1 / ||X||_1 = %.4g', constraint_rel(2));
if constraint_rel(2) > constraintTol
    fprintf('  [FLAG > %.3g]\n', constraintTol);
else
    fprintf('  [OK]\n');
end

plot_WH(What{2}, Hhat{2}, X, method_names{2}, plotAll);
if do_save
    save2pdf(sprintf('Simulated_%s_FlexMF_WH.pdf', data_file), gcf)
end

figure; plot_MR(Mhat{2}, Rhat{2}, method_names{2});
if do_save
    save2pdf(sprintf('Simulated_%s_FlexMF_MR.pdf', data_file), gcf)
end

figure; imagesc(constraint); colorbar
title(sprintf('Constraint residual — %s (rel=%.3g)', method_names{2}, constraint_rel(2)))
set(gca, 'XTickLabel', [], 'YTickLabel', [])
if do_save
    save2pdf(sprintf('Simulated_%s_FlexMF_constraint.pdf', data_file), gcf)
end

%% -------- Match to ground truth --------
% helper.similarity_WH_EMD accepts different motif lengths (L vs Lhat):
% it shifts GT profiles into the Lhat window, so no pad/trim is needed.
% eW/eH/ids are 1 x Khat (per estimated factor); ids gives the matched
% ground-truth index (0 = unmatched), so results are placed by ids into
% nMethod x K arrays indexed by ground-truth sequence.
fprintf('\n=== Matching to ground truth ===\n');
emds_W = nan(nMethod, K);
emds_H = nan(nMethod, K);
n_detected = zeros(nMethod, 1);
for m = 1:nMethod
    [eW, eH, ids] = helper.similarity_WH_EMD(W, H, What{m}, Hhat{m});
    matched = ids > 0;
    emds_W(m, ids(matched)) = eW(matched);
    emds_H(m, ids(matched)) = eH(matched);
    n_detected(m) = sum(matched);
    fprintf('  %-32s  n_detected=%d  mean EMD_W=%.4g  mean EMD_H=%.4g\n', ...
        method_names{m}, n_detected(m), mean(eW, 'omitnan'), mean(eH, 'omitnan'));
end

%% -------- Comparison plots --------
figure;
bar(1:K, emds_W');
legend(method_names, 'Location', 'best', 'Interpreter', 'none')
set(gca, 'FontSize', 12)
xlabel('Ground-truth sequence'); ylabel('EMD')
title(sprintf('EMD of W (%s)', data_tag), 'FontSize', 14)
if do_save
    save2pdf(sprintf('Simulated_%s_compare_EMD_W.pdf', data_file), gcf)
end

figure;
bar(1:K, emds_H');
legend(method_names, 'Location', 'best', 'Interpreter', 'none')
set(gca, 'FontSize', 12)
xlabel('Ground-truth sequence'); ylabel('EMD')
title(sprintf('EMD of H (%s)', data_tag), 'FontSize', 14)
if do_save
    save2pdf(sprintf('Simulated_%s_compare_EMD_H.pdf', data_file), gcf)
end

%% -------- Summary --------
fprintf('\n========== Summary (%s) ==========\n', data_tag);
fprintf('%-32s %10s %12s %12s %12s\n', ...
    'Method', 'time(s)', 'EMD_W', 'EMD_H', 'constr_rel');
for m = 1:nMethod
    fprintf('%-32s %10.2f %12.4g %12.4g %12.4g\n', ...
        method_names{m}, times(m), ...
        mean(emds_W(m,:), 'omitnan'), mean(emds_H(m,:), 'omitnan'), constraint_rel(m));
end
fprintf('====================================\n');

%% -------- Save results --------
if do_save
    save('EMD_demo_dynamics_compare.mat', ...
        'data_tag', 'temporal', 'use_burst', 'use_dynamics', 'do_normalize', ...
        'X', 'W', 'H', 'L', 'Lhat', 'K', 'Khat', 'T', ...
        'method_names', 'What', 'Hhat', 'Mhat', 'Rhat', ...
        'times', 'constraint_rel', 'emds_W', 'emds_H', 'n_detected', ...
        'lambda_seqNMF', 'lambda_FlexMF', 'lambda_M', 'lambda_R', ...
        'mu', 'muDecrement', 'maxiter', 'tolerance', 'seed');
    fprintf('Saved EMD_demo_dynamics_compare.mat\n');
end

%% -------- Local helpers --------
function name = sanitize_name(s)
name = regexprep(s, '[^a-zA-Z0-9]+', '_');
end

function plot_WH(W, H, X, name, plotAll)
figure; SimpleWHPlot_patch(W, H, 'plotAll', plotAll);
title(sprintf('%s reconstruction', name), 'Interpreter', 'none')
set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8])
figure; SimpleWHPlot_patch(W, H, 'Data', X, 'plotAll', plotAll, 'compare', true);
title(sprintf('%s factors + data', name), 'Interpreter', 'none')
set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8])
end
