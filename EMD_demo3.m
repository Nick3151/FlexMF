%% Demo: Extreme warp/jitter+noise where SeqNMF splits one sequence into two
% Large warp or jitter so SeqNMF fails (one ground-truth sequence recovered
% as multiple factors). FlexMF is warm-started from SeqNMF and compared.
%
% Supported data_type values: 'jitter+noise' | 'warp+noise'
clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))
rmpath(genpath(fullfile(root, 'seqNMF-master')));
addpath(genpath(fullfile(root, 'FlexMF')));

%% -------- User settings --------
data_type = 'warp+noise';   % 'jitter+noise' | 'warp+noise'
do_choose_lambda = true;  % true: sweep lambdas for SeqNMF before fitting
do_normalize = true;       % scale data Frobenius norm to Khat
do_save = true;            % write PDFs

K = 3;                     % number of ground-truth sequences
Khat = 5;                  % number of factors to fit (overcomplete so SeqNMF can split)
T = 4000;
Nneurons = 10*ones(K,1);
Dt = 3.*ones(K,1);
noise = 0.005;             % probability of added noise in each bin
jitter = 8*ones(K,1);      % large jitter so SeqNMF fails
participation = 1.*ones(K,1);
warp = 10;                  % large warp so SeqNMF fails
seed = 1;
maxiter = 50;

% Regularization
lambda_SeqNMF = .05;
lambda_FlexMF = .05;
lambda_M = .05;
lambda_R = 1;
lambdaL1H = 0;
tolerance = 1e-3;
reweight = 0;

%% -------- Generate data --------
base_args = {'seed', seed, 'len_burst', 1, 'dynamic', 0};
switch data_type
    case 'jitter+noise'
        gen_args = [{'noise', noise, 'jitter', jitter}, base_args];
    case 'warp+noise'
        gen_args = [{'noise', noise, 'warp', warp}, base_args];
    otherwise
        error('Unknown data_type ''%s''. Choose jitter+noise or warp+noise.', data_type);
end

[X, W, H, ~] = generate_data(T, Nneurons, Dt, gen_args{:});
L = size(W, 3);
fprintf('Data type: %s  |  size(X)=[%d %d]  L=%d  K=%d  Khat=%d\n', ...
    data_type, size(X,1), size(X,2), L, K, Khat);

plotAll = 1;
figure; SimpleWHPlot_patch(W, H, 'Data', X, 'plotAll', plotAll);
title(sprintf('Generated data (%s, extreme)', data_type), 'FontSize', 16)
set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8])
if do_save
    export_vector_pdf(sprintf('EMD_simulated_data_%s_extreme.pdf', helper.sanitize_name(data_type)));
end

%% -------- Normalize --------
if do_normalize
    frob_norm = norm(X(:));
    X = X / frob_norm * Khat;
    W = W / frob_norm * Khat;
end

%% -------- Optional SeqNMF lambda sweep --------
if do_choose_lambda
    nLambdas = 20;
    lambdas = sort(logspace(0, -4, nLambdas), 'ascend');
    regularization = zeros(1, nLambdas);
    cost = zeros(1, nLambdas);
    for li = 1:nLambdas
        [What_tmp, Hhat_tmp] = seqNMF(X, 'K', Khat, 'L', L, ...
            'lambda', lambdas(li), 'maxiter', maxiter, 'showPlot', 0);
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
    title(sprintf('SeqNMF lambda sweep (%s, extreme)', data_type))
    if do_save
        export_vector_pdf(sprintf('Simulate_%s_extreme_choose_lambda_SeqNMF', helper.sanitize_name(data_type)));
    end
end

%% -------- Run SeqNMF --------
fprintf('\n=== SeqNMF (lambda=%g, Khat=%d) ===\n', lambda_SeqNMF, Khat);
t_seq = tic;
figure;
[What_SeqNMF, Hhat_SeqNMF, ~, errors_SeqNMF, loadings_SeqNMF, power_SeqNMF] = ...
    seqNMF(X, 'K', Khat, 'L', L, 'lambda', lambda_SeqNMF, 'maxiter', maxiter, 'showPlot', 1);
time_SeqNMF = toc(t_seq);
fprintf('SeqNMF running time: %.2f s\n', time_SeqNMF);

%% -------- Run FlexMF (SeqNMF init) --------
fprintf('\n=== FlexMF (lambda=%g, lambda_M=%g, lambda_R=%g, SeqNMF init) ===\n', ...
    lambda_FlexMF, lambda_M, lambda_R);
t_flex = tic;
figure;
[What_FlexMF, Hhat_FlexMF, cost_FlexMF, errors_FlexMF, loadings_FlexMF, power_FlexMF, M, R] = ...
    FlexMF(X, 'K', Khat, 'L', L, ...
    'EMD', 1, 'lambda', lambda_FlexMF, ...
    'lambdaL1H', lambdaL1H, 'lambda_R', lambda_R, 'lambda_M', lambda_M, ...
    'maxiter', maxiter, 'tolerance', tolerance, ...
    'neg_prop', 0, 'Reweight', reweight, ...
    'W_init', What_SeqNMF, 'H_init', Hhat_SeqNMF);
time_FlexMF = toc(t_flex);
fprintf('FlexMF running time: %.2f s\n', time_FlexMF);

figure;
plot_MR(M, R)
if do_save
    export_vector_pdf(sprintf('FlexMF_%s_extreme_demo_MR_lambda=%1.1e_lambdaM=%1.1e_lambdaR=%1.1e.pdf', ...
        helper.sanitize_name(data_type), lambda_FlexMF, lambda_M, lambda_R));
end

%% -------- Match factors to ground truth --------
fprintf('\n=== Matching factors to ground truth ===\n');
t_match = tic;
[emds_W_SeqNMF, emds_H_SeqNMF, ids_SeqNMF] = helper.similarity_WH_EMD(W, H, What_SeqNMF, Hhat_SeqNMF);
[emds_W_FlexMF, emds_H_FlexMF, ids_FlexMF] = helper.similarity_WH_EMD(W, H, What_FlexMF, Hhat_FlexMF);
[What_SeqNMF_plot, Hhat_SeqNMF_plot] = helper.sort_matched_factors(What_SeqNMF, Hhat_SeqNMF, ids_SeqNMF);
[What_FlexMF_plot, Hhat_FlexMF_plot] = helper.sort_matched_factors(What_FlexMF, Hhat_FlexMF, ids_FlexMF);
[coeffs_W_SeqNMF, coeffs_H_SeqNMF, ~] = helper.similarity_WH(W, H, What_SeqNMF, Hhat_SeqNMF);
[coeffs_W_FlexMF, coeffs_H_FlexMF, ~] = helper.similarity_WH(W, H, What_FlexMF, Hhat_FlexMF);
time_match = toc(t_match);
fprintf('Matching running time: %.2f s\n', time_match);

figure; SimpleWHPlot_patch(What_SeqNMF_plot, Hhat_SeqNMF_plot, 'plotAll', plotAll);
title('SeqNMF reconstruction (ground-truth order)')
set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8])
if do_save
    export_vector_pdf(sprintf('Simulated_%s_extreme_result_SeqNMF.pdf', helper.sanitize_name(data_type)), gcf);
end

figure; SimpleWHPlot_patch(What_SeqNMF_plot, Hhat_SeqNMF_plot, 'Data', X, 'plotAll', plotAll);
title('SeqNMF factors, with raw data (ground-truth order)')
set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8])

figure; SimpleWHPlot_patch(What_FlexMF_plot, Hhat_FlexMF_plot, 'plotAll', plotAll);
title('FlexMF reconstruction (ground-truth order)')
set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8])
if do_save
    export_vector_pdf(sprintf('EMD_Simulated_%s_extreme_data_FlexMF_lambda=%1.1e_lambdaM=%1.1e_lambdaR=%1.1e.pdf', ...
        helper.sanitize_name(data_type), lambda_FlexMF, lambda_M, lambda_R), gcf);
end

figure; SimpleWHPlot_patch(What_FlexMF_plot, Hhat_FlexMF_plot, 'Data', X, 'plotAll', plotAll);
title('FlexMF factors, with raw data (ground-truth order)')
set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8])

% Match against K ground-truth factors (Khat may be larger)
emds_W_all = zeros(2, K);
matched_SeqNMF = ids_SeqNMF > 0 & ids_SeqNMF <= K;
matched_FlexMF = ids_FlexMF > 0 & ids_FlexMF <= K;
emds_W_all(1, ids_SeqNMF(matched_SeqNMF)) = emds_W_SeqNMF(matched_SeqNMF);
emds_W_all(2, ids_FlexMF(matched_FlexMF)) = emds_W_FlexMF(matched_FlexMF);

emds_H_all = zeros(2, K);
emds_H_all(1, ids_SeqNMF(matched_SeqNMF)) = emds_H_SeqNMF(matched_SeqNMF);
emds_H_all(2, ids_FlexMF(matched_FlexMF)) = emds_H_FlexMF(matched_FlexMF);

figure; bar(1:K, emds_W_all');
legend({'SeqNMF', 'FlexMF'}, 'Location', 'north')
set(gca, 'FontSize', 14)
title(sprintf('EMDs of W (%s, extreme)', data_type), 'FontSize', 16)
if do_save
    export_vector_pdf(sprintf('EMD_Simulated_%s_extreme_data_compare_W.pdf', helper.sanitize_name(data_type)), gcf);
end

figure; bar(1:K, emds_H_all');
legend({'SeqNMF', 'FlexMF'}, 'Location', 'north')
set(gca, 'FontSize', 14)
title(sprintf('EMDs of H (%s, extreme)', data_type), 'FontSize', 16)
if do_save
    export_vector_pdf(sprintf('EMD_Simulated_%s_extreme_data_compare_H.pdf', helper.sanitize_name(data_type)), gcf);
end

%% -------- Save results --------
if do_save
    results_file = sprintf('EMD_Simulated_%s_extreme_results.mat', helper.sanitize_name(data_type));
    settings = struct('data_type', data_type, 'K', K, 'Khat', Khat, 'T', T, ...
        'Nneurons', Nneurons, 'Dt', Dt, 'noise', noise, 'jitter', jitter, ...
        'participation', participation, 'warp', warp, 'seed', seed, ...
        'maxiter', maxiter, 'lambda_SeqNMF', lambda_SeqNMF, ...
        'lambda_FlexMF', lambda_FlexMF, 'lambda_M', lambda_M, ...
        'lambda_R', lambda_R, 'lambdaL1H', lambdaL1H, ...
        'tolerance', tolerance, 'reweight', reweight, ...
        'do_normalize', do_normalize);
    save(results_file, 'X', 'W', 'H', 'L', 'What_SeqNMF', 'Hhat_SeqNMF', ...
        'What_FlexMF', 'Hhat_FlexMF', 'What_SeqNMF_plot', 'Hhat_SeqNMF_plot', ...
        'What_FlexMF_plot', 'Hhat_FlexMF_plot', 'M', 'R', ...
        'ids_SeqNMF', 'ids_FlexMF', 'emds_W_SeqNMF', 'emds_H_SeqNMF', ...
        'emds_W_FlexMF', 'emds_H_FlexMF', 'coeffs_W_SeqNMF', ...
        'coeffs_H_SeqNMF', 'coeffs_W_FlexMF', 'coeffs_H_FlexMF', ...
        'emds_W_all', 'emds_H_all', 'time_SeqNMF', 'time_FlexMF', ...
        'time_match', 'settings');
    fprintf('Saved results: %s\n', results_file);
end

%% -------- Timing summary --------
fprintf('\n========== Timing summary (%s, extreme) ==========\n', data_type);
fprintf('  SeqNMF:   %8.2f s\n', time_SeqNMF);
fprintf('  FlexMF:   %8.2f s\n', time_FlexMF);
fprintf('  Matching: %8.2f s\n', time_match);
fprintf('==========================================\n');
