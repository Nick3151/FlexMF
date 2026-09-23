%% Plot results from simulation_choose_params_EMD_test
% Loads Simulation_Results/FlexMF_choose_lambda_lambdaM_*.mat and visualizes
% the lambda x lambda_M sweep (factors, M/R, EMD / detection heatmaps).
% Together with simulation_choose_lambdas_SeqNMF.m, this replaces the
% exploratory workflow formerly in EMD_test_params_2d.m.
clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))
rmpath(genpath(fullfile(root, 'seqNMF-master')));
addpath(genpath(fullfile(root, 'FlexMF')));

%% Configuration (must match a saved sweep)
data_type = 'warp+noise';   % 'jitter+noise' | 'warp+noise'
lambda_R = 1;
plotAll = 1;

out_dir = 'Simulation_EMD';
results_file = fullfile(out_dir, sprintf( ...
    'FlexMF_choose_lambda_lambdaM_%s_lambdaR=%0.3e.mat', ...
    helper.sanitize_name(data_type), lambda_R));
assert(exist(results_file, 'file') == 2, 'Missing results: %s', results_file);
S = load(results_file);
fprintf('Loaded %s\n', results_file);

lambdas = S.lambdas;
lambda_Ms = S.lambda_Ms;
nLambdas = numel(lambdas);
nMs = numel(lambda_Ms);
nSim = size(S.W_hats, 1);
assert(isfield(S, 'L1_Xtrain'), ...
    'L1_Xtrain missing; re-run simulation_choose_params_EMD_test.m');
L1_Xtrain = S.L1_Xtrain(:);

%% Aggregate grids over lambda x lambda_M (average across sims if nSim > 1)
mean_emd_W = nan(nLambdas, nMs);
mean_emd_H = nan(nLambdas, nMs);
num_detected_all = nan(nLambdas, nMs);
num_significant_all = nan(nLambdas, nMs);
M_norms = nan(nLambdas, nMs);
R_norms = nan(nLambdas, nMs);
constraints_rel = nan(nLambdas, nMs);
reg_costs = nan(nLambdas, nMs);
time_train = nan(nLambdas, nMs);
time_emd = nan(nLambdas, nMs);
time_test = nan(nLambdas, nMs);

for Li = 1:nLambdas
    for Mi = 1:nMs
        emdW = nan(nSim, 1);
        emdH = nan(nSim, 1);
        nd = nan(nSim, 1);
        ns = nan(nSim, 1);
        m1 = nan(nSim, 1);
        r1 = nan(nSim, 1);
        cons_rel = nan(nSim, 1);
        reg = nan(nSim, 1);
        tt = nan(nSim, 1);
        te = nan(nSim, 1);
        ts = nan(nSim, 1);
        for n = 1:nSim
            emdW(n) = mean(S.emds_W{n, Li, Mi}, 'omitnan');
            emdH(n) = mean(S.emds_H{n, Li, Mi}, 'omitnan');
            nd(n) = S.num_detected{n, Li, Mi};
            ns(n) = S.num_significant{n, Li, Mi};
            m1(n) = norm(S.Ms{n, Li, Mi}(:), 1) / L1_Xtrain(n);
            r1(n) = norm(S.Rs{n, Li, Mi}(:), 1) / L1_Xtrain(n);
            if isfield(S, 'constraints_rel')
                cons_rel(n) = S.constraints_rel(n, Li, Mi);
            end
            if isfield(S, 'reg_costs')
                reg(n) = S.reg_costs(n, Li, Mi);
            end
            if isfield(S, 'times')
                tt(n) = S.times{n, Li, Mi};
            end
            if isfield(S, 'times_emd')
                te(n) = S.times_emd{n, Li, Mi};
            end
            if isfield(S, 'times_test')
                ts(n) = S.times_test{n, Li, Mi};
            end
        end
        mean_emd_W(Li, Mi) = mean(emdW, 'omitnan');
        mean_emd_H(Li, Mi) = mean(emdH, 'omitnan');
        num_detected_all(Li, Mi) = mean(nd);
        num_significant_all(Li, Mi) = mean(ns);
        M_norms(Li, Mi) = mean(m1);
        R_norms(Li, Mi) = mean(r1);
        constraints_rel(Li, Mi) = mean(cons_rel, 'omitnan');
        reg_costs(Li, Mi) = mean(reg, 'omitnan');
        time_train(Li, Mi) = mean(tt);
        time_emd(Li, Mi) = mean(te);
        time_test(Li, Mi) = mean(ts);
    end
end

tick_lambda = arrayfun(@(x) sprintf('%.0e', x), lambdas, 'UniformOutput', false);
tick_lambda_M = arrayfun(@(x) sprintf('%.0e', x), lambda_Ms, 'UniformOutput', false);

plot_param_heatmap(mean_emd_W, lambdas, lambda_Ms, tick_lambda, tick_lambda_M, ...
    sprintf('mean EMD(W)  [%s]', data_type));
plot_param_heatmap(mean_emd_H, lambdas, lambda_Ms, tick_lambda, tick_lambda_M, ...
    sprintf('mean EMD(H)  [%s]', data_type));
plot_param_heatmap(num_significant_all, lambdas, lambda_Ms, tick_lambda, tick_lambda_M, ...
    sprintf('num significant  [%s]', data_type));
plot_param_heatmap(M_norms, lambdas, lambda_Ms, tick_lambda, tick_lambda_M, ...
    sprintf('||M||_1 / ||X||_1  [%s]', data_type));
plot_param_heatmap(R_norms, lambdas, lambda_Ms, tick_lambda, tick_lambda_M, ...
    sprintf('||R||_1 / ||X||_1  [%s]', data_type), [0 1]);
plot_param_heatmap(constraints_rel, lambdas, lambda_Ms, tick_lambda, tick_lambda_M, ...
    sprintf('||constraint||_1 / ||X||_1  [%s]', data_type), [0 1]);
plot_param_heatmap(reg_costs, lambdas, lambda_Ms, tick_lambda, tick_lambda_M, ...
    sprintf('regularization cost  [%s]', data_type));
    
plot_param_heatmap(time_train, lambdas, lambda_Ms, tick_lambda, tick_lambda_M, ...
    sprintf('train FlexMF time (s)  [%s]', data_type));

% plot_param_heatmap(time_emd, lambdas, lambda_Ms, tick_lambda, tick_lambda_M, ...
%     sprintf('similarity_WH_EMD time (s)  [%s]', data_type));

% plot_param_heatmap(time_test, lambdas, lambda_Ms, tick_lambda, tick_lambda_M, ...
%     sprintf('test FlexMF time (s)  [%s]', data_type));

%% Inspect one simulation
sim_idx = 1;               % which simulation replicate to inspect
Li_show = 4;               % lambda index for factor / M/R plots
Mi_show = 3;               % lambda_M index for factor / M/R plots
assert(sim_idx >= 1 && sim_idx <= nSim, 'sim_idx out of range (nSim=%d).', nSim);
assert(Li_show >= 1 && Li_show <= nLambdas, 'Li_show out of range.');
assert(Mi_show >= 1 && Mi_show <= nMs, 'Mi_show out of range.');

W = S.Ws{sim_idx};
H = S.Hs{sim_idx};
X = S.Xs{sim_idx};
T = size(X, 2);
tSplit = round(T / 2);
Htrain = H(:, 1:tSplit);
Xtrain = X(:, 1:tSplit);

W_hat = S.W_hats{sim_idx, Li_show, Mi_show};
H_hat = S.H_hats{sim_idx, Li_show, Mi_show};
M = S.Ms{sim_idx, Li_show, Mi_show};
R = S.Rs{sim_idx, Li_show, Mi_show};
ids = S.ids_match{sim_idx, Li_show, Mi_show};

figure;
SimpleWHPlot_patch(W, Htrain, 'Data', Xtrain, 'plotAll', plotAll);
title(sprintf('Generated data (%s), sim %d', data_type, sim_idx), 'FontSize', 16);
set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);

if isfield(S, 'What_SeqNMFs') && ~isempty(S.What_SeqNMFs{sim_idx})
    assert(isfield(S, 'ids_match_SeqNMF'), ...
        'ids_match_SeqNMF missing; re-run simulation_choose_params_EMD_test.m');
    [What_Seq_plot, Hhat_Seq_plot] = helper.sort_matched_factors( ...
        S.What_SeqNMFs{sim_idx}, S.Hhat_SeqNMFs{sim_idx}, ...
        S.ids_match_SeqNMF{sim_idx});
    figure;
    SimpleWHPlot_patch(What_Seq_plot, Hhat_Seq_plot, 'plotAll', plotAll);
    title(sprintf('SeqNMF warm-start (GT order), sim %d', sim_idx), 'FontSize', 16);
    set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
end

[What_plot, Hhat_plot] = helper.sort_matched_factors(W_hat, H_hat, ids);
figure;
SimpleWHPlot_patch(What_plot, Hhat_plot, 'plotAll', plotAll);
title(sprintf('FlexMF  \\lambda=%g, \\lambda_M=%g, \\lambda_R=%g (GT order), sim %d', ...
    lambdas(Li_show), lambda_Ms(Mi_show), lambda_R, sim_idx), 'FontSize', 16);
set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);

figure;
plot_MR(M, R);

fprintf('Inspected sim %d, Li=%d (lambda=%g), Mi=%d (lambda_M=%g).\n', ...
    sim_idx, Li_show, lambdas(Li_show), Mi_show, lambda_Ms(Mi_show));

function plot_param_heatmap(Z, lambdas, lambda_Ms, tick_lambda, tick_lambda_M, fig_title, clims)
figure;
if nargin < 7 || isempty(clims)
    imagesc(Z);
else
    imagesc(Z, clims);
end
colorbar
colormap(flipud(gray(256)));
title(fig_title, 'FontSize', 14, 'Interpreter', 'none');
xlabel('\lambda_M');
ylabel('\lambda');
xticks(1:numel(lambda_Ms));
xticklabels(tick_lambda_M);
yticks(1:numel(lambdas));
yticklabels(tick_lambda);
set(gcf, 'Units', 'normalized', 'Position', [0.15 0.15 0.55 0.6]);
end
