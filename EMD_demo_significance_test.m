%% Demo script: FlexMF with EMD significance test
clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))
rmpath(genpath(fullfile(root, 'seqNMF-master')));
addpath(genpath(fullfile(root, 'FlexMF')));

%% -------- User settings --------
do_choose_lambda = false;  % true: sweep lambdas for SeqNMF before fitting
do_normalize = true;      % scale train/test Frobenius norm to Khat
do_save = true;           % write PDFs

K = 3;
Khat = 5;

%% Generate some sequences with temporal warping
T = 4000; % length of data to generate
Nneurons = 5*ones(K,1); % number of neurons in each sequence
Dt = 3.*ones(K,1); % gap between each member of the sequence
noise = .001; % probability of added noise in each bin
warp = 5; % stretch should be less than Dt
gap = 100;
neg = 0;
seed = 1;
[X, W, H, ~] = generate_data(T,Nneurons,Dt, 'noise',noise, 'warp', warp, 'seed', seed, 'len_burst', 1, 'dynamic', 0);

% Split into training and test set
Xtrain = X(:,1:round(T/2));
Xtest = X(:,1+round(T/2):end);
Htrain = H(:,1:round(T/2));

figure; SimpleWHPlot_patch(W,H,'Data',X,'plotAll', 1, 'onsets', round(T/2));
title('generated data warping','Fontsize',16)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.5])
if do_save
    export_vector_pdf('EMD_demo_significance_generated_data_warp+noise.pdf');
end

%% Normalize train/test
L = size(W, 3);
if do_normalize
    frob_norm_train = norm(Xtrain(:));
    Xtrain = Xtrain / frob_norm_train * Khat;
    W = W / frob_norm_train * Khat;
    frob_norm_test = norm(Xtest(:));
    Xtest = Xtest / frob_norm_test * Khat;
    fprintf('Normalized Xtrain/Xtest to Frobenius norm %g\n', Khat);
end

%% Optional SeqNMF lambda sweep
if do_choose_lambda
    nLambdas = 20; % increase if you're patient
    lambdas = sort(logspace(0,-4,nLambdas), 'ascend');
    loadings = [];
    regularization = [];
    cost = [];
    for li = 1:length(lambdas)
        [What_SeqNMF, Hhat_SeqNMF, ~,~,loadings(li,:),power]= seqNMF(Xtrain,'K',Khat,'L',L,...
        'lambda', lambdas(li), 'maxiter', 100, 'showPlot', 0);
        [cost(li),regularization(li),~] = helper.get_seqNMF_cost(Xtrain,What_SeqNMF,Hhat_SeqNMF);
        display(['Testing lambda ' num2str(li) '/' num2str(length(lambdas))])
    end

    % plot costs as a function of lambda
    windowSize = 3;
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    Rs = filtfilt(b,a,regularization);
    minRs = prctile(regularization,10); maxRs= prctile(regularization,90);
    Rs = (Rs-minRs)/(maxRs-minRs);
    R = (regularization-minRs)/(maxRs-minRs);
    Cs = filtfilt(b,a,cost);
    minCs =  prctile(cost,10); maxCs =  prctile(cost,90);
    Cs = (Cs -minCs)/(maxCs-minCs);
    C = (cost -minCs)/(maxCs-minCs);

    figure; hold on
    plot(lambdas,Rs, 'b')
    plot(lambdas,Cs,'r')
    scatter(lambdas, R, 'b', 'markerfacecolor', 'flat');
    scatter(lambdas, C, 'r', 'markerfacecolor', 'flat');
    xlabel('Lambda'); ylabel('Cost (au)')
    set(legend('Correlation cost', 'Reconstruction cost'), 'Box', 'on')
    set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
    set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
    if do_save
        export_vector_pdf('EMD_demo_significance_choose_lambda_SeqNMF.pdf');
    end
end

%% Run SeqNMF on training data
lambda_SeqNMF = .1;
lambdaL1H = 0;
lambdaL1W = 0;
lambdaOrthoH = 0;

figure;
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
[What_SeqNMF, Hhat_SeqNMF, ~, errors_SeqNMF,loadings,power]= seqNMF(Xtrain,'K',Khat,'L',L,...
            'lambda', lambda_SeqNMF, 'maxiter', 50, 'showPlot', 1);

%% Match SeqNMF motifs to ground truth
fprintf('\n=== Matching SeqNMF factors to ground truth ===\n');
[emds_W_SeqNMF, emds_H_SeqNMF, ids_SeqNMF] = helper.similarity_WH_EMD(W, Htrain, What_SeqNMF, Hhat_SeqNMF);
[What_SeqNMF_plot, Hhat_SeqNMF_plot] = helper.sort_matched_factors(What_SeqNMF, Hhat_SeqNMF, ids_SeqNMF);
num_matched_SeqNMF = sum(ids_SeqNMF > 0);
fprintf('Matched %d / %d SeqNMF factors (K=%d ground truth)\n', num_matched_SeqNMF, Khat, K);
disp('ids_SeqNMF (estimate -> GT index, 0 = unmatched):');
disp(ids_SeqNMF);
disp('emds_W:'); disp(emds_W_SeqNMF);
disp('emds_H:'); disp(emds_H_SeqNMF);

plotAll = 1;
figure; SimpleWHPlot_patch(What_SeqNMF_plot, Hhat_SeqNMF_plot, 'plotAll', plotAll);
title('SeqNMF reconstruction (ground-truth order)')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
if do_save
    export_vector_pdf('EMD_demo_significance_SeqNMF_recon.pdf');
end

figure; SimpleWHPlot_patch(What_SeqNMF_plot, Hhat_SeqNMF_plot, 'Data', Xtrain, 'plotAll', plotAll);
title('SeqNMF factors, with raw data (ground-truth order)')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
if do_save
    export_vector_pdf('EMD_demo_significance_SeqNMF_raw.pdf');
end

%% Find sequence with FlexMF on training data (SeqNMF warm-start)
lambda_FlexMF = .05;
lambda_M = .05;
lambda_R = 1;
fprintf('\n=== FlexMF train (lambda=%g, lambda_M=%g, lambda_R=%g, SeqNMF init) ===\n', ...
    lambda_FlexMF, lambda_M, lambda_R);
figure;
t_train = tic;
[What, Hhat_train, cost_train, errors_train, loadings, power, M_train, R_train] = FlexMF(Xtrain, 'K', Khat, 'L', L, ...
    'EMD',1, 'lambda', lambda_FlexMF, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'maxiter', 50, ...
    'W_init', What_SeqNMF, 'H_init', Hhat_SeqNMF);
time_train = toc(t_train);
fprintf('FlexMF training time: %.2f s\n', time_train);

%% Match FlexMF motifs to ground truth
fprintf('\n=== Matching FlexMF factors to ground truth ===\n');
[emds_W_FlexMF, emds_H_FlexMF, ids_FlexMF] = helper.similarity_WH_EMD(W, Htrain, What, Hhat_train);
[What_plot, Hhat_train_plot] = helper.sort_matched_factors(What, Hhat_train, ids_FlexMF);
num_matched = sum(ids_FlexMF > 0);
fprintf('Matched %d / %d FlexMF factors (K=%d ground truth)\n', num_matched, Khat, K);
disp('ids_FlexMF (estimate -> GT index, 0 = unmatched):');
disp(ids_FlexMF);
disp('emds_W:'); disp(emds_W_FlexMF);
disp('emds_H:'); disp(emds_H_FlexMF);

%% Plot M, R
figure; plot_MR(M_train, R_train)
if do_save
    export_vector_pdf('EMD_demo_significance_train_MR.pdf');
end

%% Look at factors (ground-truth order)
figure; SimpleWHPlot_patch(What_plot, Hhat_train_plot, 'plotAll', plotAll);
title('FlexMF reconstruction (ground-truth order)')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
if do_save
    export_vector_pdf('EMD_demo_significance_train_recon.pdf');
end

figure; SimpleWHPlot_patch(What_plot, Hhat_train_plot, 'Data', Xtrain, 'plotAll', plotAll);
title('FlexMF factors, with raw data (ground-truth order)')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
if do_save
    export_vector_pdf('EMD_demo_significance_train_raw.pdf');
end

%% Fix What, rerun FlexMF on test data
fprintf('\n=== FlexMF test (W fixed) ===\n');
figure;
t_test = tic;
[What, Hhat_test, cost_test, errors_test, ~, ~, M_test, R_test] = FlexMF(Xtest, 'K', Khat, 'L', L, 'W_fixed', 1, 'W_init', What,...
    'EMD',1, 'lambda', lambda_FlexMF, 'lambda_R', lambda_R, 'lambda_M', lambda_M, 'maxiter', 50);
time_test = toc(t_test);
fprintf('FlexMF test time: %.2f s\n', time_test);
[What_plot, Hhat_test_plot] = helper.sort_matched_factors(What, Hhat_test, ids_FlexMF);

figure; SimpleWHPlot_patch(What_plot, Hhat_test_plot, 'plotAll', 1); title('FlexMF test recon (ground-truth order)')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
if do_save
    export_vector_pdf('EMD_demo_significance_test_recon.pdf');
end

figure; SimpleWHPlot_patch(What_plot, Hhat_test_plot, 'Data', Xtest, 'plotAll', 1); title('FlexMF test raw (ground-truth order)')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
if do_save
    export_vector_pdf('EMD_demo_significance_test_raw.pdf');
end

figure; plot_MR(M_test, R_test)
if do_save
    export_vector_pdf('EMD_demo_significance_test_MR.pdf');
end

%% test significance
fprintf('\n=== Significance test ===\n');
[pvals,is_significant,is_single,figData] = test_significance_EMD(Xtest, What, M_test, 'plot', 1);
if do_save && ~isempty(figData) && isgraphics(figData)
    export_vector_pdf('EMD_demo_significance_test_data.pdf', figData);
end

fprintf('\n========== Summary ==========\n');
fprintf('  FlexMF train time: %.2f s\n', time_train);
fprintf('  FlexMF test time:  %.2f s\n', time_test);
fprintf('  SeqNMF matched:    %d / %d (K=%d GT)\n', num_matched_SeqNMF, Khat, K);
fprintf('  FlexMF matched:    %d / %d (K=%d GT)\n', num_matched, Khat, K);
fprintf('  Significant factors: %d / %d\n', sum(is_significant), Khat);
fprintf('  p-values:      '); fprintf('%g ', pvals); fprintf('\n');
fprintf('  is_significant:'); fprintf('%d ', is_significant); fprintf('\n');
fprintf('  is_single:     '); fprintf('%d ', is_single); fprintf('\n');
fprintf('==============================\n');

%% Save results
if do_save
    results_file = 'EMD_demo_significance_results.mat';
    save(results_file, 'X', 'W', 'H', 'Htrain', 'Xtrain', 'Xtest', 'L', ...
        'What_SeqNMF', 'Hhat_SeqNMF', 'What_SeqNMF_plot', 'Hhat_SeqNMF_plot', ...
        'ids_SeqNMF', 'emds_W_SeqNMF', 'emds_H_SeqNMF', ...
        'What', 'Hhat_train', 'Hhat_test', 'What_plot', ...
        'Hhat_train_plot', 'Hhat_test_plot', ...
        'ids_FlexMF', 'emds_W_FlexMF', 'emds_H_FlexMF', ...
        'M_train', 'R_train', 'M_test', 'R_test', ...
        'pvals', 'is_significant', 'is_single', ...
        'time_train', 'time_test', 'K', 'Khat', ...
        'lambda_SeqNMF', 'lambda_FlexMF', 'lambda_M', 'lambda_R');
    fprintf('Saved results: %s\n', results_file);
end

