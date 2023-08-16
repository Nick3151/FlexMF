clear all
close all
clc
root = fileparts(pwd);
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))

%% Plot sparsity of W and H
nLambdas = 5;
lambdas = sort([logspace(-1,-3,nLambdas)], 'ascend'); 
for n = 1:nLambdas
    lambda = lambdas(n);
    
    load(sprintf('sparsity_alpha_lambda=%0.3e.mat', lambda))
    
    nWs = length(alpha_Ws);
    nHs = length(alpha_Hs);
    reg_crosses = cell2mat(reg_crosses);
    recon_errors = cell2mat(recon_errors);
    regs = reg_crosses*lambda;
    sparsity_Ws = cell2mat(sparsity_Ws);
    sparsity_Hs = cell2mat(sparsity_Hs);
    scores_W = cellfun(@mean, scores_W);
    scores_H = cellfun(@mean, scores_H);
    
    marginX = .1;
    marginY = .1;
    pos_Ws = [zeros(nWs,1), (nWs-1:-1:0)'*(1-marginY)/nWs,... 
        marginX*ones(nWs,1), (1-marginY)/nWs*ones(nWs,1)];
    pos_Hs = [marginX+(0:nHs-1)'*(1-marginX)/nHs, (1-marginY)*ones(nHs,1),...
        (1-marginX)/nHs*ones(nHs,1), marginY*ones(nHs,1)];
    
    f1 = figure;
    for Wi = 1:nWs
        for Hi = 1:nHs
            ax = subplot('Position', [marginX+(Hi-1)*(1-marginX)/nHs, (nWs-Wi)*(1-marginY)/nWs, (1-marginX)/nHs, (1-marginY)/nWs]);
            boxplot(ax, [sparsity_Ws(:,Wi,Hi),sparsity_Hs(:,Wi,Hi)], 'Colors', 'br')
            set(gca,'YLim',[0,1],'XTick',[],'YTick',[], ...
                'Position', [marginX+(Hi-1)*(1-marginX)/nHs, (nWs-Wi)*(1-marginY)/nWs, (1-marginX)/nHs, (1-marginY)/nWs])
        end
    end
    
    for Wi = 1:nWs
        annotation('textbox',pos_Ws(Wi,:),'String',sprintf('\\alpha_W=%0.3e',alpha_Ws(Wi)),...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
    end
    for Hi = 1:nHs
        annotation('textbox',pos_Hs(Hi,:),'String',sprintf('\\alpha_H=%0.3e',alpha_Hs(Hi)),...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
    end
    set(f1, 'Position',[100 100 900 900])
    export_fig(sprintf('sparsity_alpha_lambda=%0.3e.pdf', lambda))

    % Coefficient scores to ground truth
    f2 = figure;
    for Wi = 1:nWs
        for Hi = 1:nHs
            ax = subplot('Position', [marginX+(Hi-1)*(1-marginX)/nHs, (nWs-Wi)*(1-marginY)/nWs, (1-marginX)/nHs, (1-marginY)/nWs]);
            boxplot(ax, [scores_W(:,Wi,Hi),scores_H(:,Wi,Hi)], 'Colors', 'br')
            set(gca,'YLim',[0,1],'XTick',[],'YTick',[], ...
                'Position', [marginX+(Hi-1)*(1-marginX)/nHs, (nWs-Wi)*(1-marginY)/nWs, (1-marginX)/nHs, (1-marginY)/nWs])
        end
    end
    
    for Wi = 1:nWs
        annotation('textbox',pos_Ws(Wi,:),'String',sprintf('\\alpha_W=%0.3e',alpha_Ws(Wi)),...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
    end
    for Hi = 1:nHs
        annotation('textbox',pos_Hs(Hi,:),'String',sprintf('\\alpha_H=%0.3e',alpha_Hs(Hi)),...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
    end
    set(f2, 'Position',[100 100 900 900])
    save2pdf(sprintf('scores_lambda=%0.3e.pdf', lambda), f2)

    lim= max(max(regs(:)),max(recon_errors(:)));

    f3 = figure;
    for Wi = 1:nWs
        for Hi = 1:nHs
            ax = subplot('Position', [marginX+(Hi-1)*(1-marginX)/nHs, (nWs-Wi)*(1-marginY)/nWs, (1-marginX)/nHs, (1-marginY)/nWs]);
            boxplot(ax, [regs(:,Wi,Hi),recon_errors(:,Wi,Hi)], 'Colors', 'br')
            set(gca,'YLim',[0,lim],'XTick',[],'YTick',[], 'YScale', 'log', ...
                'Position', [marginX+(Hi-1)*(1-marginX)/nHs, (nWs-Wi)*(1-marginY)/nWs, (1-marginX)/nHs, (1-marginY)/nWs])
        end
    end
    for Wi = 1:nWs
        annotation('textbox',pos_Ws(Wi,:),'String',sprintf('\\alpha_W=%0.3e',alpha_Ws(Wi)),...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
    end
    for Hi = 1:nHs
        annotation('textbox',pos_Hs(Hi,:),'String',sprintf('\\alpha_H=%0.3e',alpha_Hs(Hi)),...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
    end
    set(f3, 'Position',[100 100 900 900])
    export_fig(sprintf('errors_alpha_lambda=%0.3e.pdf', lambda))
    
    % Number of detected motifs
    f4 = figure;
    K = 5;
    imagesc(squeeze(mean(cell2mat(num_detected), 1))/K, [0.8,1])
    set(gca, 'Position', [0.1 0.1 0.7 0.7], 'XTickLabel', [], 'YTickLabel', []);
    colorbar('Position', [0.85 0.1 0.05 0.7]);
    
    pos_Ws = [0.1-0.7/nHs*ones(nWs,1), 0.1+(nWs-1:-1:0)'*0.7/nWs,... 
    0.7/nHs*ones(nWs,1), 0.7/nWs*ones(nWs,1)];
    pos_Hs = [0.1+(0:nHs-1)'*0.7/nHs, 0.8*ones(nHs,1),...
        0.7/nHs*ones(nWs,1), 0.7/nWs*ones(nWs,1)];
    for Wi = 1:nWs
        annotation('textbox',pos_Ws(Wi,:),'String',sprintf('\\alpha_W=%0.2e',alpha_Ws(Wi)),...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
    end
    for Hi = 1:nHs
        annotation('textbox',pos_Hs(Hi,:),'String',sprintf('\\alpha_H=%0.2e',alpha_Hs(Hi)),...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
    end
    annotation('textbox',[0.1 0.85 0.8 0.1],'String', 'Detected motifs(%)', 'FontSize', 16, 'FontWeight', 'bold',...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
    set(f4, 'Position',[100 100 900 900])

    % Number of successfully detected motifs
    f5 = figure;
    K = 5;
    imagesc(squeeze(mean(cell2mat(num_success), 1))/K, [0,1])
    set(gca, 'Position', [0.1 0.1 0.7 0.7], 'XTickLabel', [], 'YTickLabel', []);
    colorbar('Position', [0.85 0.1 0.05 0.7]);
    
    for Wi = 1:nWs
        annotation('textbox',pos_Ws(Wi,:),'String',sprintf('\\alpha_W=%0.2e',alpha_Ws(Wi)),...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
    end
    for Hi = 1:nHs
        annotation('textbox',pos_Hs(Hi,:),'String',sprintf('\\alpha_H=%0.2e',alpha_Hs(Hi)),...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
    end
    annotation('textbox',[0.1 0.85 0.8 0.1],'String', 'Successfully detected motifs(%)', 'FontSize', 16, 'FontWeight', 'bold',...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
    set(gcf, 'Position',[100 100 900 900])
    save2pdf(sprintf('nMotifs_success_lambda=%0.3e.pdf', lambda), f5)
end

%% Plot X_hat
n = 3;
lambda = lambdas(n);
load(sprintf('sparsity_alpha_lambda=%0.3e.mat', lambda))
Wi = 5;
Hi = 8;

plotAll = 1;
f6 = figure; SimpleWHPlot_patch(W_hats{1,Wi,Hi},H_hats{1,Wi,Hi},[],[],[],[],plotAll); 
set(f6,'position',[200,200,1200,900])
title('FlexMF reconstruction')
save2pdf(sprintf('Simulated_results_lambda=%0.3e_alphaW=%0.2e_alphaH=%0.2e.pdf', lambda, alpha_Ws(Wi), alpha_Hs(Hi)), f6)

f7 = figure; SimpleWHPlot_patch(Ws{1},Hs{1},[],[],[],[],plotAll); 
set(f7,'position',[200,200,1200,900])
title('Ground Truth')
save2pdf(sprintf('Simulated_results_GT.pdf'), f7)

plotAll = 0;
Wi = 1;
Hi = 8;
What = W_hats{1,Wi,Hi};
Hhat = H_hats{1,Wi,Hi};
N = size(What,1);
% Sort motifs by active neuron ids
[~,kSort] = sort((1:N)*sum(What,3)./sum(What,[1,3]));
f8 = figure; SimpleWHPlot_patch(What(:,kSort,:),Hhat(kSort,:),[],[],[],[],plotAll); 
set(f8,'position',[200,200,1200,900])
title('FlexMF reconstruction')
save2pdf(sprintf('Simulated_results_lambda=%0.3e_alphaW=%0.2e_alphaH=%0.2e.pdf', lambda, alpha_Ws(Wi), alpha_Hs(Hi)), f8)

Wi = 9;
Hi = 1;
What = W_hats{1,Wi,Hi};
Hhat = H_hats{1,Wi,Hi};
N = size(What,1);
% Sort motifs by active neuron ids
[~,kSort] = sort((1:N)*sum(What,3)./sum(What,[1,3]));
f9 = figure; SimpleWHPlot_patch(What(:,kSort,:),Hhat(kSort,:),[],[],[],[],plotAll);
set(f9,'position',[200,200,1200,900])
title('FlexMF reconstruction')
save2pdf(sprintf('Simulated_results_lambda=%0.3e_alphaW=%0.2e_alphaH=%0.2e.pdf', lambda, alpha_Ws(Wi), alpha_Hs(Hi)), f9)