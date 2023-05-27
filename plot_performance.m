clear all
close all
clc
%% plot costs as a function of lambda and alpha
load choose_lambda_clean.mat;
nLambdas = length(lambdas); 
nAlphas = length(alphas);
nSim = 100;

reg_crosses = cell2mat(reg_crosses);
recon_errors = cell2mat(recon_errors);

regs = reg_crosses.*repmat(lambdas, nSim, 1, nAlphas);
maxregs= max(regs(:));

f1 = figure;
marginX = 0.1;
marginY = 0.1;
for li = 1:nLambdas
    for ai = 1:nAlphas
        ax = subplot('Position', [marginX+(ai-1)*(1-marginX)/nAlphas, (nLambdas-li)*(1-marginY)/nLambdas, (1-marginX)/nAlphas*.8, (1-marginY)/nLambdas*.8]);
%         b = bar([regs(li,ai), recon_errors(li,ai)],'FaceColor','flat');
%         b.CData(2,:) = [1 0 0];\
        boxplot(ax, [regs(:,li,ai),recon_errors(:,li,ai)], 'Colors', 'br')
        set(gca,'YLim',[0,maxregs],'XTick',[],'YTick',[], ...
            'Position', [marginX+(ai-1)*(1-marginX)/nAlphas, (nLambdas-li)*(1-marginY)/nLambdas, (1-marginX)/nAlphas, (1-marginY)/nLambdas])
    end
end
pos_Lambdas = [zeros(nLambdas,1), (nLambdas-1:-1:0)'*(1-marginY)/nLambdas,... 
    (1-marginX)/nAlphas*ones(nLambdas,1), (1-marginY)/nLambdas*ones(nLambdas,1)];
pos_Alphas = [marginX+(0:nAlphas-1)'*(1-marginX)/nAlphas, (1-marginY)*ones(nAlphas,1),...
    (1-marginX)/nAlphas*ones(nAlphas,1), (1-marginY)/nLambdas*ones(nAlphas,1)];
for li = 1:nLambdas
    annotation('textbox',pos_Lambdas(li,:),'String',sprintf('\\lambda=%0.3e',lambdas(li)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
for ai = 1:nAlphas
    annotation('textbox',pos_Alphas(ai,:),'String',sprintf('\\alpha=%0.3e',alphas(ai)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end  
set(f1, 'Units', 'normalized', 'Position',[.1 .1 .7 .7])
save2pdf('choose_lambda_clean.pdf', f1)

% % Plot costs for single alpha
% windowSize = 3; 
% b = (1/windowSize)*ones(1,windowSize);
% a = 1;
% ai = 3;
% alpha_sel = alphas(ai);
% 
% r = filtfilt(b,a,regs(:,ai)); 
% c = filtfilt(b,a,recon_errors(:,ai)); 
% 
% figure; clf; hold on
% plot(lambdas,r, 'b')
% plot(lambdas,c,'r')
% scatter(lambdas, regs(:,ai), 'b', 'markerfacecolor', 'flat');
% scatter(lambdas, recon_errors(:,ai), 'r', 'markerfacecolor', 'flat');
% xlabel('Lambda'); ylabel('Cost')
% set(legend('Correlation cost', 'Reconstruction cost'), 'Box', 'on')
% set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
% set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
% title(['alpha=', num2str(alpha_sel)])
% 
% % Plot costs for single lambda
% windowSize = 3; 
% b = (1/windowSize)*ones(1,windowSize);
% a = 1;
% li = 7;
% lambda_sel = lambdas(li);
% 
% r = filtfilt(b,a,regs(li,:)); 
% c = filtfilt(b,a,recon_errors(li,:)); 
% 
% figure; clf; hold on
% plot(alphas,r, 'b')
% plot(alphas,c,'r')
% scatter(alphas, regs(li,:), 'b', 'markerfacecolor', 'flat');
% scatter(alphas, recon_errors(li,:), 'r', 'markerfacecolor', 'flat');
% xlabel('Alpha'); ylabel('Cost')
% set(legend('Correlation cost', 'Reconstruction cost'), 'Box', 'on')
% set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
% set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
% title(['lambda=', num2str(lambda_sel)])

%% Regularize costs and remake plots
% minRsel = prctile(reg_crosses(li,:),10); maxRsel= prctile(reg_crosses(li,:),90);
% RegSel = (reg_crosses(li,:)-minRsel)/(maxRsel-minRsel); 
% minCsel =  prctile(recon_errors(li,:),10); maxCsel =  prctile(recon_errors(li,:),90); 
% CostSel = (recon_errors(li,:) -minCsel)/(maxCsel-minCsel); 
% R = filtfilt(b,a,RegSel); 
% C = filtfilt(b,a,CostSel); 
% 
% figure; clf; hold on
% plot(alphas,R, 'b')
% plot(alphas,C,'r')
% scatter(alphas, RegSel, 'b', 'markerfacecolor', 'flat');
% scatter(alphas, CostSel, 'r', 'markerfacecolor', 'flat');
% xlabel('Alpha'); ylabel('Cost (au)')
% set(legend('Regularized Correlation cost', 'Regularized Reconstruction cost'), 'Box', 'on')
% set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
% set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
% title(['lambda=', num2str(lambda_sel)])
% 
% minRsel = prctile(reg_crosses(:,ai),10); maxRsel= prctile(reg_crosses(:,ai),90);
% RegSel = (reg_crosses(:,ai)-minRsel)/(maxRsel-minRsel); 
% minCsel =  prctile(recon_errors(:,ai),10); maxCsel =  prctile(recon_errors(:,ai),90); 
% CostSel = (recon_errors(:,ai) -minCsel)/(maxCsel-minCsel); 
% R = filtfilt(b,a,RegSel); 
% C = filtfilt(b,a,CostSel); 
% 
% figure; clf; hold on
% plot(lambdas,R, 'b')
% plot(lambdas,C,'r')
% scatter(lambdas, RegSel, 'b', 'markerfacecolor', 'flat');
% scatter(lambdas, CostSel, 'r', 'markerfacecolor', 'flat');
% xlabel('Lambda'); ylabel('Cost (au)')
% set(legend('Regularized Correlation cost', 'Regularized Reconstruction cost'), 'Box', 'on')
% set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
% set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
% title(['alpha=', num2str(alpha_sel)])

minRs = min(reg_crosses(:)); maxregs= max(reg_crosses(:));
Rs = (reg_crosses-minRs)/(maxregs-minRs);  
minCs =  min(recon_errors(:)); maxCs =  max(recon_errors(:)); 
Cs = (recon_errors-minCs)/(maxCs-minCs); 
nAlphas = length(alphas);
nLambdas = length(lambdas);

figure;
marginX = 0.1;
marginY = 0.1;
for li = 1:nLambdas
    for ai = 1:nAlphas
        ax = subplot('Position', [marginX+(ai-1)*(1-marginX)/nAlphas, (nLambdas-li)*(1-marginY)/nLambdas, (1-marginX)/nAlphas, (1-marginY)/nLambdas]);
%         b = bar([Rs(li,ai), Cs(li,ai)],'FaceColor','flat');
%         b.CData(2,:) = [1 0 0];
        boxplot(ax, [Rs(:,li,ai),Cs(:,li,ai)], 'Colors', 'br')
        set(gca,'YLim',[0,1],'XTick',[],'YTick',[], ...
            'Position', [marginX+(ai-1)*(1-marginX)/nAlphas, (nLambdas-li)*(1-marginY)/nLambdas, (1-marginX)/nAlphas, (1-marginY)/nLambdas])
    end
end

for li = 1:nLambdas
    annotation('textbox',pos_Lambdas(li,:),'String',sprintf('\\lambda=%0.3e',lambdas(li)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
for ai = 1:nAlphas
    annotation('textbox',pos_Alphas(ai,:),'String',sprintf('\\alpha=%0.3e',alphas(ai)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
set(gcf, 'Units', 'normalized', 'Position',[.1 .1 .7 .7])

%% Plot running times and scores
f2 = figure;
imagesc(squeeze(mean(cell2mat(times), 1)))
set(gca, 'Position', [0.1 0.1 0.7 0.7], 'XTickLabel', [], 'YTickLabel', []);
colorbar('Position', [0.85 0.1 0.05 0.7]);
pos_Lambdas = [0.1-0.7/nAlphas*ones(nLambdas,1), 0.1+(nLambdas-1:-1:0)'*0.7/nLambdas,... 
    0.7/nAlphas*ones(nLambdas,1), 0.7/nLambdas*ones(nLambdas,1)];
pos_Alphas = [0.1+(0:nAlphas-1)'*0.7/nAlphas, 0.8*ones(nAlphas,1),...
    0.7/nAlphas*ones(nLambdas,1), 0.7/nLambdas*ones(nLambdas,1)];
for li = 1:nLambdas
    annotation('textbox',pos_Lambdas(li,:),'String',sprintf('\\lambda=%0.2e',lambdas(li)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
for ai = 1:nAlphas
    annotation('textbox',pos_Alphas(ai,:),'String',sprintf('\\alpha=%0.2e',alphas(ai)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
annotation('textbox',[0.1 0.85 0.8 0.1],'String', 'Running time(s)', 'FontSize', 16, 'FontWeight', 'bold',...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
set(gcf, 'Units', 'normalized', 'Position',[.1 .1 .7 .7])
save2pdf('times_clean.pdf', f2)

% Number of detected motifs
f3 = figure;
K = 5;
imagesc(squeeze(mean(cell2mat(num_detected), 1))/K, [0.8,1])
set(gca, 'Position', [0.1 0.1 0.7 0.7], 'XTickLabel', [], 'YTickLabel', []);
colorbar('Position', [0.85 0.1 0.05 0.7]);

for li = 1:nLambdas
    annotation('textbox',pos_Lambdas(li,:),'String',sprintf('\\lambda=%0.2e',lambdas(li)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
for ai = 1:nAlphas
    annotation('textbox',pos_Alphas(ai,:),'String',sprintf('\\alpha=%0.2e',alphas(ai)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
annotation('textbox',[0.1 0.85 0.8 0.1],'String', 'Detected motifs(%)', 'FontSize', 16, 'FontWeight', 'bold',...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
set(gcf, 'Units', 'normalized', 'Position',[.1 .1 .7 .7])
save2pdf('nMotifs_clean.pdf', f3)

% Coefficient scores to ground truth
f4 = figure;
K = 5;
imagesc(squeeze(mean(cellfun(@mean,scores), 1)), [0.5,1])
set(gca, 'Position', [0.1 0.1 0.7 0.7], 'XTickLabel', [], 'YTickLabel', []);
colorbar('Position', [0.85 0.1 0.05 0.7]);

for li = 1:nLambdas
    annotation('textbox',pos_Lambdas(li,:),'String',sprintf('\\lambda=%0.2e',lambdas(li)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
for ai = 1:nAlphas
    annotation('textbox',pos_Alphas(ai,:),'String',sprintf('\\alpha=%0.2e',alphas(ai)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
annotation('textbox',[0.1 0.85 0.8 0.1],'String', 'Coefficients', 'FontSize', 16, 'FontWeight', 'bold',...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
set(gcf, 'Units', 'normalized', 'Position',[.1 .1 .7 .7])
save2pdf('scores_clean.pdf', f4)

%% Ratio of Reg/Recon
f5 = figure;
ratios = zeros(li,ai);
for li = 1:nLambdas
    for ai = 1:nAlphas
        ratios(li,ai) = median(regs(:,li,ai))/median(recon_errors(:,li,ai));
    end
end
imagesc(ratios)
set(gca, 'Position', [0.1 0.1 0.7 0.7], 'XTickLabel', [], 'YTickLabel', [], 'ColorScale', 'log');
colorbar('Position', [0.85 0.1 0.05 0.7]);

for li = 1:nLambdas
    annotation('textbox',pos_Lambdas(li,:),'String',sprintf('\\lambda=%0.3e',lambdas(li)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
for ai = 1:nAlphas
    annotation('textbox',pos_Alphas(ai,:),'String',sprintf('\\alpha=%0.3e',alphas(ai)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end  
annotation('textbox',[0.1 0.85 0.8 0.1],'String', 'Reg Error/Recon Error', 'FontSize', 16, 'FontWeight', 'bold',...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
set(f5, 'Units', 'normalized', 'Position',[.1 .1 .7 .7])
save2pdf('ratios_clean.pdf', f5)

%% Plot X_hat
li = 9;
ai = 7;
plotAll = 1;
figure; SimpleWHPlot_patch(W_hats{1,li,ai},H_hats{1,li,ai},[],[],[],[],plotAll); 
set(gcf,'position',[200,200,1200,900])
title('FlexMF reconstruction')

figure; SimpleWHPlot_patch(Ws{1},Hs{1},[],[],[],[],plotAll); 
set(gcf,'position',[200,200,1200,900])
title('Ground Truth')

% %% Similarity of the two results
% for li = 1:length(lambdas)
%     for ai = 1:length(alphas)
%         scores(li,ai) = helper.similarity(W_hats{li,ai}, H_hats{li,ai}, W_hats_dev{li,ai}, H_hats_dev{li,ai});
%     end
% end
% 
% figure;
% set(gcf,'position',[200,100,1000,1000])
% imagesc(scores)
% set(gca, 'Position', [0.1 0.1 0.7 0.7], 'XTickLabel', [], 'YTickLabel', []);
% colorbar('Position', [0.85 0.1 0.05 0.7]);
% pos_Lambdas = [0.1-0.7/nAlphas*ones(nLambdas,1), 0.1+(nLambdas-1:-1:0)'*0.7/nLambdas,... 
%     0.7/nAlphas*ones(nLambdas,1), 0.7/nLambdas*ones(nLambdas,1)];
% pos_Alphas = [0.1+(0:nAlphas-1)'*0.7/nAlphas, 0.8*ones(nAlphas,1),...
%     0.7/nAlphas*ones(nLambdas,1), 0.7/nLambdas*ones(nLambdas,1)];
% for li = 1:nLambdas
%     annotation('textbox',pos_Lambdas(li,:),'String',sprintf('\\lambda=%0.2e',lambdas(li)),...
%         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
% end
% for ai = 1:nAlphas
%     annotation('textbox',pos_Alphas(ai,:),'String',sprintf('\\alpha=%0.2e',alphas(ai)),...
%         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
% end
% annotation('textbox',[0.1 0.85 0.8 0.1],'String', 'Similarity score', 'FontSize', 16, 'FontWeight', 'bold',...
%         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')