clear all
close all
clc
%% plot costs as a function of lambda and alpha
load choose_lambda_dev.mat;
W_hats_dev = W_hats;
H_hats_dev = H_hats;
nLambdas = length(lambdas); 
nAlphas = length(alphas);
regularization = regularization.*repmat(lambdas', 1, nAlphas);
maxRs= max(regularization(:));

figure;
set(gcf,'position',[200,200,900,900])
marginX = 0.1;
marginY = 0.1;
for li = 1:nLambdas
    for ai = 1:nAlphas
        subplot('Position', [marginX+(ai-1)*(1-marginX)/nAlphas, (nLambdas-li)*(1-marginY)/nLambdas, (1-marginX)/nAlphas, (1-marginY)/nLambdas]);
        b = bar([regularization(li,ai), costs(li,ai)],'FaceColor','flat');
        b.CData(2,:) = [1 0 0];
        set(gca,'YLim',[0,maxRs],'XTick',[],'YTick',[])
    end
end
dim_Lambdas = [zeros(nLambdas,1), (nLambdas-1:-1:0)'*(1-marginY)/nLambdas,... 
    (1-marginX)/nAlphas*ones(nLambdas,1), (1-marginY)/nLambdas*ones(nLambdas,1)];
dim_Alphas = [marginX+(0:nAlphas-1)'*(1-marginX)/nAlphas, (1-marginY)*ones(nAlphas,1),...
    (1-marginX)/nAlphas*ones(nAlphas,1), (1-marginY)/nLambdas*ones(nAlphas,1)];
for li = 1:nLambdas
    annotation('textbox',dim_Lambdas(li,:),'String',sprintf('\\lambda=%0.3e',lambdas(li)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
for ai = 1:nAlphas
    annotation('textbox',dim_Alphas(ai,:),'String',sprintf('\\alpha=%0.3e',alphas(ai)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end     

% Plot costs for single alpha
windowSize = 3; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
alpha_sel = 1e-3;
sel = logical(alphas==alpha_sel);
Rs = filtfilt(b,a,regularization(:,sel)); 
Cs = filtfilt(b,a,costs(:,sel)); 

figure; clf; hold on
plot(lambdas,Rs, 'b')
plot(lambdas,Cs,'r')
scatter(lambdas, regularization(:,sel), 'b', 'markerfacecolor', 'flat');
scatter(lambdas, costs(:,sel), 'r', 'markerfacecolor', 'flat');
xlabel('Lambda'); ylabel('Cost')
set(legend('Correlation cost', 'Reconstruction cost'), 'Box', 'on')
set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
title(['alpha=', num2str(alpha_sel)])

% Regularize costs and remake plots
minRs = prctile(regularization(:,sel),10); maxRs= prctile(regularization(:,sel),90);
Rs = (Rs-minRs)/(maxRs-minRs); 
R = (regularization(:,sel)-minRs)/(maxRs-minRs); 
minCs =  prctile(costs(:,sel),10); maxCs =  prctile(costs(:,sel),90); 
Cs = (Cs -minCs)/(maxCs-minCs); 
C = (costs(:,sel) -minCs)/(maxCs-minCs); 

figure; clf; hold on
plot(lambdas,Rs, 'b')
plot(lambdas,Cs,'r')
scatter(lambdas, R, 'b', 'markerfacecolor', 'flat');
scatter(lambdas, C, 'r', 'markerfacecolor', 'flat');
xlabel('Lambda'); ylabel('Cost (au)')
set(legend('Regularized Correlation cost', 'Regularized Reconstruction cost'), 'Box', 'on')
set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
title(['alpha=', num2str(alpha_sel)])

minRs = min(regularization(:)); maxRs= max(regularization(:));
R = (regularization-minRs)/(maxRs-minRs);  
minCs =  min(costs(:)); maxCs =  max(costs(:)); 
C = (costs-minCs)/(maxCs-minCs); 
nAlphas = length(alphas);
nLambdas = length(lambdas);

figure;
set(gcf,'position',[200,200,900,900])
marginX = 0.1;
marginY = 0.1;
for li = 1:nLambdas
    for ai = 1:nAlphas
        subplot('Position', [marginX+(ai-1)*(1-marginX)/nAlphas, (nLambdas-li)*(1-marginY)/nLambdas, (1-marginX)/nAlphas, (1-marginY)/nLambdas]);
        b = bar([R(li,ai), C(li,ai)],'FaceColor','flat');
        b.CData(2,:) = [1 0 0];
        set(gca,'YLim',[0,1],'XTick',[],'YTick',[])
    end
end
dim_Lambdas = [zeros(nLambdas,1), (nLambdas-1:-1:0)'*(1-marginY)/nLambdas,... 
    (1-marginX)/nAlphas*ones(nLambdas,1), (1-marginY)/nLambdas*ones(nLambdas,1)];
dim_Alphas = [marginX+(0:nAlphas-1)'*(1-marginX)/nAlphas, (1-marginY)*ones(nAlphas,1),...
    (1-marginX)/nAlphas*ones(nAlphas,1), (1-marginY)/nLambdas*ones(nAlphas,1)];
for li = 1:nLambdas
    annotation('textbox',dim_Lambdas(li,:),'String',sprintf('\\lambda=%0.3e',lambdas(li)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
for ai = 1:nAlphas
    annotation('textbox',dim_Alphas(ai,:),'String',sprintf('\\alpha=%0.3e',alphas(ai)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end

%% Plot running times and scores
figure;
set(gcf,'position',[200,100,1000,1000])
imagesc(times)
set(gca, 'Position', [0.1 0.1 0.7 0.7], 'XTickLabel', [], 'YTickLabel', []);
colorbar('Position', [0.85 0.1 0.05 0.7]);
dim_Lambdas = [0.1-0.7/nAlphas*ones(nLambdas,1), 0.1+(nLambdas-1:-1:0)'*0.7/nLambdas,... 
    0.7/nAlphas*ones(nLambdas,1), 0.7/nLambdas*ones(nLambdas,1)];
dim_Alphas = [0.1+(0:nAlphas-1)'*0.7/nAlphas, 0.8*ones(nAlphas,1),...
    0.7/nAlphas*ones(nLambdas,1), 0.7/nLambdas*ones(nLambdas,1)];
for li = 1:nLambdas
    annotation('textbox',dim_Lambdas(li,:),'String',sprintf('\\lambda=%0.2e',lambdas(li)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
for ai = 1:nAlphas
    annotation('textbox',dim_Alphas(ai,:),'String',sprintf('\\alpha=%0.2e',alphas(ai)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
annotation('textbox',[0.1 0.85 0.8 0.1],'String', 'Running time(s)', 'FontSize', 16, 'FontWeight', 'bold',...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
    
figure;
set(gcf,'position',[200,100,1000,1000])
imagesc(scores)
set(gca, 'Position', [0.1 0.1 0.7 0.7], 'XTickLabel', [], 'YTickLabel', []);
colorbar('Position', [0.85 0.1 0.05 0.7]);
dim_Lambdas = [0.1-0.7/nAlphas*ones(nLambdas,1), 0.1+(nLambdas-1:-1:0)'*0.7/nLambdas,... 
    0.7/nAlphas*ones(nLambdas,1), 0.7/nLambdas*ones(nLambdas,1)];
dim_Alphas = [0.1+(0:nAlphas-1)'*0.7/nAlphas, 0.8*ones(nAlphas,1),...
    0.7/nAlphas*ones(nLambdas,1), 0.7/nLambdas*ones(nLambdas,1)];
for li = 1:nLambdas
    annotation('textbox',dim_Lambdas(li,:),'String',sprintf('\\lambda=%0.2e',lambdas(li)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
for ai = 1:nAlphas
    annotation('textbox',dim_Alphas(ai,:),'String',sprintf('\\alpha=%0.2e',alphas(ai)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
annotation('textbox',[0.1 0.85 0.8 0.1],'String', 'Similarity score(to ground truth)', 'FontSize', 16, 'FontWeight', 'bold',...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')

%% Load another group of results
load choose_lambda.mat;
nLambdas = length(lambdas); 
nAlphas = length(alphas);
regularization = regularization.*repmat(lambdas', 1, nAlphas);
maxRs= max(regularization(:));

figure;
set(gcf,'position',[200,200,900,900])
marginX = 0.1;
marginY = 0.1;
for li = 1:nLambdas
    for ai = 1:nAlphas
        subplot('Position', [marginX+(ai-1)*(1-marginX)/nAlphas, (nLambdas-li)*(1-marginY)/nLambdas, (1-marginX)/nAlphas, (1-marginY)/nLambdas]);
        b = bar([regularization(li,ai), costs(li,ai)],'FaceColor','flat');
        b.CData(2,:) = [1 0 0];
        set(gca,'YLim',[0,maxRs],'XTick',[],'YTick',[])
    end
end
dim_Lambdas = [zeros(nLambdas,1), (nLambdas-1:-1:0)'*(1-marginY)/nLambdas,... 
    (1-marginX)/nAlphas*ones(nLambdas,1), (1-marginY)/nLambdas*ones(nLambdas,1)];
dim_Alphas = [marginX+(0:nAlphas-1)'*(1-marginX)/nAlphas, (1-marginY)*ones(nAlphas,1),...
    (1-marginX)/nAlphas*ones(nAlphas,1), (1-marginY)/nLambdas*ones(nAlphas,1)];
for li = 1:nLambdas
    annotation('textbox',dim_Lambdas(li,:),'String',sprintf('\\lambda=%0.3e',lambdas(li)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
for ai = 1:nAlphas
    annotation('textbox',dim_Alphas(ai,:),'String',sprintf('\\alpha=%0.3e',alphas(ai)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end     

% Plot costs for single alpha
windowSize = 3; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
sel = logical(alphas==alpha_sel);
Rs = filtfilt(b,a,regularization(:,sel)); 
Cs = filtfilt(b,a,costs(:,sel)); 

figure; clf; hold on
plot(lambdas,Rs, 'b')
plot(lambdas,Cs,'r')
scatter(lambdas, regularization(:,sel), 'b', 'markerfacecolor', 'flat');
scatter(lambdas, costs(:,sel), 'r', 'markerfacecolor', 'flat');
xlabel('Lambda'); ylabel('Cost')
set(legend('Correlation cost', 'Reconstruction cost'), 'Box', 'on')
set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
title(['alpha=', num2str(alpha_sel)])

% Regularize costs and remake plots
minRs = prctile(regularization(:,sel),10); maxRs= prctile(regularization(:,sel),90);
Rs = (Rs-minRs)/(maxRs-minRs); 
R = (regularization(:,sel)-minRs)/(maxRs-minRs); 
minCs =  prctile(costs(:,sel),10); maxCs =  prctile(costs(:,sel),90); 
Cs = (Cs -minCs)/(maxCs-minCs); 
C = (costs(:,sel) -minCs)/(maxCs-minCs); 

figure; clf; hold on
plot(lambdas,Rs, 'b')
plot(lambdas,Cs,'r')
scatter(lambdas, R, 'b', 'markerfacecolor', 'flat');
scatter(lambdas, C, 'r', 'markerfacecolor', 'flat');
xlabel('Lambda'); ylabel('Cost (au)')
set(legend('Regularized Correlation cost', 'Regularized Reconstruction cost'), 'Box', 'on')
set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
title(['alpha=', num2str(alpha_sel)])

minRs = min(regularization(:)); maxRs= max(regularization(:));
R = (regularization-minRs)/(maxRs-minRs);  
minCs =  min(costs(:)); maxCs =  max(costs(:)); 
C = (costs-minCs)/(maxCs-minCs); 
nAlphas = length(alphas);
nLambdas = length(lambdas);

figure;
set(gcf,'position',[200,200,900,900])
marginX = 0.1;
marginY = 0.1;
for li = 1:nLambdas
    for ai = 1:nAlphas
        subplot('Position', [marginX+(ai-1)*(1-marginX)/nAlphas, (nLambdas-li)*(1-marginY)/nLambdas, (1-marginX)/nAlphas, (1-marginY)/nLambdas]);
        b = bar([R(li,ai), C(li,ai)],'FaceColor','flat');
        b.CData(2,:) = [1 0 0];
        set(gca,'YLim',[0,1],'XTick',[],'YTick',[])
    end
end
dim_Lambdas = [zeros(nLambdas,1), (nLambdas-1:-1:0)'*(1-marginY)/nLambdas,... 
    (1-marginX)/nAlphas*ones(nLambdas,1), (1-marginY)/nLambdas*ones(nLambdas,1)];
dim_Alphas = [marginX+(0:nAlphas-1)'*(1-marginX)/nAlphas, (1-marginY)*ones(nAlphas,1),...
    (1-marginX)/nAlphas*ones(nAlphas,1), (1-marginY)/nLambdas*ones(nAlphas,1)];
for li = 1:nLambdas
    annotation('textbox',dim_Lambdas(li,:),'String',sprintf('\\lambda=%0.3e',lambdas(li)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
for ai = 1:nAlphas
    annotation('textbox',dim_Alphas(ai,:),'String',sprintf('\\alpha=%0.3e',alphas(ai)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end

%% Plot running times and scores
figure;
set(gcf,'position',[200,100,1000,1000])
imagesc(times)
set(gca, 'Position', [0.1 0.1 0.7 0.7], 'XTickLabel', [], 'YTickLabel', []);
colorbar('Position', [0.85 0.1 0.05 0.7]);
dim_Lambdas = [0.1-0.7/nAlphas*ones(nLambdas,1), 0.1+(nLambdas-1:-1:0)'*0.7/nLambdas,... 
    0.7/nAlphas*ones(nLambdas,1), 0.7/nLambdas*ones(nLambdas,1)];
dim_Alphas = [0.1+(0:nAlphas-1)'*0.7/nAlphas, 0.8*ones(nAlphas,1),...
    0.7/nAlphas*ones(nLambdas,1), 0.7/nLambdas*ones(nLambdas,1)];
for li = 1:nLambdas
    annotation('textbox',dim_Lambdas(li,:),'String',sprintf('\\lambda=%0.2e',lambdas(li)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
for ai = 1:nAlphas
    annotation('textbox',dim_Alphas(ai,:),'String',sprintf('\\alpha=%0.2e',alphas(ai)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
annotation('textbox',[0.1 0.85 0.8 0.1],'String', 'Running time(s)', 'FontSize', 16, 'FontWeight', 'bold',...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
    
figure;
set(gcf,'position',[200,100,1000,1000])
imagesc(scores)
set(gca, 'Position', [0.1 0.1 0.7 0.7], 'XTickLabel', [], 'YTickLabel', []);
colorbar('Position', [0.85 0.1 0.05 0.7]);
dim_Lambdas = [0.1-0.7/nAlphas*ones(nLambdas,1), 0.1+(nLambdas-1:-1:0)'*0.7/nLambdas,... 
    0.7/nAlphas*ones(nLambdas,1), 0.7/nLambdas*ones(nLambdas,1)];
dim_Alphas = [0.1+(0:nAlphas-1)'*0.7/nAlphas, 0.8*ones(nAlphas,1),...
    0.7/nAlphas*ones(nLambdas,1), 0.7/nLambdas*ones(nLambdas,1)];
for li = 1:nLambdas
    annotation('textbox',dim_Lambdas(li,:),'String',sprintf('\\lambda=%0.2e',lambdas(li)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
for ai = 1:nAlphas
    annotation('textbox',dim_Alphas(ai,:),'String',sprintf('\\alpha=%0.2e',alphas(ai)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
annotation('textbox',[0.1 0.85 0.8 0.1],'String', 'Similarity score(to ground truth)', 'FontSize', 16, 'FontWeight', 'bold',...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')

%% Similarity of the two results
for li = 1:length(lambdas)
    for ai = 1:length(alphas)
        scores(li,ai) = helper.similarity(W_hats{li,ai}, H_hats{li,ai}, W_hats_dev{li,ai}, H_hats_dev{li,ai});
    end
end

figure;
set(gcf,'position',[200,100,1000,1000])
imagesc(scores)
set(gca, 'Position', [0.1 0.1 0.7 0.7], 'XTickLabel', [], 'YTickLabel', []);
colorbar('Position', [0.85 0.1 0.05 0.7]);
dim_Lambdas = [0.1-0.7/nAlphas*ones(nLambdas,1), 0.1+(nLambdas-1:-1:0)'*0.7/nLambdas,... 
    0.7/nAlphas*ones(nLambdas,1), 0.7/nLambdas*ones(nLambdas,1)];
dim_Alphas = [0.1+(0:nAlphas-1)'*0.7/nAlphas, 0.8*ones(nAlphas,1),...
    0.7/nAlphas*ones(nLambdas,1), 0.7/nLambdas*ones(nLambdas,1)];
for li = 1:nLambdas
    annotation('textbox',dim_Lambdas(li,:),'String',sprintf('\\lambda=%0.2e',lambdas(li)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
for ai = 1:nAlphas
    annotation('textbox',dim_Alphas(ai,:),'String',sprintf('\\alpha=%0.2e',alphas(ai)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
annotation('textbox',[0.1 0.85 0.8 0.1],'String', 'Similarity score', 'FontSize', 16, 'FontWeight', 'bold',...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')