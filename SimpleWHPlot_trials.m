function SimpleWHPlot_trials(W, H, is_significant, Data, plotAll) 
% plots factors W and H with trials
% W: N x K x L
% H: K x Trials
% Data: N x L x Trials
% Also plots Data if provided, and reconstruction if data is not provided
% plotAll=1 means plot all data
% plotAll=0 means crop data so it's on the same timescale as W's 
% Emily Mackevicius 2/3/2018
% Adapted by Liang Xiang

clf

% set(gcf, 'color', 'w');
if nargin < 5
    plotAll = 0;
end
if nargin < 4 || isempty(Data)
    plotData = 0;
else
    plotData = 1; 
end

[N,K,L] = size(W); 
[~,Trials] = size(H);

if isempty(is_significant)
    is_significant = zeros(1,K);
else
    assert(length(is_significant)==K, 'Dimensions of is_significant do not match!')
end

color_palet = [[0 .6 .3]; [.7 0 .7]; [1 .6 0];  [.1 .3 .9];  [1 .1 .1];  [0 .9 .3]; [.4 .2 .7]; [.7 .2 .1]; [.1 .8 1 ]; [1 .3 .7]; [.2 .8 .2]; [.7 .4 1]; [.9 .6 .4]; [0 .6 1]; [1 .1 .3]]; 
color_palet = repmat(color_palet, ceil(K/size(color_palet,1)),1); 
kColors = color_palet(1:K,:); 
epsilon = 1e-4;
%% set widths of subplots
m = .05; % margin
ww = min(.05*K, .25); % width of W plot
% wwflat = .05; % width of Wflat plot
hh = max(.05*K, .2); % height of H plot
hdata = 1-hh-2*m; 
wdata = 1-ww-2*m; 
sep = ceil(L*.1); 

%% crop data, unless plotAll
if plotAll
    indplot = 1:Trials*L;
else
    indplot = (1:ceil((K*(L+sep))/ww*wdata)); % so that data and W's are on same scale
    indplot(indplot>Trials*L) = [];
end


%% plot W's
axW = subplot('Position', [m m ww hdata]);
hold on
set(gca, 'ColorOrder', kColors); 

dnW = prctile(W(:),100); 
XsEdge = zeros(3,K); 
YsEdge = [zeros(1,K); dnW*N*ones(2,K)];

for ki=1:K
    for ni =1:N
        Xs = [(ki-1)*(L+sep)+1, (ki-1)*(L+sep)+1:(ki-1)*(L+sep)+L, (ki-1)*(L+sep)+L];
        Ys = [dnW*(N-ni) dnW*(N-ni) + squeeze(W(ni,ki,:))' dnW*(N-ni)];
        patch(Xs,Ys, 'k', 'edgecolor', 'none')
        hold on
    end
    XsEdge(:,ki) = [(L+sep)*(ki-1)+1 (L+sep)*(ki-1)+1 (L+sep)*(ki-1)+L];
end

plot(XsEdge,YsEdge, 'LineWidth', 2);
xlim([1 K*(L+sep)]);ylim([0 dnW*N])

axis off

%% Plot significance of each factor
pos = [(m+(0:K-1)/K*ww)', repmat(m+hdata,K,1), repmat(ww/K,K,1), repmat(0.05,K,1)];
for k=1:K
    if is_significant(k)
        annotation('textbox', pos(k,:), 'string', '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom',...
            'Color', kColors(k,:), 'LineStyle', 'none', 'FontWeight', 'bold', 'FontSize', 14)
    end
end

%% plot data
axIm = subplot('Position', [m+ww m wdata hdata]);
Xs = [1 1:length(indplot) length(indplot)]; 

for ni=1:N
    Dhat = zeros(N,Trials*L);
    if plotData
        assert(isequal(size(Data), [N,L,Trials]), 'Dimensions of X do not match!')
        for t=1:Trials
            Dhat(:,(t-1)*L+1:t*L) = squeeze(Data(:,:,t));
        end
    else
        for t=1:Trials
            for k=1:K
                Dhat(:,(t-1)*L+1:t*L) = Dhat(:,(t-1)*L+1:t*L) + squeeze(W(:,k,:))*H(k,t);
            end
        end
    end

    dnX = prctile(Dhat(:),100); 
    Ys = [dnW*(N-ni) dnW*(N-ni)+Dhat(ni,indplot)/dnX*dnW dnW*(N-ni)];
    patch(Xs,Ys, 'k', 'edgecolor', 'none')
    hold on
end

hold on
trial_start = ceil(indplot(1)/L);
trial_end = floor(indplot(end)/L);
for trial=trial_start:trial_end
    xline(trial*L-indplot(1)+1, 'LineWidth', 2);
end

plot([0 0 length(indplot)+1], [0 dnW*N dnW*N], 'k', 'LineWidth', 2)
xlim([0 length(indplot)+1]);ylim([0 dnW*N])
axis off

%% plot H's
axH = subplot('Position', [m+ww m+hdata wdata hh]);
HToPlot = zeros(K,Trials*L);
HToPlot(:,1:L:Trials*L) = H;
dn = max(HToPlot(:)); 
for ki = K:-1:1
    Xs = [1 1:length(indplot) length(indplot)]; 
    Ys = [dn*ki HToPlot(K-ki+1,1:length(indplot))+dn*ki dn*ki]-dn/2;
    patch(Xs,Ys, kColors(K-ki+1,:), 'edgecolor', kColors(K-ki+1,:))
    hold on
end
ylim([0 dn*(K+1)+epsilon]);xlim([0 length(indplot)+1])
axis off
%%
if plotAll
      linkaxes([axIm axW], 'y'); linkaxes([axIm axH], 'x');
end

