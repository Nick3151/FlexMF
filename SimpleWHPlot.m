function SimpleWHPlot(W, H, trials, frames, is_significant, neg, Data, plotAll) 
% plots factors W and H with trials
% Also plots Data if provided, and reconstruction if data is not provided
% plotAll=1 means plot all data
% plotAll=0 means crop data so it's on the same timescale as W's 
% Emily Mackevicius 2/3/2018
% Adapted by Liang Xiang

clf

% set(gcf, 'color', 'w');
if nargin < 8
    plotAll = 0;
end
if nargin < 7 || size(Data,2)==0
    plotData = 0;
else
    plotData = 1; 
end

if nargin < 6 || neg
    cmap_red = [ones(128,1),linspace(1,0,128)',linspace(1,0,128)'];
    cmap_blue = [linspace(0,1,128)',linspace(0,1,128)',ones(128,1)];
    cmap = [cmap_blue; cmap_red];
else
    cmap = flipud(gray);
end

[N,K,L] = size(W); 
[~,T] = size(H);
if ~isempty(trials) || ~isempty(frames)
    assert(trials*frames==T, 'Dimensions of trials do not match!')
    plotTrials = 1;
else
    plotTrials = 0;
end

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
    indplot = 1:T;
else
    indplot = 2*L+(1:ceil((K*(L+sep))/ww*wdata)); % so that data and W's are on same scale
    indplot(indplot>T) = [];
end


%% plot W's
axW = subplot('Position', [m m ww hdata]);
hold on
set(gca, 'ColorOrder', kColors); 

WsToPlot = zeros(N,K*(L+sep)); 
XsToPlot = zeros(3,K); 
YsToPlot = [N*ones(1,K); zeros(2,K)]+.5;
for ki = 1:K
    WsToPlot(:,((L+sep)*(ki-1)+1):((L+sep)*(ki-1)+L)) = squeeze(W(:,ki,:));
    XsToPlot(:,ki) = [(L+sep)*(ki-1)+1 (L+sep)*(ki-1)+1 (L+sep)*(ki-1)+L];
end

maxValue = prctile(abs(WsToPlot),99,'all')+epsilon;
if neg
    clims = [-maxValue, maxValue]; 
else 
    clims = [0, maxValue]; 
end

imagesc(WsToPlot, clims); 
colormap(cmap)

plot(XsToPlot,YsToPlot);
xlim([1 K*(L+sep)]);ylim([0 N+.1]+.5)
set(gca, 'ydir', 'reverse')
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
if plotData
    maxValue = prctile(abs(Data),99,'all')+epsilon;
    if neg
        clims = [-maxValue, maxValue]; 
    else 
        clims = [0, maxValue]; 
    end
    imagesc(Data(:,indplot), clims);
    colormap(cmap)
else
    toplot = helper.reconstruct(W,H(:,indplot));
    maxValue = prctile(abs(toplot),99.9,'all')+epsilon;
    if neg
        clims = [-maxValue, maxValue]; 
    else 
        clims = [0, maxValue]; 
    end
    imagesc(toplot,clims);
    colormap(cmap)
end

hold on
if plotTrials
    trial_start = ceil(indplot(1)/frames);
    trial_end = floor(indplot(end)/frames);
    for trial=trial_start:trial_end
        xline(trial*frames-indplot(1)+1, 'LineWidth', 2);
    end
end
set(gca,'ydir','reverse')
plot([0 0 length(indplot)+1], [N 0 0]+.5, 'k')
xlim([0 length(indplot)+1]);ylim([0 N+.1]+.5)
axis off
% %% plot Wflat (collapse out L dimension of W)
% axWflat = subplot('Position', [m+ww+wdata m wwflat hdata]);
% hold on
% set(gca, 'ColorOrder', kColors); 
% plot(squeeze(sum(W,3)), 1:N,'>', 'markersize', 2.5);
% ylim([0 N+.1]+.5)
% axis tight
% % xlims = xlim; 
% % xlim([xlims(2)*.1 xlims(2)])
% set(gca, 'ydir', 'reverse')
% axis off
%% plot H's
axH = subplot('Position', [m+ww m+hdata wdata hh]);
Hrescaled = repmat(squeeze(sum(sum(abs(W),1),3))',1,T).*H; % rescale by approximate loading
dn = prctile(Hrescaled(:),100)/2; 
for ki = K:-1:1
    Xs = [1 1:length(indplot) length(indplot)]; 
    Ys = [dn*ki (dn*ki + Hrescaled(K-ki+1,indplot)) dn*ki]-dn/2;
    patch(Xs,Ys, kColors(K-ki+1,:), 'edgecolor', kColors(K-ki+1,:))
    hold on
end
ylim([0 dn*K+dn*3+epsilon]);xlim([0 length(indplot)+1])
axis off
%%
if plotAll
      linkaxes([axIm axW], 'y'); linkaxes([axIm axH], 'x');
end

