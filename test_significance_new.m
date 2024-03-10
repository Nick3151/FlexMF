function [pvals,is_significant,is_single] = test_significance_new(TestData, W,plot,p,nnull)
%
% USAGE: 
%
% [pvals,is_significant] = test_significance(TestData,W,0.01)
%
% ------------------------------------------------------------------------
% DESCRIPTION:
%
% Tests each factor in W for significance using a held out test dataset at
% a p-value of p using Bonferroni correction. 
%  
% ------------------------------------------------------------------------
%
% INPUTS:
%
% Name              Default                 Description
% TestData                                  Held out data matrix (NxT) 
% W                                         NxKxL tensor containing factor exemplars
% plot              0                       If plot distributions
% p                 0.05                    Desired p-value to test
% nnull             ceil(K/p)*2             Number of null datasets to use
%
% ------------------------------------------------------------------------
% OUTPUTS:
%
% pvals                      A vector (1xK) containing the p-value of each factor
% is_significant             A boolean vector (1xK) which is 1 if a factor
%                            is significant using the specified pvalue with correction
%
% ------------------------------------------------------------------------
% CREDITS:
%   Emily Mackevicius and Andrew Bahle, 2/1/2018
%
%   Please cite our paper: 
%       XXXXXXXXXXXXXXXXXXXXX
% Remove factors where there is obviously no sequence
% That is, W is empty, or has one neuron with >99.9% of the power

indempty = sum(sum(W>0,1),3)==0; % W is literally empty
Wflat = sum(W,3); 
is_single = (max(Wflat,[],1).^2> .999*sum(Wflat.^2,1)); % one neuron has >99.9% of the power
W(:,indempty|is_single,:) = []; % Delete factors that meet the above critera

[N,K,L] = size(W);
[~,T] = size(TestData);

if nargin < 3 || isempty(plot)
    plot = 0;
end

if nargin < 4
    p = 0.05;
end

if nargin < 5
    nnull = ceil(K/p)*5;
end

WTX = helper.transconv(W,TestData);
thresh = prctile(WTX', 95);
WTX_large_mean = zeros(K,1);
for k=1:K
    WTX_tmp = WTX(k,:);
    WTX_large_mean(k) = mean(WTX_tmp(WTX_tmp>thresh(k)));
end

if plot
    k_plot = 1;
    nbins = 100;
    figure; histogram(WTX(k_plot,:), nbins, 'FaceColor','k'); 
    hold on
    xline(thresh(k_plot), 'r', 'LineWidth',1)
    title('WTX distribution')
    set(gca,'yscale','log')
end

% make nnull shifted datasets 
thresh_null = zeros(K,nnull);
WTX_large_mean_null = zeros(K,nnull);

for n = 1:nnull
    Wnull = zeros(N,K,L);
    for k = 1:K
        for ni = 1:N
            Wnull(ni,k,:) = circshift(W(ni,k,:),[0,0,randi(L)]);
        end
    end
    WTX = helper.transconv(Wnull,TestData);
    thresh_null(:,n) = prctile(WTX', 95);
    for k=1:K
        WTX_tmp = WTX(k,:);
        WTX_large_mean_null(k,n) = mean(WTX_tmp(WTX_tmp>thresh_null(k,n)));
    end
    if plot
        if n<5
            figure; histogram(WTX(k_plot,:), nbins, 'FaceColor','k'); 
            hold on
            xline(thresh_null(k_plot), 'r', 'LineWidth',1)
            title('WTX null distribution')
            set(gca,'yscale','log')
        end
    end
end

pvals = zeros(K,1);
for k = 1:K
    % Assign pvals from mean of top 5% values
    pvals(k) = (1+sum(WTX_large_mean_null(k,:)>WTX_large_mean(k)))/nnull;
end
allpvals(indempty|is_single) = Inf; 
allpvals(~(indempty|is_single)) = pvals; 
pvals = allpvals;
is_significant = (pvals <= p/K);

if plot
    figure;
    histogram(WTX_large_mean_null(k_plot,:))
    hold on
    xline(WTX_large_mean(k_plot), 'Color', 'r', 'Linewidth', 2);
end
end