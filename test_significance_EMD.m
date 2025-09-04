function [pvals,is_significant,is_single] = test_significance_EMD(TestData, W, M, varargin)
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
% M                                         Motion field matrix (NxT)
% plot              0                       If plot distributions
% pValue            0.05                    Desired p-value to test
% nNull             ceil(K/p)*2             Number of null datasets to use
% useSkew           1                       Use skewness of WTX or mean of top 5% WTX
%
% ------------------------------------------------------------------------
% OUTPUTS:
%
% pvals                      A vector (1xK) containing the p-value of each factor
% is_significant             A boolean vector (1xK) which is 1 if a factor
%                            is significant using the specified pvalue with correction
% is_single                  A boolean vector (1xK) which is 1 if a factor 
%                            only contains single neuron
%
% ------------------------------------------------------------------------
% CREDITS:
%   Adapted by Liang Xiang
%   Original script developed by Emily Mackevicius and Andrew Bahle, 2/1/2018
%
%   Please cite our paper: 
%       XXXXXXXXXXXXXXXXXXXXX

p = inputParser;
addOptional(p, 'nNull', []);
addOptional(p, 'plot', 0);
addOptional(p, 'pValue', .05);
addOptional(p, 'useSkew', 1);

parse(p, varargin{:});
nNull = p.Results.nNull;
plot = p.Results.plot;
pValue = p.Results.pValue;
useSkew = p.Results.useSkew;

indempty = sum(sum(W>0,1),3)==0; % W is literally empty
Wflat = sum(W,3); 
is_single = (max(Wflat,[],1).^2> .999*sum(Wflat.^2,1)); % one neuron has >99.9% of the power
W(:,indempty|is_single,:) = []; % Delete factors that meet the above critera

[N,K,L] = size(W);
[~,T] = size(TestData);
assert(isequal(size(M), [N,T]), 'Dimension of M should match TestData')

if isempty(nNull)
    nnull = ceil(K/pValue)*5;
end

% Correct the temporal warpped/jittered part of TestData
D = eye(T) - diag(ones(T-1,1),-1);
TestData_corr = M*D'+TestData;

WTX = helper.transconv(W,TestData_corr);
thresh = prctile(WTX', 95);
if useSkew
    skew = skewness(WTX,1,2);
else  
    WTX_large_mean = zeros(K,1);
    for k=1:K
        WTX_tmp = WTX(k,:);
        WTX_large_mean(k) = mean(WTX_tmp(WTX_tmp>thresh(k)));
    end
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
skewnull = zeros(K,nnull);

for n = 1:nnull
    Wnull = zeros(N,K,L);
    for k = 1:K
        for ni = 1:N
            Wnull(ni,k,:) = circshift(W(ni,k,:),[0,0,randi(L)]);
        end
    end
    WTX_null = helper.transconv(Wnull,TestData_corr);
    thresh_null(:,n) = prctile(WTX_null', 95);
    if useSkew
        skewnull(:,n) = skewness(WTX_null,1,2);
    else        
        for k=1:K
            WTX_tmp = WTX_null(k,:);
            WTX_large_mean_null(k,n) = mean(WTX_tmp(WTX_tmp>thresh_null(k,n)));
        end
    end
    if plot
        if n<5
            figure; histogram(WTX_null(k_plot,:), nbins, 'FaceColor','k'); 
            hold on
            xline(thresh_null(k_plot), 'r', 'LineWidth',1)
            title('WTX null distribution')
            set(gca,'yscale','log')
        end
    end
end

pvals = zeros(K,1);
for k = 1:K
    if useSkew
        % Assign pvals from skewness
        pvals(k) = (1+sum(skewnull(k,:)>skew(k)))/nnull;
    else
        % Assign pvals from mean of top 5% values
        pvals(k) = (1+sum(WTX_large_mean_null(k,:)>WTX_large_mean(k)))/nnull;
    end
end
allpvals(indempty|is_single) = Inf; 
allpvals(~(indempty|is_single)) = pvals; 
pvals = allpvals;
is_significant = (pvals <= pValue/K);

if plot
    figure;
    if useSkew
        histogram(skewnull(k_plot,:))
        hold on
        xline(skew(k_plot), 'Color', 'r', 'Linewidth', 2);
    else
        histogram(WTX_large_mean_null(k_plot,:))
        hold on
        xline(WTX_large_mean(k_plot), 'Color', 'r', 'Linewidth', 2);
    end
end
