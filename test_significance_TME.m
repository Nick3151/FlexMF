function [pvals,is_significant] = test_significance_TME(TestData,W,p,nnull)
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
indempty = indempty | (max(Wflat,[],1).^2> .999*sum(Wflat.^2,1)); % or one neuron has >99.9% of the power
W(:,indempty,:) = []; % Delete factors that meet the above critera

[N,K,L] = size(W);
[~,T] = size(TestData);

if nargin < 3
    p = 0.05;
end

if nargin < 4
    nnull = ceil(K/p)*2;
end

rng('shuffle', 'twister') % randomize the seed
surrogate_type = 'surrogate-TNC';

% make nnull shifted datasets 
skewnull = zeros(K,nnull);

% Sample surrogate data with Tensor-Maximum-Entropy (TME) 
X = log(TestData' + eps);    % Transform data to log space
[targetSigmaT, targetSigmaN, targetSigmaC, M] = extractFeatures(X);
figure; imagesc(targetSigmaT)
figure; imagesc(targetSigmaN)

params = [];
if strcmp(surrogate_type, 'surrogate-T')
    params.margCov{1} = targetSigmaT;
    params.margCov{2} = [];
    params.margCov{3} = [];
    params.meanTensor = M.T;
elseif strcmp(surrogate_type, 'surrogate-TN')
    params.margCov{1} = targetSigmaT;
    params.margCov{2} = targetSigmaN;
    params.margCov{3} = [];
    params.meanTensor = M.TN;
elseif strcmp(surrogate_type, 'surrogate-TNC')
    params.margCov{1} = targetSigmaT;
    params.margCov{2} = targetSigmaN;
    params.margCov{3} = targetSigmaC;
    params.meanTensor = M.TNC; 
else
    error('please specify a correct surrogate type') 
end

maxEntropy = fitMaxEntropy(params);             % fit the maximum entropy distribution

for n = 1:nnull
    [Xnull] = sampleTME(maxEntropy);
    [SurrSigmaT, SurrSigmaN, SurrSigmaC, ~] = extractFeatures(Xnull);
    SurrData = exp(Xnull');
    WTX = helper.transconv(W,SurrData);
    if n<5
        maxValue = prctile(abs(SurrData),99,'all')+eps;
        clims = [0, maxValue]; 
        cmap = flipud(gray); 
        figure; imagesc(SurrData, clims);
        colormap(cmap)
        set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
        figure; histogram(WTX(1,:),'FaceColor','k')
        set(gca,'yscale','log')
    end
    % Get skewness of each
    skewnull(:,n) = skewness(WTX,1,2);
end


WTX = helper.transconv(W,TestData);
skew = skewness(WTX,1,2);
for k = 1:K
    % Assign pvals from skewness
    pvals(k) = (1+sum(skewnull(k,:)>skew(k)))/nnull;
end
figure;
histogram(skewnull(1,:))
hold on
xline(skew(1), 'Color', 'r', 'Linewidth', 2);

allpvals(indempty) = Inf; 
allpvals(~indempty) = pvals; 
pvals = allpvals;
is_significant = (pvals <= p/K);
end