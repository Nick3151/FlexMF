function [pvals,is_significant] = test_significance_trials(TestData, trials, frames, W,plot,p,nnull)
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
indempty = indempty | (max(Wflat,[],1).^2> .999*sum(Wflat.^2,1)); % or one neuron has >99.9% of the power
W(:,indempty,:) = []; % Delete factors that meet the above critera

[N,K,L] = size(W);
[~,T] = size(TestData);
assert(trials*frames==T, 'Dimensions of trials do not match!')

if nargin < 5
    plot = 0;
end

if nargin < 6
    p = 0.05;
end

if nargin < 7
    nnull = ceil(K/p)*2;
end

% make nnull shifted datasets 
skewnull = zeros(K,nnull);

% Normalize each row
% X = bsxfun(@rdivide, TestData, sum(TestData,2));
X = TestData;

for i = 1:nnull
    % Make a null dataset, randomly circuilarly shift each row in each
    % trial, yet keep transient structures.
    Xnull = zeros(N,T);
    for trial = 1:trials
        X_tmp = X(:,(trial-1)*frames+1:trial*frames);  
        [ids_start, ids_end, ~] = locate_true_transients(X_tmp,6,L);
        for n=1:N
            if isempty(ids_start{n}) || isempty(ids_end{n})
                Xnull(n,(trial-1)*frames+1:trial*frames) = circshift(X_tmp(n,:), randi(L));
            else
                Xnull(n,(trial-1)*frames+1:trial*frames) = circshift(X_tmp(n,:), randi([1-ids_start{n}(1), L-ids_end{n}(end)]));
            end
        end
    end
    WTX = helper.transconv(W,Xnull);
    if plot
        if i<5
            SimpleXplot_patch(Xnull, trials, frames);
            set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
            figure; histogram(WTX(1,:),'FaceColor','k'); title('WTX distribution')
            set(gca,'yscale','log')
        end
    end
    % Get skewness of each
    skewnull(:,n) = skewness(WTX,1,2);
end

WTX = helper.transconv(W,X); 
skew = skewness(WTX,1,2);
for k = 1:K
    % Assign pvals from skewness
    pvals(k) = (1+sum(skewnull(k,:)>skew(k)))/nnull;
end
if plot
    figure;
    histogram(skewnull(1,:)); title('skewness distribution')
    hold on
    xline(skew(1), 'Color', 'r', 'Linewidth', 2);
end

allpvals(indempty) = Inf; 
allpvals(~indempty) = pvals; 
pvals = allpvals;
is_significant = (pvals <= p/K);
end