function C = max_xcorr(X, lag, varargin)
% Max correlation between rows of X
[N,T] = size(X);
p = inputParser;
addOptional(p, 'bin', 0);
addOptional(p, 'option', 'none');
parse(p, varargin{:});
bin = p.Results.bin;
option = p.Results.option;

if bin
    X = (X>0);
end

% nNull = 100;    % Number of shuffled data
% shifts = randi(T, nNull, N);
% XNull = zeros(nNull, N, T);
% for k=1:nNull
%     for n = 1:N
%         XNull(k,n,:) = circshift(X(n,:), [0, shifts(k,n)]);
%     end
% end

C = zeros(N);
for i=1:N
    for j=i+1:N
        C(i,j) = max(xcorr(X(i,:), X(j,:), lag, option));
        C(j,i) = C(i,j);
    end
end

% CNull = zeros(nNull, N, N);
% for k=1:nNull
%     for i=1:N
%         for j=i+1:N
%             CNull(k,i,j) = max(xcorr(squeeze(XNull(k,i,:)), squeeze(XNull(k,j,:)), lag, option));
%             CNull(k,j,i) = CNull(k,i,j);
%         end
%     end
% end
% 
% pval = eye(N);
% for k=1:nNull
%     pval = pval + (squeeze(CNull(k,:,:))>C)/nNull;
% end
