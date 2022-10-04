function score = similarity(W1,H1,W2,H2)
% Measure the similarity between two factorizations 
% by matching each factor
%%
% W1 = W1.*(W1 > 0.001);
% H1 = H1.*(H1 > 0.001);

X1 = {};
X2 = {};

jj = 1;
for ii = 1:size(W1,2)
    X = helper.reconstruct(W1(:,ii,:),H1(ii,:));
    if(any(X(:)))
        X1{jj} = X;
        jj = jj+1;
    end
end

jj = 1;
for ii = 1:size(W2,2)
    X = helper.reconstruct(W2(:,ii,:),H2(ii,:));
    if(any(X(:)))
        X2{jj} = X;
        jj = jj+1;
    end
end
%%
for ii = 1:size(X1,2)
    for jj = 1:size(X2,2)
        t1 = X1{ii};
        t2 = X2{jj};
        S(ii,jj) = (t1(:)'*t2(:))/((sqrt(t1(:)'*t1(:))*sqrt(t2(:)'*t2(:)))+eps);
    end
end
% S(isnan(S)) = 0;
% S(S<0) = eps;
%%
temp = S;
num = 0;
K = min(size(X1,2), size(X2,2));
% Matching each non-zero factor to all the factors of another reconstruction
for ii = 1:K
    [r,c]= find(temp == max(temp(:)));
    maximum = temp(r(1), c(1));
    num = num + maximum;
    temp(r(1),:) = 0;
    temp(:,c(1)) = 0;

end
score = num/(K+eps);
% score(isnan(score)) = 0;
end