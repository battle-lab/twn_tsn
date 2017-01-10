function [X] = standardize(X)
[n p] = size(X);
X = X - repmat(nanmean(X,2),1,p);
m = mean(X,2);
X(isnan(X)) = 0;
st = sqrt(sum(X.^2,2));
X = X./repmat(nanstd(X,0,2)+1e-100,1,p);
