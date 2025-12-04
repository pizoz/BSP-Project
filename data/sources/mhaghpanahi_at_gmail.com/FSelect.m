function [Forder,avgKurt] = FSelect(X,blocksize)

[nrow,ncol]=size(X); % make sure X is a column matrix 

if nrow < ncol
    X=X';
end

nblocks = max(nrow,ncol)/blocksize;
% range=randi(nblocks,1,checkno);

% maxInd= zeros(1,checkno);
% kurtMeas= zeros(1,checkno);

maxInd= zeros(1,nblocks);
kurtMeas= zeros(1,nblocks);

for j=1:nblocks
    ku = kurtosis(X((j-1)*blocksize+1:j*blocksize,:));
    [kurtMeas(j), maxInd(j)] = max(ku);
end

compare= hist(maxInd,1:min(nrow,ncol));
[dummy, order]=sort(compare,'descend');

Forder=order(1);

avgKurt=mean(kurtMeas(maxInd==Forder));