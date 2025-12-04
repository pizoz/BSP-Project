function final_order = SigSelect(X,fs,checkno)

sub_order=[];
[nrow,ncol]=size(X);

if nrow > ncol
    X=X';
end


blocksize=2000;
nblocks = max(nrow,ncol)/blocksize;
n=min(nrow,ncol);

range=randi(nblocks,1,checkno);
% e=zeros(n,checkno);


F=0:150;
interval=1:38;
maxInd= zeros(1,checkno);

for j=1:checkno
    for i=1:n    
        [S(:,i),F,T,P(:,i)] = spectrogram(X(i,(range(j)-1)*blocksize+1:range(j)*blocksize),blocksize,[],F,fs);
    end
    
    p_ratio =  sum(P(interval,:))./sum(P);

    
    k_meas=kurtosis(X(:,(range(j)-1)*blocksize+1:range(j)*blocksize)');
    [dummy, maxInd(j)]=max(k_meas.*p_ratio);

    clear S T P;
end

compare= hist(maxInd,1:n);
[csort, order]=sort(compare,'descend');

if sum(csort==0)<3 
    final_order = order;
    return;
else
    sub_order = SigSelect(X(order(csort==0),:),fs,checkno);
end

w=order(csort==0);
final_order = [order(csort~=0) w(sub_order)];
    
end
