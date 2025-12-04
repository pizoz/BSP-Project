function [pcM, Ucomp] = findPrinComp(Y,blocksize)


[nrow,ncol]=size(Y);

if ncol < nrow % make sure Y is a row matrix
    Y=Y';
    temp=nrow;
    nrow=ncol;
    ncol=temp;
end
    
nblocks=ncol/blocksize;

pcM=zeros(nrow,ncol);
Yred=zeros(ncol,nrow);
Ucomp=zeros(ncol,nrow);
kset=1:4;
epsilon=1e-5;

Ynormalized=Y';

i=1;
[U,S,V]=svd(Ynormalized((i-1)*blocksize+1:i*blocksize,:));
Yrot= U(:,1:4)'*Ynormalized((i-1)*blocksize+1:i*blocksize,:);
[aa,bb]=max(abs(Yrot),[],2);

sg=ones(1,4);
for j=1:4
   sg(j)=sign(Yrot(j,bb(j)));
end

Ucomp((i-1)*blocksize+1:i*blocksize,:)=U(:,1:4)*diag(sg);
Yrot= Yrot*diag(1./(diag(S)+epsilon));
Yred((i-1)*blocksize+1:i*blocksize,:)=U(:,kset)*Yrot(kset,:);
[pc,score,latent] = princomp(Y(:,(i-1)*blocksize+1:i*blocksize));
pcM(1:4,(i-1)*blocksize+1:i*blocksize)=pc(:,1:4)';


for i = 2:nblocks
    [U,S,V]=svd(Ynormalized((i-1)*blocksize+1:i*blocksize,:));
    Yrot= U(:,1:4)'*Ynormalized((i-1)*blocksize+1:i*blocksize,:)*diag(1./(diag(S)+epsilon));
    
    for j=1:4
        sg(j)=sign(Yrot(j,bb(j)));
    end
    
    Ucomp((i-1)*blocksize+1:i*blocksize,:)=U(:,1:4)*diag(sg);
    
    
    Yred((i-1)*blocksize+1:i*blocksize,:)=U(:,kset)*Yrot(kset,:);
    pc = princomp(Y(:,(i-1)*blocksize+1:i*blocksize));
    pcM(1:4,(i-1)*blocksize+1:i*blocksize)=pc(:,1:4)';
end

Ucomp=Ucomp';
end