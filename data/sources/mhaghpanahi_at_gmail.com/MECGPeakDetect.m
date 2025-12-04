function [Y,accY,peaks,goodint] = MECGPeakDetect(Y,blocksize,msample,fs)

[nrow,ncol]=size(Y);

if nrow>ncol % make sure Y is a row matrix
    Y=Y';
    temp=nrow;
    nrow=ncol;
    ncol=temp;
end

nblocks=ncol/blocksize;

[factor,strength,intmax]=CheckSign(Y,blocksize);

%check for anomalies:
intmed=median(intmax')';
goodint=sum(intmax>4*intmed*ones(1,nblocks))==0;

anomCheck=intmax>4*intmed*ones(1,nblocks);

for i=1:nrow
    anomind=find(anomCheck(i,:));
    if ~isempty(anomind)
        anomnum=length(anomind);
        for k=1:anomnum
            ff=find(abs(Y(i,(anomind(k)-1)*blocksize+1:anomind(k)*blocksize))>2*intmed(i));
            Y(i,(anomind(k)-1)*blocksize+ff)=0;
        end
    end
end

if(sum(abs(strength)>25)<2)
    [factor,strength]=CheckSign(Y,blocksize);
end

[sm,si]=sort(abs(strength));
acc_set=si(end-1:end);

accY=(factor(acc_set).*intmed(acc_set)'/sum(intmed(acc_set)))*Y(acc_set,:);
peaks = PeakDetection3(accY,fs,msample,.2,3,goodint);

end