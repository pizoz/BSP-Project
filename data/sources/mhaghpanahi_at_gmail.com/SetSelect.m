function S = SetSelect(FECG,fs)

[nrow,ncol]=size(FECG);


blocksize=2000;
nblocks=ncol/blocksize;
subset=[];

% F=0:150;
% interval=21:38;
maxInd= zeros(1,nblocks);

for j=1:nblocks
%     for i=1:nrow    
%         [S(:,i),F,T,P(:,i)] = spectrogram(FECG(i,(j-1)*blocksize+1:j*blocksize),blocksize,[],F,fs);
%     end
% 
%     p_ratio =  sum(P(interval,:))./sum(P);
    
    
    k_meas=kurtosis(FECG(:,(j-1)*blocksize+1:j*blocksize)');
    
%     [dummy, maxInd(j)]=max(k_meas.*p_ratio);
    [dummy, maxInd(j)]=max(k_meas);

    clear S T P;
end

compare= hist(maxInd,1:nrow);
[csort, set]=sort(compare,'descend');

if sum(csort==0)<2
    S = set;
    return;
else
    sub_set = SetSelect(FECG(set(csort==0),:),fs);
end

w=set(csort==0);
S = [set(csort~=0) w(sub_set)];

end