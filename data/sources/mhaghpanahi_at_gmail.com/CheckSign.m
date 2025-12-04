function [y1,y2,y3]=CheckSign(F,blocksize)

% randomly selects intevals of F and check the sign of the peak value.
% final decision on the sign is made based on the majority of the sign of
% the peaks.

% F is assumed to be a row vecor. If F is a matrix, the process is done
% independently on each row. 
margin=0; % if sum(factor)<=margin, then the signal is too noisy to be trustworthy.
[nleads,nsamples]=size(F);
nblocks=nsamples/blocksize;
s=zeros(1,nleads);
strength=zeros(1,nleads);
intmax=zeros(nleads,nblocks);
% interval=1:nblocks;
for j=1:nleads
%     interval=randi(nblocks,1,floor(nblocks/3));
    factor=zeros(1,nblocks);
    for i=1:nblocks
        [absmax ind]=max(abs(F(j,(i-1)*blocksize+1:i*blocksize)));
        factor(i)=sign(F(j,ind+(i-1)*blocksize));
        intmax(j,i)=absmax;
    end
    strength(j)=sum(factor);
    if abs(strength(j))> margin
        s(j)=sign(strength(j));
    end
end

y1=s;
if nargout>1
    y2=strength;
end
if nargout==3
    y3=intmax;
end

    

