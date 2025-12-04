function estimated_data=fillMissed(ns,ne,fwdata,bwdata)
% estimates the missing values in the signal using a forward and backward
% AR model
% Inputs:
% ns: first index of missing data interval
% ne: last index of missing data interval
% fwdata: samples before the missing part
% bwdata: samples after the missing part 
% Output: estimated_data
% This code is written based on the estimation algorithm of (Esquef et al., 2006)


N=length(fwdata);
p=floor(N/2);
assert(length(fwdata)==length(bwdata),'predata and postdata sequences must have equal lengthes');
G=ne-ns+1; %length of missing data interval

% Forward data estimation
fwcoeffs=arburg(fwdata,p);
fwexcit=filter(fwcoeffs,1,fwdata);
fwexcit=[fwexcit fwexcit(N:-1:N-G+1)];
fwdata_extend=filter(1,fwcoeffs,fwexcit);


% Backward data estimation
bwcoeffs=arburg(bwdata,p);
bwexcit=filter(bwcoeffs,1,bwdata);
bwexcit=[bwexcit bwexcit(N:-1:N-G+1)];
bwdata_extend=filter(1,bwcoeffs,bwexcit);


% Creating the crossfading window
alpha=2;
U=([ns:ne]-ns)./(G-1);
W=zeros(size(U));
    for i=1:G
        if U(i)<=0.5
            W(i)=1-0.5*(2*U(i))^alpha;
        else
            W(i)=0.5*(2-2*U(i))^alpha;
        end
    end

% Crossfading forward and backkward estimates    
estimated_data=W.*fwdata_extend(end-G+1:end)+(1-W).*bwdata_extend(end-G+1:end);

end