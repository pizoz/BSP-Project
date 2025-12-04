function [FECG,fetal_peaks_round1,fetal_ind_cleaned1,psample] = PCMPeakDetect(x,peaks,fsample,goodint,fs,blocksize)

[nrow,ncol]=size(x);

if nrow > ncol % make sure x is a row vector
    x=x';
    temp=nrow;
    nrow=ncol;
    ncol=temp;
end

ind1=find(peaks==1);
for i=1:length(ind1)
    x(max(1,ind1(i)-60):min(ind1(i)+60,60000))=0;
end

xL=x-LPFilter(x,15/1000);
xB=BPFilter(x,15/1000,45/1000);



choice=[x' xL' xB'];
F_order = FSelect(choice,blocksize);
    
switch F_order
    case 1
        FECG = x;
    case 2
        FECG = xL;
    case 3
        FECG = xB;
end




[fpeaks1,mn1,r1] = PeakDetection3(FECG,fs,fsample,.1,4,goodint);
[fpeaks2,mn2,r2] = PeakDetection3(FECG,fs,-fsample,.1,4,goodint);

if sum(fpeaks1)>=sum(fpeaks2)
    find1=find(fpeaks1==1);
    pp=PhaseCalculation(fpeaks1);
    sign=1;
else
    find1=find(fpeaks2==1);
    pp=PhaseCalculation(fpeaks2);
    sign=-1;
end

if length(find1) < 70
    fetal_peaks_round1=[];
    psample=0;
    fetal_ind_cleaned1=[];
    return;
end



bin=250;
Fmean=MeanECGExtraction(FECG,pp,bin,1);
[dummy, FI] = max(abs(Fmean));
psample=sign*Fmean(max(0,FI-30):min(bin,FI+29));



dfetal=diff(find1);
fetal_ind_test1=dfetal>=600;
fetal_ind_test2= fetal_ind_test1+[0 fetal_ind_test1(1:end-1)]+[fetal_ind_test1(2:end) 0];
fetal_ind_test3=filter([0.5 0 0.5],1,dfetal);
fetal_ind_test3=[fetal_ind_test3(3) fetal_ind_test3(3:end) fetal_ind_test3(end)];
ind_set=find(abs(fetal_ind_test3-dfetal).*~fetal_ind_test2>100)+1;

dd=find(dfetal<max(200,mean(dfetal(~fetal_ind_test1))/1.5));
ind_set=union(ind_set,union(dd,dd+1));


fetal_ind_cleaned1=setdiff(find1,find1(ind_set));
fetal_peaks_round1=FillMissedPeaks3(fetal_ind_cleaned1,fs);

end
