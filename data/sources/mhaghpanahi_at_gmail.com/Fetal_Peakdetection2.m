function [fetal_peaks_round1,Ftest_round1,sample2,find1,fetal_ind_cleaned1,qflag]=Fetal_Peakdetection2(FECG,sample,fs,factorcheat,test_flag)
qflag=1;
blocksize=2000;
[nleads,nsamples]=size(FECG);
nblocks=nsamples/blocksize;
bin=400;


factor=factorcheat;


[Ftest_round1,good_sig1]=AddSig(FECG,factor,sample,fs);
fpeaks1 = PeakDetection3(Ftest_round1,fs,sample,.1,4);
find1=find(fpeaks1==1);
if length(find1)<30
    clear find1 fpeaks1 mn1 r1
    sig_power1=sum(good_sig1);
    [mpow,mpowi]=max(sig_power1);
    test_sample=Ftest_round1((mpowi-1)*blocksize+1:mpowi*blocksize);
    [mt,mti]=max(abs(test_sample(26:blocksize-25)));
    mti=mti+25;
    fpeaks1 = PeakDetection3(Ftest_round1,fs,test_sample(mti-25:mti+24),.1,4);
    find1=find(fpeaks1==1);
end
% % % % assert(length(find1)>30,'sample does not match with Ftest_round1 quite well. Very few peaks found.');
if length(find1) < 50
    qflag=0;
    fetal_peaks_round1=[];
    Ftest_round1=[];
    sample2=0;
    find1=[];
    fetal_ind_cleaned1=[];
    return;
end


%checking for untrustworthy peaks
dfetal=diff(find1);
fetal_ind_test1=dfetal>=600;
fetal_ind_test2= fetal_ind_test1+[0 fetal_ind_test1(1:end-1)]+[fetal_ind_test1(2:end) 0];
fetal_ind_test3=filter([0.5 0 0.5],1,dfetal);
fetal_ind_test3=[fetal_ind_test3(3) fetal_ind_test3(3:end) fetal_ind_test3(end)];
ind_set=find(abs(fetal_ind_test3-dfetal).*~fetal_ind_test2>100)+1;

dd=find(dfetal<max(200,mean(dfetal(~fetal_ind_test1))/1.5));
ind_set=union(ind_set,union(dd,dd+1));


fetal_ind_cleaned1=setdiff(find1,find1(ind_set));

if ~test_flag

    fetal_peaks_round1=FillMissedPeaks3(fetal_ind_cleaned1,fs);

    fetal_peaks_round1=sort(fetal_peaks_round1);

    clear dfetal fetal_ind_test1 fetal_ind_test2 fetal_ind_test3 ind_set dd


    finalpeaks1=zeros(1,60000);
    finalpeaks1(fetal_peaks_round1)=1;
    phase = PhaseCalculation(finalpeaks1);


    nanind=find(sum(good_sig1)==0);
    if length(nanind)>=1
        FF=Ftest_round1(~isnan(Ftest_round1));
        pp=phase(~isnan(Ftest_round1));
    %     pp=pphase(1:(nanind(1)-1)*blocksize);
    else
        FF=Ftest_round1;
        pp=phase;
    end

    Fmean = MeanECGExtraction(FF,pp,bin,1);

    [Fsigmax FI] = max(abs(Fmean));
    sample2=Fmean(FI-30:FI+29);
%     sample2=Fmean(FI-30:FI+99);
else
    fetal_peaks_round1=0;
    FF=Ftest_round1(~isnan(Ftest_round1));
    ll=PeakDetection3(FF,fs,sample,0.1,4);
    pp=PhaseCalculation(ll);
    Fmean=MeanECGExtraction(FF,pp,bin,1);
    [dummy, FI] = max(abs(Fmean));
    sample2=Fmean(FI-30:FI+29);
%     sample2=Fmean(FI-30:FI+99);
end


end




function [Ftest, good_sig]=AddSig(FECG,factor,sample,fs)

blocksize=2000;
leadset=factor~=0;
[nleads,nsamples]=size(FECG);
nblocks=nsamples/blocksize;
pinterval=25:41;
interval=21:38;
F=0:150;
Ftest=zeros(1,nsamples);
maxFtest=zeros(1,nblocks);
psd_peak=zeros(1,nblocks);
good_sig=zeros(4,nblocks);


for j=1:nblocks
    for i=1:4    
        [S(:,i),F,T,P(:,i)] = spectrogram(FECG(i,(j-1)*blocksize+1:j*blocksize),blocksize,[],F,fs);
    end

    peakf=PsdPeak2(P,pinterval,leadset);
    psd_peak(j)=peakf;
    Fqindex=logical(zeros(1,4));
    Tqindex=logical(zeros(1,4));
    qindex=logical(zeros(1,4));
    Pqindex=logical(zeros(1,4));

    for i=1:4
%         [dummy FI] = max(abs(FECG(i,(j-1)*blocksize+1:j*blocksize)));
        if leadset(i)
            peaktest=PeakDetection3(FECG(i,(j-1)*blocksize+1:j*blocksize),fs,factor(i)*sample,.01,4);
            peakdiff=diff(find(peaktest==1));
            numpeak=length(peakdiff);
            Tqindex(i)=sum(peakdiff>250 & peakdiff<600) >= 2/3 * numpeak & numpeak>=2 & numpeak <=5 ;
            Pqindex(i)= all(peakdiff>250 & peakdiff<600) & sum(peakdiff)>=0.65*blocksize & all(abs(peakdiff-sum(peakdiff)/numpeak)<100);%& (numpeak==3 | numpeak==4);

        end
    end
    if peakf>0
        Larea=sum(P(1:interval(1)-1,:));
        Rarea=sum(P(interval(end)+1:end,:));
        area_test= Larea./Rarea>=1.2 | sum(P(interval,:))./sum(P)>0.4;% & Rarea./sum(P)<0.275);
        peak_testL=sum(P(peakf-5:peakf+5,:))./Larea>=0.21;
        peak_testR=sum(P(peakf-5:peakf+5,:))./Rarea>=0.21;
        Fqindex=area_test & peak_testL & peak_testR;
        set=setdiff(peakf-2:peakf+2,peakf);
    else
        Larea(Pqindex)=sum(P(1:interval(1)-1,Pqindex));
        Rarea(Pqindex)=sum(P(interval(end)+1:end,Pqindex));
        area_test(Pqindex)= Larea(Pqindex)./Rarea(Pqindex)>=1.2 | sum(P(interval,Pqindex))./sum(P(:,Pqindex))>0.4;

        peak_testL(Pqindex)=sum(P(interval,Pqindex))./Larea(Pqindex)>=0.21;
        peak_testR(Pqindex)=sum(P(interval,Pqindex))./Rarea(Pqindex)>=0.21;
        Fqindex(Pqindex)=area_test(Pqindex) & peak_testL(Pqindex) & peak_testR(Pqindex);
    end
    if peakf>0
        qindex=Tqindex & Fqindex | Pqindex;
    else
        qindex(Pqindex)=Fqindex(Pqindex);
    end
  
    if any(qindex)
        Ftest((j-1)*blocksize+1:j*blocksize)= factor(qindex)*FECG(qindex,(j-1)*blocksize+1:j*blocksize)/sum(qindex);
        good_sig(:,j)=qindex;
        maxFtest(j)=max(Ftest((j-1)*blocksize+1:j*blocksize));
    else
        Ftest((j-1)*blocksize+1:j*blocksize)=NaN;
    end
    
    clear S T P Larea Rarea;
end
end





    
