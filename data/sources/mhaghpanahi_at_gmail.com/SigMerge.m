function [y, set] = SigMerge(FECG,fsample,goodint,fs)

set=SetSelect(FECG,fs);

j=1;

factorcheat=[0 0 0 0];
test_flag=1;
sig_sign=[0 1 -1];

factorcheat(set(j))=1;
[dummy,dummy,s0,dummy,findc0,qflag0]=Fetal_Peakdetection2(FECG,fsample,fs,factorcheat,test_flag);
factorcheat(set(j))= qflag0*sign(sum(fsample.*s0));

if factorcheat(set(j))==-1
    [dummy,dummy,s0,dummy,findc0,qflag0]=Fetal_Peakdetection2(FECG,fsample,fs,factorcheat,test_flag);
end


if qflag0==0
    fpeaks1 = PeakDetection3(FECG(set(j),:),fs,fsample,.1,4,goodint);
    fpeaks2 = PeakDetection3(FECG(set(j),:),fs,-fsample,.1,4,goodint);

    if sum(fpeaks1)>=sum(fpeaks2)
        find1=find(fpeaks1==1);
        pp=PhaseCalculation(fpeaks1);
        sgn=1;
    else
        find1=find(fpeaks2==1);
        pp=PhaseCalculation(fpeaks2);
        sgn=-1;
    end
    dfetal=diff(find1);
    fetal_ind_test1=dfetal>=600;
    fetal_ind_test2= fetal_ind_test1+[0 fetal_ind_test1(1:end-1)]+[fetal_ind_test1(2:end) 0];
    fetal_ind_test3=filter([0.5 0 0.5],1,dfetal);
    fetal_ind_test3=[fetal_ind_test3(3) fetal_ind_test3(3:end) fetal_ind_test3(end)];
    ind_set=find(abs(fetal_ind_test3-dfetal).*~fetal_ind_test2>100)+1;
    dd=find(dfetal<max(200,mean(dfetal(~fetal_ind_test1))/1.5));
    ind_set=union(ind_set,union(dd,dd+1));
    findc0=setdiff(find1,find1(ind_set));
    
    
    bin=250;
    Fmean=MeanECGExtraction(FECG(set(j),:),pp,bin,1);
    [dummy, FI] = max(abs(Fmean));
    s0=sgn*Fmean(max(0,FI-30):min(bin,FI+29));
    factorcheat(set(j))=sgn;
    
    if length(findc0)<60
        qflag0=0;
    else
        qflag0=1;
    end
    
end

corrcoef=length(findc0);



j=j+1;
if qflag0
    
    while j<=4
        factorcheat(set(j))=1;
        [dummy,dummy,s1,dummy,findc1]=Fetal_Peakdetection2(FECG,fsample,fs,factorcheat,test_flag);
        factorcheat(set(j))=-1;
        [dummy,dummy,s2,dummy,findc2]=Fetal_Peakdetection2(FECG,fsample,fs,factorcheat,test_flag);
    %     [total_peak, ind]=max([total_peak, length(findc1), length(findc2)]);
        [corrcoef,ind]=max([corrcoef (sum(s0.*s1)/(norm(s0)*norm(s1)))*length(findc1) (sum(s0.*s2)/(norm(s0)*norm(s2)))*length(findc2)]);
        factorcheat(set(j))=sig_sign(ind);
        j=j+1;
    end

end

y=factorcheat;

end

% % % % function S = SetSelect(FECG,fs,checkno)
% % % % 
% % % % [nrow,ncol]=size(FECG);
% % % % 
% % % % 
% % % % blocksize=2000;
% % % % nblocks=ncol/blocksize;
% % % % range=randi(nblocks,1,checkno);
% % % % sub_set=[];
% % % % 
% % % % 
% % % % F=0:150;
% % % % interval=21:38;
% % % % maxInd= zeros(1,checkno);
% % % % 
% % % % for j=1:checkno
% % % %     for i=1:nrow    
% % % %         [S(:,i),F,T,P(:,i)] = spectrogram(FECG(i,(range(j)-1)*blocksize+1:range(j)*blocksize),blocksize,[],F,fs);
% % % %     end
% % % % 
% % % %     p_ratio =  sum(P(interval,:))./sum(P);
% % % %     
% % % %     
% % % %     k_meas=kurtosis(FECG(:,(range(j)-1)*blocksize+1:range(j)*blocksize)');
% % % %     
% % % %     [dummy, maxInd(j)]=max(k_meas.*p_ratio);
% % % % 
% % % %     clear S T P;
% % % % end
% % % % 
% % % % compare= hist(maxInd,1:nrow);
% % % % [csort, set]=sort(compare,'descend');
% % % % 
% % % % if sum(csort==0)<2
% % % %     S = set;
% % % %     return;
% % % % else
% % % %     sub_set = SetSelect(FECG(set(csort==0),:),fs,checkno);
% % % % end
% % % % 
% % % % w=set(csort==0);
% % % % S = [set(csort~=0) w(sub_set)];
% % % % 
% % % % end




