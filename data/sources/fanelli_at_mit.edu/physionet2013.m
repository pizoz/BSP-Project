%% Algorithm for FECG extraction from abdominal recordings
%The algorithm extracts FECG traces from abdominal recordings from pregnancies.


%Andrea Fanelli
%May 2013

%MODIFICATO


%% Definition of parameters

function [cleanedBestRR,QT_Interval] = physionet2013(tm,ECG)


load filt3


SMOOTH_PARAM=50;
SIG_NUMBER=size(ECG,2);

UPPER_THRESH=0.4;
LOWER_THRESH=0.2;

			
SAMPL_FREQ=1000;


START_POINT=1;      

MAT_QRS_AVG=10; %number of maternal QRS to be averaged

MAT_QRS_TIME=0.12;	 %duration of maternal QRS in sec
MAT_QRS_SAMPLE=floor(SAMPL_FREQ*MAT_QRS_TIME);		 %duration of maternal QRS in samples


%definition of left and right maternal QRS windows
if rem(MAT_QRS_SAMPLE,2)==0,
    L_MQRS=MAT_QRS_SAMPLE/2-1;
    R_MQRS=MAT_QRS_SAMPLE/2;
else
    L_MQRS=(MAT_QRS_SAMPLE-1)/2;
    R_MQRS=(MAT_QRS_SAMPLE-1)/2;    
end


FET_QRS_TIME=0.06; %duration of fetal QRS in sec
FET_QRS_SAMPLE=floor(SAMPL_FREQ*FET_QRS_TIME); %duration of fetal QRS in samples

%definition of left and right fetal QRS windows
if rem(FET_QRS_SAMPLE,2)==0,
    L_FQRS=FET_QRS_SAMPLE/2-1;
    R_FQRS=FET_QRS_SAMPLE/2;
else
    L_FQRS=(FET_QRS_SAMPLE-1)/2;
    R_FQRS=(FET_QRS_SAMPLE-1)/2;    
end
 	
DATA_LENGTH=length(ECG);

if DATA_LENGTH < 3500,
    error('Data is too short to be processed. It should be at least 3500 samples long.');
end


%% data correction

for l=1:size(ECG,2)
    wrongValues=find(isnan(ECG(:,l)));
    ECG(wrongValues,l)=0;
end

%% Preprocessing

postSmooth=zeros(DATA_LENGTH ,SIG_NUMBER);
noMean=zeros(DATA_LENGTH,SIG_NUMBER);
noStd=zeros(DATA_LENGTH,SIG_NUMBER);


% filtering

abdomRecB=filtfilt(Hbs.numerator,1, ECG);
abdomRec=filtfilt(Hbp3.numerator,1, abdomRecB);




for h=1:SIG_NUMBER,

    %mean subtraction
    noMean(:,h)=abdomRec(:,h)-mean(abdomRec(:,h));
    %std subtraction
    noStd(:,h)=noMean(:,h)./(sqrt(sum(noMean(:,h).^2)));

end


%% Select channel for maternal Qrs detection

matMatrix=noStd;
abdQuality=sum(abs(diff(matMatrix,1)));
[~,ch]=min(abdQuality);

matSignal=matMatrix(:,ch);

%% Mat QRS detection

qrsMat=[];

i=1;

while i - 1 + SAMPL_FREQ + R_MQRS <= DATA_LENGTH,

    %detection of Maternal QRS template for each window
    
    sigWindow=matSignal(i:i+SAMPL_FREQ-1);
    [~,time]=max(abs(sigWindow));

    if time - L_MQRS < 1 && i==1, 
           winStart=1;
    else
           winStart=time - L_MQRS;
    end
    
    mqrsTempl=matSignal(i -1 +winStart:i - 1 + time + R_MQRS);    
    
     
    %avoid mistakes in the last window
            
    if (i + SAMPL_FREQ + MAT_QRS_SAMPLE <= DATA_LENGTH),

        
        %compute of cross correlation for each window      
        crossTempl=conv(matSignal(i:i+SAMPL_FREQ+MAT_QRS_SAMPLE),rot90(rot90(mqrsTempl)));
        crossTempl=crossTempl(MAT_QRS_SAMPLE:end-MAT_QRS_SAMPLE-1);
        
        %normalize cross correlation signal between 0 and 1        
        crossTempl= crossTempl ./ sum(mqrsTempl.^2);
   
        %detection of all the peaks in the signal
       	overThresh=0;
        actualMax=0;

        for k=1:length(crossTempl),       
            if (overThresh==0 && (crossTempl(k)) > UPPER_THRESH),
                overThresh=1;
            end
            if (overThresh==1),
                 if (crossTempl(k) > actualMax), 
                     actualMax = crossTempl(k);
                     peakTime=k;
                 end

            end
    
    
            if (overThresh==1 && crossTempl(k) <  LOWER_THRESH), 
                overThresh=0;
                
                [~,realTime]=max(abs(matSignal(i + peakTime: i + peakTime + MAT_QRS_SAMPLE)));
                
                %if the interbeat distance is > 1 sec, the algorithm identify the peak between the actual beat and the previous detected beat 
                if (isempty(qrsMat)==0) && (i-1 + peakTime + realTime - qrsMat(end,1)) > SAMPL_FREQ,                    
                    [~,timePast]=max(abs(matSignal(qrsMat(end,1) + MAT_QRS_SAMPLE : i + peakTime)));
                    qrsMat=[qrsMat; [qrsMat(end,1) + MAT_QRS_SAMPLE + timePast, -1]];
                end
                       
                %simple control to avoid double identification
                if ((isempty(qrsMat)==1) || (i-1 + peakTime + realTime - qrsMat(end,1) > 5))
                    qrsMat=[qrsMat; [i-1 + peakTime + realTime, actualMax]];
                end
                  
                actualMax=0;


            end
        

        end
        

    end

i=i+SAMPL_FREQ-1;    
end

F_AVG_MAT=mean(diff(qrsMat));

%% Clean extracted mecg
[mQrs]=clearFHR(F_AVG_MAT(:,1)/1000,qrsMat,1000); 


% figure, plot(matSignal);
% hold on
% plot(qrsMat(:,1),matSignal(qrsMat(:,1)),'.g');
% plot(mQrs,matSignal(mQrs),'.r');

%% Smoothing filter to increase fetal to maternal SNR

abdomRec=filtfilt(Hlp.numerator,1, abdomRecB);
for h=1:SIG_NUMBER,

    %mean subtraction
    noMean(:,h)=abdomRec(:,h)-mean(abdomRec(:,h));
    %std subtraction
    noStd(:,h)=noMean(:,h)./(sqrt(sum(noMean(:,h).^2)));

end

for h=1:SIG_NUMBER,
postSmooth(:,h)=noStd(:,h)-smooth(noStd(:,h),SMOOTH_PARAM);
end

%% Mat QRS removal

fEcgMatrix= zeros(DATA_LENGTH,SIG_NUMBER); 
mQrsTempl=zeros(MAT_QRS_SAMPLE,1);

for z=1:SIG_NUMBER,
   
    fEcg=postSmooth(:,z);
    
    if exist('mQrsMatrix','var')==1,
        clear mQrsMatrix
    end
    
    
    for j=1:length(mQrs),
       
        if (j~=1 || mQrs(j) - L_MQRS >= 1),            
            
            mBeat=fEcg(mQrs(j)-L_MQRS:mQrs(j)+R_MQRS);
            
            if exist('mQrsMatrix','var')==0,
                
                mQrsMatrix=mBeat;
                mQrsTempl=mBeat;
               
            else
                
                rmsDiff=sqrt(sum((mQrsTempl-mBeat).^2)/MAT_QRS_SAMPLE);   
                mQrsMatrix=cat(2,mQrsMatrix,mBeat); 
                   
                if size(mQrsMatrix,2)> MAT_QRS_AVG, 
                    
                       mQrsMatrix=mQrsMatrix(:,2:MAT_QRS_AVG+1);
                end
                       mQrsTempl=mean(mQrsMatrix,2);

                       a=((mQrsTempl'*mQrsTempl)^(-1))*(mQrsTempl')* mBeat;
                       mBeatCap=mQrsTempl*a;

                       fEcg(mQrs(j)-L_MQRS:mQrs(j)+R_MQRS)=mBeat-mBeatCap;

                   
            end
  
        else
            
           fEcg(1:mQrs(j)+R_MQRS)=0; 
            
        end % j != 1
        
    end % j > length(mQrs) 
    
    fEcgMatrix(:,z)=fEcg;
    
end % z > SIG_NUMBER
% 
% for ll=1:SIG_NUMBER,
%     figure, plot(fEcgMatrix(:,ll));
% end



%% Detezione QRS fetali

for z=1:SIG_NUMBER,

    depFecg=fEcgMatrix(:,z)-mean(fEcgMatrix(:,z));
    normFecg=depFecg/sqrt(sum(depFecg.^2)); 
    
    fQrs=[];
    
    for i=1:SAMPL_FREQ/2:DATA_LENGTH - SAMPL_FREQ/2 -1,

        %detection of Fetal QRS template for each window

        sigWindow=normFecg(i:i+ SAMPL_FREQ/2 -1);
        [peak,time]=max(abs(sigWindow));

        if time - L_FQRS < 1 && i==1, 
               winStart=1;
        else
               winStart=time - L_FQRS;
        end

        fqrsTempl=normFecg(i -1 +winStart:i - 1 + time + R_FQRS);    


        %avoid mistakes in the last window

        if (i + SAMPL_FREQ/2 + FET_QRS_SAMPLE <= DATA_LENGTH),


            %compute of cross correlation for each window      
            crossTempl=conv(normFecg(i:i+SAMPL_FREQ/2+FET_QRS_SAMPLE),rot90(rot90(fqrsTempl)));
            crossTempl=crossTempl(FET_QRS_SAMPLE:end-FET_QRS_SAMPLE-1);

            %normalize cross correlation signal between 0 and 1        
            crossTempl= crossTempl ./ sum(fqrsTempl.^2);

            %detection of all the peaks in the signal
            overThresh=0;
            actualMax=0;

            for k=1:length(crossTempl),

                if (overThresh==0 && (crossTempl(k)) > UPPER_THRESH),
                    overThresh=1;
                end

                if (overThresh==1),

                     if (crossTempl(k) > actualMax), 

                         actualMax = crossTempl(k);
                         peakTime=k;

                     end

                end


                if (overThresh==1 && crossTempl(k) <  LOWER_THRESH), 

                    overThresh=0;
                    [realPeak,realTime]=max(abs(normFecg(i + peakTime: i + peakTime + FET_QRS_SAMPLE)));
                    
                    %if the interbeat distance is > 1 sec, the algorithm identify the peak between the actual beat and the previous detected beat 
%                     if (isempty(fQrs)==0) && (i-1 + peakTime + realTime - fQrs(end,1)) > 500, 
% 
%                         [pastPeak,timePast]=max(abs(matSignal(fQrs(end,1) + FET_QRS_SAMPLE : i + peakTime)));
%                         fQrs=[fQrs; [fQrs(end,1) + FET_QRS_SAMPLE + timePast, -1]];
%    
%                     end
                    
                    
                    %simple control to avoid double identification
                    if ((isempty(fQrs)==1) || (i-1 + peakTime + realTime - fQrs(end,1) > 5))
                        fQrs=[fQrs; [i-1 + peakTime + realTime, actualMax]];
                    end
                    

                    actualMax=0;


                end


            end


        end


    end
    
    eval(['qrsFet',num2str(z),'=fQrs;']);

end




%% Clean extracted fecg

quality=zeros(4,1);

for k=1:SIG_NUMBER,
    if eval(['isempty(qrsFet',num2str(k),')==0,']);
         eval(['cleanQrsFet',num2str(k),'=clearFHR(mean(diff(qrsFet',num2str(k),'(:,1)))/1000,qrsFet',num2str(k),',1000);']);   
         eval(['quality(k)=mean(abs(diff(diff(cleanQrsFet',num2str(k),'))));']);
    else  
        quality(k)=NaN;
    end
end

%% Identify HRV series of best quality

[~,bestQuality]=min(quality);

eval(['bestRR=qrsFet',num2str(bestQuality),'(:,1);']);
eval(['cleanedBestRR=cleanQrsFet',num2str(bestQuality),';']);

%% Visualization

% figure, plot(ECG(:,bestQuality));
% hold on
% % plot(bestRR,ECG(bestRR,bestQuality),'.g');
% plot(cleanedBestRR,ECG(cleanedBestRR,bestQuality),'.r');

%% 

QT_Interval=NaN;




