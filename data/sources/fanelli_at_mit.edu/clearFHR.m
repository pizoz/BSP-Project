%% clean RR

function [cleanQrsMat]=clearFHR(F_AVG_TIME,qrsMat,SAMPL_FREQ)

%matrix of coeffient
coeff=[0 0.2 0.4 0.6 0.8 1];

%expected average frequency

F_AVG=F_AVG_TIME/SAMPL_FREQ;

%parameters
FUTURE_BEATS=7; %number of future beats to consider for the correction

%Inizialization;
i=2;
j=2;
qrsAvg=F_AVG;
cleanQrsMat=qrsMat(1,1);
beats=zeros(7,3);



    while j<=length(qrsMat),

        beats(:,1:2)=qrsMat(j:j+FUTURE_BEATS-1,:);
        beats(:,3)=(1:FUTURE_BEATS);

        %adaptation process for the first 6 beats
        if i<8,
        qrsAvg=coeff(8-i)*F_AVG + coeff(i-1)*qrsAvg;
        end

        %distance between last correctly detected beat and the following N=FUTURE_BEATS 
        %beats    
        dist=beats(:,1)-cleanQrsMat(i-1);

        %distance between interbeat distances and expected interbeat distance
        ranking=abs(dist-qrsAvg);

        %detection of beats close to the avg value
        okBeats=find(ranking < 0.15*qrsAvg);

        %if none sutisfy the previous condition, the best one is selected
        if isempty(okBeats),        
           [~,okBeats]= min(ranking);        
        end

        %identification of the best beat according to the cross correlation
        %value
        data=beats(okBeats,:);    
        [~,b]=max(data(:,2));
        index=data(b,3);    

        cleanQrsMat(i)=beats(index,1);

        %construction of the rr series;
        rr(i-1)=cleanQrsMat(i)-cleanQrsMat(i-1);


        j=j+index;

        %update of the average interbeat distance
        if i < 12,
            qrsAvg=(qrsAvg+(rr(i-1)/(i-1)))*((i-1)/i);
            else    
            qrsAvg=mean(rr(i-10:i-1));    
        end

        i=i+1; 

        if j+FUTURE_BEATS - 1 > length(qrsMat),          
            FUTURE_BEATS=length(qrsMat)-j+1;
            beats=zeros(FUTURE_BEATS,3);    
        end

    end %while


end %function







    
    