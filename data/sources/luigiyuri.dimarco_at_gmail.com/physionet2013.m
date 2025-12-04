function [fetal_QRSAnn_est,QT_Interval] = physionet2013(tm,ECG)
%
% Uses: 
%       - Signal Processing Toolbox
%       - Statistics Toolbox
%       - m-modules: 
%                   BuildChain.m 
%                   FindEMG.m 
%                   MyAlignQRS.m
%                   MyBeatDect.m
%                   MyBuildSeq.m
%                   MyFBeatDect.m 
%                   MyRRTemplate.m
%
%

% ---- check size of ECG ----
if size(ECG,2)>size(ECG,1)
    ECG = ECG';
end

% impute missing data
for j=1:size(ECG,1)
    for k=1:size(ECG,2)
        if(isnan(ECG(j, k)))
            if(j==1)
                ECG(j, k) = 0;
            else
                ECG(j, k) = ECG(j-1, k);
            end
        end
    end
end

Fs      = 1000;             % sampling frequency

% notch filter - 50Hz
for k=1:size(ECG,2)
    fP = FindEMG(ECG(:,k), Fs);
    if(fP>0)
        wo    = 50/(Fs/2);  bw = wo/45;
        [b,a] = iirnotch(wo,bw);
        ECG(:,k)   = filtfilt(b,a,ECG(:,k));
    end
end

% HP filter
[b,a] = butter(2, 8/(Fs/2), 'high');
ECGh  = filtfilt(b,a,ECG);

% QRS dect in abdominal ECG 
mqrs = MyBeatDect(ECG, Fs);

% correct QRS series, subtracrt maternal ECG
rECG = MyRRTemplate(ECGh, mqrs, Fs, []);

% fQRS dect 
fqrs = MyFBeatDect(rECG, Fs, [], mqrs, 0);

fetal_QRSAnn_est    = round(1000*fqrs'/Fs);
QT_Interval         = 0;
