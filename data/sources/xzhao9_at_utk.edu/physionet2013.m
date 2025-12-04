function [fetal_QRSAnn_est,QT_Interval] = physionet2013(tm,ECG)
% Template algorithm for Physionet/CinC competition 2013. This function can
% be used for events 1 and 2. Participants are free to modify any
% components of the code. However the function prototype must stay the
% same:
%
% [fetal_QRSAnn_est,QT_Interval] = physionet2013(tm,ECG) where the inputs and outputs are specified
% below.
%
% inputs:
%   ECG: 4x60000 (4 channels and 1min of signal at 1000Hz) matrix of
%   abdominal ECG channels.
%   tm : Nx1 vector of time in milliseconds
% output:
%   FQRS: FQRS markers in seconds. Each marker indicates the position of one
%   of the FQRS detected by the algorithm.
%   QT_Interval:   1x1 estimated fetal QT duration (enter NaN or 0 if you do wish to calculate)
%
%
% Author: Joachim Behar - IPMG Oxford (joachim.behar@eng.ox.ac.uk)
% updated: March 3, 2013 Ikaro Silva
% update: March 22, 2013 Xiaopeng Zhao
%

% ---- check size of ECG ----
if size(ECG,2)>size(ECG,1)
    ECG = ECG';
end

fs      = 1000;             % sampling frequency
N_chn   = size(ECG,2);      % number of abdominal channels
debug   = 0;                % enter debug mode?
chn     = 2;                % use channel 2 to detect fQRS

% ---- preprocessing ----
gflag=1;nflag=1;bflag=1;
[X] = preprocessing(tm,ECG,gflag,nflag,bflag);

% ---- find maternal R peaks
wtScale_qrs=2:4;
qrs_coeff = cwt(X(:,chn),wtScale_qrs,'mexh');
mEnergy = nanmean(qrs_coeff.^2,1);
mPeaks = PeakDetection(mEnergy,1/fs,1);

% ---- MECG cancellation ----
fEnergy = MECGcancellation(mPeaks,mEnergy,fs);

% ---- find fetal R peaks ----
fPeaks = PeakDetection(fEnergy,2/fs,1);

% ---- fill in gaps in R peaks ----
% the minimum normal heart rate is 120 bpm, i.e. 500 ms per beat
fRR=diff(fPeaks);
id= fRR<=500;
avgRR=mean(fRR(id)); %average fetal RR interval excluding gaps

% fill the large intervals with average heartbeats
fold=[0,fPeaks,size(X,1)];
dRR=diff(fold);
fGaps=[];%gaps will be filled
for i=1:length(dRR)
    for j=1:round(dRR(i)/avgRR)-1
        fGaps=[fGaps,fold(i)+round(j*avgRR)];
    end
end

% ---- output of fetal QRS ----
fetal_QRSAnn_est=sort([fPeaks,fGaps]);
QT_Interval         = 0;

end



