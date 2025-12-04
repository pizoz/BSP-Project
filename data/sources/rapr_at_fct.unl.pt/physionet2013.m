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

% ---- check size of ECG ----
if size(ECG,1)>size(ECG,2)
    ECG = ECG';
end

fs      = 1000;             % sampling frequency
N       = size(ECG,1);      % number of abdominal channels
debug   = 0;                % enter debug mode?

% ---- preprocessing ----
[FilteredECG] = preprocessing_v2(ECG,fs);


mqrs=detectqrsmultichannel2(FilteredECG,fs);

% ---- MECG cancellation ----

samplingrate=fs;

qrsL=floor(130*samplingrate/1000);#130ms 
modelorder=91;
olsdata=preparedata2supressmqrs(FilteredECG, mqrs, qrsL, modelorder, fs);


for i=1:N               % run algorithm for each channel
  supressFECG(i,:) =  supressmqrsuseall(FilteredECG,i,mqrs,olsdata,modelorder);
end


% ---- FQRS detection ----

FQRS = cell(N,1);

 for i=1:N
   FQRS{i,1}=detectfqrsupressnew(supressFECG(i,:),fs);
 end




% ---- channel selection based on statistics obtained on the fqrs and mqrs detected
numstats=2;
stats=zeros(N,numstats);
for i=1:N

  if(length(FQRS{i,1})>2 &&  length(mqrs>2))
    aux=[length(FQRS{i,1}),std(diff(FQRS{i,1}))];
    stats(i,:)=aux;
  endif

end

channel=choosechannelnew2(stats);
printf("                        chosen channel is %d\r",channel);fflush(stdout);


fetal_QRSAnn_est    = round(1000*FQRS{channel,1}'/fs);
QT_Interval         = 0;


end





