function residual = MECGcancellation(peaks,ECG,fs)
% MECG cancellation using blanking
%
% inputs:
%   fs: sampling frequency
%   ECG: matrix of abdominal ECG channels.
%   peaks: MQRS markers in seconds. Each marker corresponds to the
%   position of a MQRS.
%
% output:
%   residual: residual containing the FECG.
% ---- constants ----
Qstart              = floor(0.1*fs);
Sstop               = floor(0.1*fs);
N                   = length(peaks);                  % number of MECG QRS
residual            = ECG;

% ---- MECG cancellation ----
for i=1:N
    left=min(peaks(i)-1,Qstart);
    right=min(length(ECG)-peaks(i),Sstop);
    residual(peaks(i)-left:peaks(i)+right) = 0;
end

% % the starting and end have "edge" effect
% residual(1:10)=0;
% residual(end-9:end)=0;

end
