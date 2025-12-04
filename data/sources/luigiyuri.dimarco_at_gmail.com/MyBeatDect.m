function [qrs, ECGp, pc] = MyBeatDect(ECG, Fs)
%
% function: computes beat detection
% 
% IN: 
%     ECG: raw data Nx4 matrix, N=Fs*60s
%     Fs:  sample rate
%
% OUT:
%     qrs:  detected beat time series [smpls] 
%           Nx4 matrix. Columns contain detected beats, padded with 0 to fill to N  
%     ECGp: PCA-ECG is Nx4 matrix
%     pc:   selected PC
%

% LP filter
[b,a] = butter(4, 25/(Fs/2));
ECGl  = filtfilt(b,a,ECG);
% HP filter
[b,a] = butter(4, 5/(Fs/2), 'high');
ECG  = filtfilt(b,a,ECGl);

N = 32+1;
b = ones(1,N)/N; a = 1;
Le  = size(ECG,1); 


for k=1:size(ECG,2)
    D = diff(ECG(:,k)).^2;
    Sk = filtfilt(b,a,D);
    S(:,k) = Sk/std(Sk);
end

Sx = sqrt( S(:,1).^2+S(:,2).^2+S(:,3).^2+S(:,4).^2 );

[pks, locs] = findpeaks(Sx);

DW  = 1*Fs; qrs = [];
for q=1:numel(locs)
    i1  = max(1, locs(q)-DW);
    i2  = min(i1+2*DW, Le);
    pkM = max(pks(locs>=i1 & locs<=i2));
    if( pks(q)>(0.5*pkM) )
        qrs = [qrs locs(q)];
    end
end
