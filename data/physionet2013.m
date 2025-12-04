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
% Last updated: March 3, 2013 Ikaro Silva
%               April 15, 2013 Joachim Behar

try
% ---- check size of ECG ----
if size(ECG,2)>size(ECG,1)
    ECG = ECG';
end

fs      = 1000;             % sampling frequency
N       = size(ECG,2);      % number of abdominal channels
debug   = 0;                % enter debug mode?

% ---- preprocessing ----
[FilteredECG] = preprocessing(ECG,fs);

MQRS = PeakDetection(FilteredECG(:,1),1/fs);

% ---- MECG cancellation ----
for i=1:N               % run algorithm for each channel
    FECG(:,i) = MECGcancellation(MQRS,FilteredECG(:,i)',fs,20);
end

% ---- Post processing of the residual ----
FECG = ResidualPostProcessing(FECG);

% ---- channel selection or combination prior FQRS detection ----
[SelectedResidual,ChannelNb] = ChannelSelectionOrCombination(FECG);

% ---- FQRS detection ----
FQRS    = PeakDetection(SelectedResidual,2/fs);

if debug
    m = length(SelectedResidual);
    plot(tm,FilteredECG(:,1),'LineWidth',2);
    hold on, plot(tm,SelectedResidual,'r',...
        tm(FQRS),SelectedResidual(FQRS),'+r',...
        tm,FilteredECG(:,ChannelNb)-SelectedResidual,'--k','LineWidth',2);
    
    title('Extracted FECG and detected FQRS');
    xlabel('Time (sec)'); ylabel('Amplitude (NU)');
end

fetal_QRSAnn_est    = round(1000*FQRS'/fs);
QT_Interval         = 0;

catch ME
    ME
    fetal_QRSAnn_est    = [];
    QT_Interval         = [];
end

end

function [FilteredECG] = preprocessing(ECG,fs)
% ---- preprocess the data ----
FilteredECG = ECG;

end

function residual = MECGcancellation(peaks,ECG,fs,nbCycles)
% MECG cancellation algorithm inspired from [1].
%
% inputs:
%   fs: sampling frequency
%   nbCycles: number of cycles on which to build the mean MECG template
%   ECG: matrix of abdominal ECG channels.
%   peaks: MQRS markers in seconds. Each marker corresponds to the
%   position of a MQRS.
%
% output:
%   residual: residual containing the FECG.
%
% Author: Joachim Behar - IPMG Oxford (joachim.behar@eng.ox.ac.uk)
% Last updated: 03_02_2013
%
% [1] Martens S et al. A robust fetal ECG detection method for
% abdominal recordings. Physiol. Meas. (2007) 28(4) 373�388

% ---- constants ----
r                   = nbCycles;
ECG_last_r_cycles   = zeros(0.7*fs,r);
Pstart              = 0.25*fs-1;
Tstop               = 0.45*fs;
N                   = length(peaks);                    % number of MECG QRS
ECG_temp            = zeros(1,length(ECG));

% ---- ECG template ----
for i=1:r
    peak_nb = peaks(i+1);   % +1 to unsure full cycles
    ECG_last_r_cycles(:,i) = ECG(peak_nb-Pstart:peak_nb+Tstop)';
end
ECG_mean = mean(ECG_last_r_cycles,2);

% ---- MECG cancellation ----
for i=1:N
    if peaks(i)>Pstart && length(ECG)-peaks(i)>Tstop
        M  = zeros (0.7*fs,3);
        M(1:0.2*fs,1)           = ECG_mean(1:Pstart-0.05*fs+1);
        M(0.2*fs+1:0.3*fs,2)    = ECG_mean(Pstart-0.05*fs+2:Pstart+0.05*fs+1);
        M(0.3*fs+1:end,3)       = ECG_mean(Pstart+2+0.05*fs:Pstart+1+Tstop);
        a = (M'*M)\M'*ECG(peaks(i)-Pstart:peaks(i)+Tstop)';
        ECG_temp(peaks(i)-Pstart:peaks(i)+Tstop) = a(1)*M(:,1)'+a(2)*M(:,2)'+a(3)*M(:,3)';
    end
end

% compute residual
residual = ECG - ECG_temp;
end


function [SelectedResidual,ChannelNb] = ChannelSelectionOrCombination(FECG)
% This function is used to select one of the four abdominal channels
% that are available or to combine information from these channels
% (e.g. using PCA) before FQRS detection
ChannelNb = 1;
SelectedResidual = FECG(:,ChannelNb); % channel 1 is arbitrarily selected here
end

function FECG = ResidualPostProcessing(FECG)
% if postprocessing is performed on the residuals.
end




function peaks = PeakDetection(x,ff,varargin)
%
% peaks = PeakDetection(x,f,flag),
% R-peak detector based on max search
%
% inputs:
% x: vector of input data
% f: approximate ECG beat-rate in Hertz, normalized by the sampling frequency
% flag: search for positive (flag=1) or negative (flag=0) peaks. By default
% the maximum absolute value of the signal, determines the peak sign.
%
% output:
% peaks: vector of R-peak impulse train
%
% Notes:
% - The R-peaks are found from a peak search in windows of length N; where
% N corresponds to the R-peak period calculated from the given f. R-peaks
% with periods smaller than N/2 or greater than N are not detected.
% - The signal baseline wander is recommended to be removed before the
% R-peak detection
%
%
% Open Source ECG Toolbox, version 1.0, November 2006
% Released under the GNU General Public License
% Copyright (C) 2006  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- GIPSA-Lab, INPG, Grenoble, France
% reza.sameni@gmail.com

% Last modified 03_02_2013: Joachim Behar, IPMG Oxford.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

N = length(x);
peaks = zeros(1,N);

th = .5;
rng = floor(th/ff);

if(nargin==3),
    flag = varargin{1};
else
    flag = abs(max(x))>abs(min(x));
end

if(flag)
    for j = 1:N,
        %         index = max(j-rng,1):min(j+rng,N);
        if(j>rng && j<N-rng)
            index = j-rng:j+rng;
        elseif(j>rng)
            index = N-2*rng:N;
        else
            index = 1:2*rng;
        end
        
        if(max(x(index))==x(j))
            peaks(j) = 1;
        end
    end
else
    for j = 1:N,
        %         index = max(j-rng,1):min(j+rng,N);
        if(j>rng && j<N-rng)
            index = j-rng:j+rng;
        elseif(j>rng)
            index = N-2*rng:N;
        else
            index = 1:2*rng;
        end
        
        if(min(x(index))==x(j))
            peaks(j) = 1;
        end
    end
end


% remove fake peaks
I = find(peaks);
d = diff(I);
% z = find(d<rng);
peaks(I(d<rng))=0;


peaks = find(peaks);
end



