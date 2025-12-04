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
% Last updated: 2. 9. 2013, Jakub Kuzilek (jakub.kuzilek@gmail.com)
% load test.mat
addpath(fullfile(pwd,'QRS_detection'));

% ---- check size of ECG ----
if size(ECG,2)>size(ECG,1)
    ECG = ECG';
end

fs      = 1000;             % sampling frequency
N       = size(ECG,2);      % number of abdominal channels
debug   = 0;                % enter debug mode?

% ---- preprocessing ----
[FilteredECG] = preprocessing(ECG,fs);

% MQRS = PeakDetection(FilteredECG(:,1),1/fs);
[~,indx] = sort(kurtosis(FilteredECG));
MQRS = QRSdetection(FilteredECG(:,sort(indx(1:2))),fs,0);

% ---- MECG cancellation ----
for i=1:N               % run algorithm for each channel
    FECG(:,i) = MECGcancellation(MQRS,FilteredECG(:,i)',fs,20,0); % 5
end

% ---- FQRS detection ----
FQRS = {};
stdfqrs = [];
fqrsnum = [];

for m = 1:N
    SelectedResidual = FECG(:,m);
    FQRS{end+1} = PeakDetection(SelectedResidual,2/fs);
    stdfqrs(end+1) = std(diff(FQRS{end}));
    fqrsnum(end+1) = length(FQRS{end});
    FQRS{end+1} = correctFQRS(SelectedResidual,FQRS{end});
    stdfqrs(end+1) = std(diff(FQRS{end}));
    fqrsnum(end+1) = length(FQRS{end});
    FQRS{end+1} = findFQRS(SelectedResidual);
    stdfqrs(end+1) = std(diff(FQRS{end}));
    fqrsnum(end+1) = length(FQRS{end});
    FQRS{end+1} = correctFQRS(SelectedResidual,FQRS{end});
    stdfqrs(end+1) = std(diff(FQRS{end}));
    fqrsnum(end+1) = length(FQRS{end});
end


FQRS(fqrsnum < 20) = [];
stdfqrs(fqrsnum < 20) = [];
fqrsnum(fqrsnum < 20) = [];

[v, lok] = min(stdfqrs);
FQRS = FQRS{lok};

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

rmpath(fullfile(pwd,'QRS_detection'));
end

function [FilteredECG] = preprocessing(ECG,fs)
% ---- preprocess the data ----
time = 1:length(ECG);
for m = 1:min(size(ECG))
    if any(isnan(ECG(:,m)))
        ECG(:,m) = pchip(time(~isnan(ECG(:,2))),ECG(~isnan(ECG(:,2)),2),time);
    end
end
FilteredECG = ECG;
for m = 1:size(ECG,2)
    FilteredECG(:,m) = filterIsoline(ECG(:,m),fs);
end
aa = filter(ones(1,50)/50,1,[FilteredECG; zeros(25,size(ECG,2))]);
FilteredECG = FilteredECG - aa(26:end,:);

end

function residual = MECGcancellation(peaks,ECG,fs,nbCycles, rl)
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
% Last updated: 2. 9. 2013, Jakub Kuzilek (jakub.kuzilek@gmail.com)
%
% [1] Martens S et al. A robust fetal ECG detection method for
% abdominal recordings. Physiol. Meas. (2007) 28(4) 373�388

% ---- constants ----
r                   = nbCycles;
Pstart              = 0.25*fs-1;
Tstop               = 0.45*fs;
N                   = length(peaks);                    % number of MECG QRS
ECG_temp            = zeros(1,length(ECG));

% ---- ECG template ----
if mod(N,r)
    iters = floor(N/r)+1;
else
    iters = floor(N/r);
end
for m = 0:iters
    ECG_last_r_cycles   = zeros(0.7*fs,r);
    ECG_mean = zeros(0.7*fs,1);
    cnt = 1;
    for i=m*r+1:min([m*r+r N])
        if peaks(i)>Pstart && length(ECG)-peaks(i)>Tstop && (i<N && peaks(i+1)-peaks(i) > Tstop + 50)
            ECG_last_r_cycles(:,cnt) = ECG(peaks(i)-Pstart:peaks(i)+Tstop)';
        elseif  peaks(i)<=Pstart
            ECG_last_r_cycles(:,cnt) = [zeros(1,Pstart-peaks(i)+1), ECG(1:peaks(i)+Tstop)]';
        elseif i == N && length(ECG)-peaks(i)<Tstop
            ECG_last_r_cycles(:,cnt) = [ECG(peaks(i)-Pstart:end), zeros(1,Tstop+peaks(i)-length(ECG))]';
        elseif (i<N && peaks(i+1)-peaks(i) < Tstop + 50)
            lengthReduction = Tstop - peaks(i+1) + peaks(i)+50;
            if lengthReduction < 50
                lengthReduction = 50;
            end
            ECG_last_r_cycles(:,cnt) = [ ECG(peaks(i)-Pstart:peaks(i)+Tstop-lengthReduction) zeros(1, lengthReduction)]';
        end
        cnt = cnt+1;
    end
    
    ECG_mean =mean(ECG_last_r_cycles,2);
    
    % ---- MECG cancellation ----
    for i= m*r+1:min([m*r+r N])%1:N
        if peaks(i)>Pstart && length(ECG)-peaks(i)>Tstop && (i<N && peaks(i+1)-peaks(i) > Tstop + 50)
            M  = zeros (0.7*fs,3);
            M(1:0.2*fs,1)           = ECG_mean(1:Pstart-0.05*fs+1);
            M(0.2*fs+1:0.3*fs,2)    = ECG_mean(Pstart-0.05*fs+2:Pstart+0.05*fs+1);
            M(0.3*fs+1:end,3)       = ECG_mean(Pstart+2+0.05*fs:Pstart+1+Tstop);
            out = bestFit(peaks(i), ECG, M(:,1)'+M(:,2)'+M(:,3)', Pstart, Tstop);
            if out-Pstart<1
                out = peaks(i);
            end
            a = (M'*M)\M'*ECG(out-Pstart:out+Tstop)';
            ECG_temp(out-Pstart:out+Tstop) = a(1)*M(:,1)'+a(2)*M(:,2)'+a(3)*M(:,3)';
        elseif peaks(i)<=Pstart
            M  = zeros (0.7*fs,3);
            M(1:0.2*fs,1)           = ECG_mean(1:Pstart-0.05*fs+1);
            M(0.2*fs+1:0.3*fs,2)    = ECG_mean(Pstart-0.05*fs+2:Pstart+0.05*fs+1);
            M(0.3*fs+1:end,3)       = ECG_mean(Pstart+2+0.05*fs:Pstart+1+Tstop);
            out = bestFit(peaks(i), ECG, M(:,1)'+M(:,2)'+M(:,3)', Pstart, Tstop);
            if Pstart < out
                template = ECG(out-Pstart:out+Tstop);
                a = (M'*M)\M'*template';
                ECG_temp(out-Pstart:out+Tstop) = a(1)*M(:,1)'+a(2)*M(:,2)'+a(3)*M(:,3)';
            else
                template = [zeros(1,Pstart-out+1), ECG(1:out+Tstop)];
                a = (M'*M)\M'*template';
                ECG_temp(1:out+Tstop) = a(1)*M(Pstart-out+2:end,1)'+a(2)*M(Pstart-out+2:end,2)'+a(3)*M(Pstart-out+2:end,3)';
            end
        elseif (i<N && peaks(i)-peaks(i+1) < Tstop + 50)
            lengthReduction = Tstop - peaks(i+1) + peaks(i)+50;
            if lengthReduction < 50
                lengthReduction = 50;
            end
            M  = zeros (0.7*fs-lengthReduction,3);
            M(1:0.2*fs,1)           = ECG_mean(1:Pstart-0.05*fs+1);
            M(0.2*fs+1:0.3*fs,2)    = ECG_mean(Pstart-0.05*fs+2:Pstart+0.05*fs+1);
            M(0.3*fs+1:end,3)       = ECG_mean(Pstart+2+0.05*fs:Pstart+1+Tstop-lengthReduction);
            out = bestFit(peaks(i), ECG, M(:,1)'+M(:,2)'+M(:,3)', Pstart, Tstop-lengthReduction);
            if Pstart < out
                template = ECG(out-Pstart:out+Tstop-lengthReduction);
                a = (M'*M)\M'*template';
                ECG_temp(out-Pstart:out+Tstop-lengthReduction) = a(1)*M(:,1)'+a(2)*M(:,2)'+a(3)*M(:,3)';
            else
                template = [zeros(1,Pstart-out+1), ECG(1:out+Tstop-lengthReduction)];
                a = (M'*M)\M'*template';
                ECG_temp(1:out+Tstop-lengthReduction) = a(1)*M(Pstart-out+2:end,1)'+a(2)*M(Pstart-out+2:end,2)'+a(3)*M(Pstart-out+2:end,3)';
            end
        end
    end
end

% compute residual
residual = ECG - ECG_temp;
x = 0;
for m = 1:length(peaks)
    x = x + sum(residual(max([peaks(m)-Pstart 1]):min([peaks(m)+Tstop length(residual)])).^2)/700;
end
x = x/length(peaks);
if x > 15 && rl < 100
    rl = rl+1;
    residual = MECGcancellation(peaks,residual,fs,20,rl);
end
end

function out = bestFit(peak, ECG, template, Pstart, Tstop)

err = zeros(1,101);
for m = peak-50:peak+50
    if m-Pstart<1
        err(m-peak+51) = ([zeros(1, Pstart-m+1), ECG(1:m+Tstop)]-template)*([zeros(1, Pstart-m+1), ECG(1:m+Tstop)]-template)';
    else
        err(m-peak+51) = (ECG(m-Pstart:m+Tstop)-template)*(ECG(m-Pstart:m+Tstop)-template)';
    end
end

[v, pos] = min(err);

out = peak-51+pos;

end


function out = computeECGMEAN(ECG, fs, r, Tstop, Pstart,peaks)
ECG_last_r_cycles   = zeros(0.7*fs,r);
for i=1:r
    peak_nb = peaks(i+1);
    ECG_last_r_cycles(:,i) = ECG(peak_nb-Pstart:peak_nb+Tstop)';
end
out = mean(ECG_last_r_cycles,2);
end

function [SelectedResidual,ChannelNb] = ChannelSelectionOrCombination(FECG)
% This function is used to select one of the four abdominal channels
% that are available or to combine information from these channels
% (e.g. using PCA) before FQRS detection
ChannelNb = 1;
SelectedResidual = FECG(:,ChannelNb); % channel 1 is arbitrarily selected here
end

function FECG = ResidualPostProcessing(FECG)
for m = 1:size(FECG,2)
    for n = 30:35
        FECG(:,m) = adaptiveANC(FECG(:,m),1000,n);
    end
end
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



