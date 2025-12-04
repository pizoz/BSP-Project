function [fetal_QRSAnn_est,QT_Interval] = physionet2013(tm,ECG)
% This algorithm for Physionet/CinC competition 2013. Entry 1 No PCA
% [fetal_QRSAnn_est,QT_Interval] = physionet2013(tm,ECG) where the inputs and outputs are specified
% below.
% inputs:
%   ECG: 4x60000 (4 channels and 1min of signal at 1000Hz) matrix of
%   abdominal ECG channels.
%   tm : Nx1 vector of time in milliseconds
% output:
%   fetal_QRSAnn_est: FQRS Annotations each value indicates the position of
%   one  of the FQRS detected by the algorithm.
%   QT_Interval:   1x1 estimated fetal QT duration 
%
% Modified By: Aruna Deogire, MET's Institute of Engineering,Nashik,India
% email - (arunadeogire@gmail.com)
% Last updated: August 24th, 2013 by Aruna Deogire


% ---- check size of ECG ----
if size(ECG,2)>size(ECG,1)
    ECG = ECG';
end

% Remove NaN
A=ECG(:,1);
A(isnan(A))=0;
ECG(:,1)=A;

B=ECG(:,2);
B(isnan(B))=0;
ECG(:,2)=B;

C=ECG(:,3);
C(isnan(C))=0;
ECG(:,3)=C;

D=ECG(:,4);
D(isnan(D))=0;
ECG(:,4)=D;

% Baseline Wander Removal
% *************************************************************************

Smooth_Period=149; wav_typ='db6'; level=7;
fs      = 1000;             % sampling frequency
N       = size(ECG,2);      % number of abdominal channels
debug   = 0;                % enter debug mode?
% _________________________________________________________________________

% Preprocessing 
% *************************************************************************
    [FilteredECG] = preprocessing(ECG,fs);
% _________________________________________________________________________


% Mother QRS Detection 
% *************************************************************************
e1=FilteredECG(:,1); 
MQRS1=PeakDetection(FilteredECG(:,1),1/fs);Mpoint1=e1(MQRS1);

e2=FilteredECG(:,2); 
MQRS2=PeakDetection(FilteredECG(:,2),1/fs);Mpoint2=e2(MQRS2);

e3=FilteredECG(:,3); 
MQRS3=PeakDetection(FilteredECG(:,3),1/fs);Mpoint3=e3(MQRS3);

e4=FilteredECG(:,4); 
MQRS4=PeakDetection(FilteredECG(:,4),1/fs);Mpoint4=e4(MQRS4);
% _________________________________________________________________________

% Mother QRS Cancellation
% *************************************************************************

FECG1 = MECGcancellation(MQRS1,FilteredECG(:,1)',fs,20);
FECG2 = MECGcancellation(MQRS2,FilteredECG(:,2)',fs,20);
FECG3 = MECGcancellation(MQRS3,FilteredECG(:,3)',fs,20);
FECG4 = MECGcancellation(MQRS4,FilteredECG(:,4)',fs,20);
% _________________________________________________________________________


% Post Processing of Residual Signal
% *************************************************************************
[c,l]=wavedec(FECG1,level,wav_typ); 
FECG12 = wden(c,l,'minimaxi','s','sln',level,wav_typ);  

[c,l]=wavedec(FECG2,level,wav_typ); 
FECG22 = wden(c,l,'minimaxi','s','sln',level,wav_typ);

[c,l]=wavedec(FECG3,level,wav_typ); 
FECG32 = wden(c,l,'minimaxi','s','sln',level,wav_typ);

[c,l]=wavedec(FECG4,level,wav_typ); 
FECG42 = wden(c,l,'minimaxi','s','sln',level,wav_typ);

%Array of post processed Fetal ECG
Fetal_ECG = [FECG12 FECG22 FECG32 FECG42]; 
% _________________________________________________________________________

% Fetal QRS Detection by PeakDetection Algorithm
%*************************************************************************

FQRS1=PeakDetection(FECG12,2/fs);Fpoint1=FECG12(FQRS1);

FQRS2=PeakDetection(FECG22,2/fs);Fpoint2=FECG22(FQRS2);

FQRS3=PeakDetection(FECG32,2/fs);Fpoint3=FECG32(FQRS3);

FQRS4=PeakDetection(FECG42,2/fs);Fpoint4=FECG42(FQRS4);
% _________________________________________________________________________

% Fetal QRS Selection Outof four Leads
%*************************************************************************
FT1=FQRS1'; FT2=FQRS2'; FT3=FQRS3'; FT4=FQRS4';

LA = [length(FT1) length(FT2) length(FT3) length(FT4)];

T = max(LA);

if ( T == LA(:,1))
    FQRS=FT1';
    S = FECG12;
else
    if ( T == LA(:,2))
    FQRS=FT2';
    S = FECG22;
    else
       if ( T == LA(:,3))
          FQRS=FT3';
          S = FECG32;
       else
           FQRS=FT4';
           S = FECG42;
       end
    end
end
% _________________________________________________________________________

% QT Point Detections 
% Q Point Detection
%*************************************************************************
LQ = length(FQRS);
 
stepq=fs*0.02; % interval to search Q 
for i=2:LQ-1
j=FQRS(i);
[Q,tq]=min(S(j-stepq:j));
j=j-stepq+tq-1;
Qindex(i-1)=j;
Qpoint(i-1)=S(j);
end


% Plot ECG signal with QRS onset marked
% figure;
% plot(S), hold on, plot(FQRS,Fpoint,'r.',Qindex,Qpoint,'k.'),xlim([10000 15000]);grid on; hold off;
% _________________________________________________________________________

% S Point Detection
%*************************************************************************
LS = length(FQRS);

steps=fs*0.02; % interval to search Q 
for i=2:LS-1
j=FQRS(i);
[Q,ts]=min(S(j:j+steps));
j=j+ts;
Sindex(i)=j;
Spoint(i)=S(j);
end


% Plot ECG signal with QRS onset marked
% figure;
% plot(S), hold on, plot(FQRS,Fpoint,'r.',Qindex,Qpoint,'k.',Sindex,Spoint,'k.'),xlim([20000 30000]);grid on; hold off;
% _________________________________________________________________________

% T Point Detection
%*************************************************************************
LT = length(Sindex);

stept=fs*0.02; % interval to search T 
for i=2:LT-1
j=Sindex(i);
[Q,tt]=max(S(j:j+stept));
j=j+tt;
Tindex(i)=j;
Tpoint(i)=S(j);
end


% T End Detection
%*************************************************************************
LTe = length(Sindex);
% for i=2:LTe-1
% j=Tindex(i);
% while 1
% j=j+1;
% a=S(j+1)-S(j);
% b=S(j)-S(j-1);
% if a>=0 & b<0
% Tendpoint(i)=S(j);
% Tendindex(i)=j;
% break;
% end
% end
% end

stepte=fs*0.01; % interval to search T 
for i=2:LTe-1
j=Tindex(i);
[Q,tte]=min(S(j:j+stepte));
j=j+tte;
Tendindex(i)=j;
Tendpoint(i)=S(j);
end

QT = Qindex - Tendindex;
QT_Interval = median(QT);


% If Debug is Activated
%*************************************************************************
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
% QT_Interval         = QTavg;
end
% _________________________________________________________________________

%*************************************************************************
                                % Functions Started
%*************************************************************************

% Preprocessing Function
% *************************************************************************
function [FilteredECG] = preprocessing(ECG,fs)


    fs=1000;Smooth_Period=149; wav_typ='db6'; level=6;

    [c,l]=wavedec(ECG(:,1),level,wav_typ); AD1 = wden(c,l,'minimaxi','s','sln',level,wav_typ);   
    Pattern=smooth(AD1,Smooth_Period);AF1= AD1-Pattern;
    
    [c,l]=wavedec(ECG(:,2),level,wav_typ); AD2 = wden(c,l,'minimaxi','s','sln',level,wav_typ);   
    Pattern=smooth(AD2,Smooth_Period);AF2= AD2-Pattern;
    
    [c,l]=wavedec(ECG(:,3),level,wav_typ); AD3 = wden(c,l,'minimaxi','s','sln',level,wav_typ);   
    Pattern=smooth(AD3,Smooth_Period);AF3= AD3-Pattern;
   
    [c,l]=wavedec(ECG(:,4),level,wav_typ); AD4 = wden(c,l,'minimaxi','s','sln',level,wav_typ);   
    Pattern=smooth(AD4,Smooth_Period);AF4= AD4-Pattern;

    FilteredECG=[AF1 AF2 AF3 AF4]; %Vector of Filtered Signals
end
% _________________________________________________________________________

% MECGcancellation Function
% *************************************************************************
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
% _________________________________________________________________________

%PeakDetection Function
% *************************************************************************

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

% _________________________________________________________________________

