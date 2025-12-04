function ccQRS = smooth_rr_lots(QRS,fs,seglength,debug)
% this function is used to smooth the rr time series based on
% robust smoothing of HR by taking the reciprocal of f:x->60/x. It will return a 
% meaningfull RR time series HOWEVER it can sometimes not correspond to 
% the actual QRS positions depending on initialisation/prediction. 
% (so do not use that for building a template ecg for eg!). Using robust
% smoothing does remove unwanted outliers and create a nice and smooth HR,
% however there might be some trouble with highly variable FHR depending on 'how much 
% the time series is smoothed'. In practical application using this function is fine but the
% raw HR (or HR smoothed with smooth_rr.m only) should be displayed in case
% the function fails- this is, as an example, done in some comercial CTG devices.
% 
% inputs
%   QRS:        observed QRS position in ms (required)
%   fs:         sampling frequency (optional, default: 1KHz)
%   seglength:  length of the segment in sec (optional, default: 60s)
% 
% output
%   ccQRS:       corrected QRS location (in ms)
%   (nargout:    if no output requested the smoothed time series is plotted)
%
%
% FECG extraction toolbox, version 1.0, Sept 2013
% Released under the GNU General Public License
%
% Copyright (C) 2013  Joachim Behar
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2013
% joachim.behar@eng.ox.ac.uk
%
% Last updated : 02-08-2013
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

% == managing inputs
if nargin<1; error('smooth_rr_lots: wrong number of input arguments \n'); end;
if nargin<2; fs = 1000; seglength = 60; debug = 0; end;
if nargin<3; seglength = 60; debug = 0; end;
if nargin<4; debug = 0; end;

% == general
HR = 60./diff(QRS);
START_QRS = 2;
DIST_MAX = 0.020; % max distance tolerated between predicted R-peak and observed R-peak (in ms)
STOPTIME = QRS(end); % stop at the last QRS
MAX_NB_ITER = 1000; % this is to avoid infinite loops
MAX_TIME_SEGMENT = seglength; % length of the segment (in sec)
comp = 1;

% == smooth HR
HRinterp = smooth(QRS(1:end-1),HR,0.05,'rloess'); 
X = QRS(1:end-1);

% == now try to figure out where the corresponding QRS lie
cQRS = []; % corrected QRS table
QRSb = QRS(START_QRS); % QRS buffer
QRStmp = 0;
while QRStmp<STOPTIME && comp<MAX_NB_ITER
    % instantaneous HR at current location?
    instHR = interp1(X,HRinterp,QRSb,'linear');
    % predict where the next QRS should be
    QRStmp = QRSb+60/instHR;
    % is there a detected QRS at this location?
    [k,dis] = dsearchn(QRS',QRStmp);
    
    % decide whether to trust the observation or rely on prediction
    if dis>DIST_MAX
       %then we are missing a point
       cQRS = [cQRS QRStmp];
       QRSb = QRStmp;
    else
       % else keep the point where it is
       cQRS = [cQRS QRS(k)];
       QRSb = QRS(k);
    end
    comp = comp+1;
end

% make sure no duplicated and no prediction is made after MAX_TIME_SEGMENT
ccQRS = unique([QRS(1:START_QRS) cQRS(cQRS<MAX_TIME_SEGMENT)]);

% == debug plots
if debug || isempty(nargout)
   HRnew  = 60./diff(ccQRS);
   figure(1); hold all;
   plot(QRS(1:end-1),HR,'LineWidth',2); 
   plot(QRS(1:end-1),HRinterp+1,'--k','LineWidth',2);
   plot(ccQRS(1:end-1),HRnew+2,'r','LineWidth',2);
   legend('HR old','HR smooth','HR reconstructed from RR output');
   xlabel('Time [sec]'); ylabel('HR [bpm]');
end

% convert to ms
ccQRS = round(ccQRS*fs);

end