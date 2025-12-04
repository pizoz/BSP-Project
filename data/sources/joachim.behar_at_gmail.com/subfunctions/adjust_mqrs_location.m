function cMQRS = adjust_mqrs_location(ecg,MQRS,fs,debug)
% the point of this function is to adjust the MQRS location since
% MQRS are not at the exact same location from a lead to another 
% due to conduction delay and/or quantification limitation. 
% This can cause the template substraction technique to be
% imprecise (a small lack of alignement can cause poor results) if not
% performed. This function showed to improve the results of the PCinC2013.
%
% inputs
%   ecg:    the ecg signals
%   MQRS:   MQRS location (in ms)
%   fs:     sampling frequency
%
% outputs
%   cMQRS: corrected MQRS location
%   (if nargout is empty then plot the initial and corrected time series)
%
%
% FECG extraction toolbox, version 1.0, Sept 2013
% Released under the GNU General Public License
%
% Copyright (C) 2013  Joachim Behar
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2013
% joachim.behar@eng.ox.ac.uk
%
% Last updated : 28-08-2013
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
if nargin<2; error('adjust_mqrs_location: wrong number of input arguments \n'); end;
if nargin<3; fs=1000; end;
if nargin<4; debug=0; end;
if size(ecg,2)>size(ecg,1); ecg=ecg'; end;

% == general
cMQRS = zeros(length(MQRS),1); % corrected MQRS vector
NB_MQRS = length(MQRS);
WINDOW  = ceil(fs*0.030); % allow a maximum of 30ms shift
NB_SAMPLES = length(ecg);
signVec = zeros(30,1);

% == core function
try
    % == estimates if we are looking for local max or min on the first few
    % peaks
    for i=1:30
        startLooking = MQRS(i)-WINDOW;
        stopLooking = MQRS(i)+WINDOW;
        if startLooking<=0; startLooking=1; end;
        if stopLooking>=NB_SAMPLES; stopLooking=NB_SAMPLES; end;
        [~,indm] = max(abs(ecg(startLooking:stopLooking)));
        if startLooking==1; signVec(i) = MQRS(i)+indm-1;
        else signVec(i) = MQRS(i)+indm-WINDOW; end;
    end
    SIGN = sign(median(ecg(signVec)));

    % == local refinement to correct for sampling freq resolution
    for i=1:NB_MQRS
        if MQRS(i)>WINDOW && MQRS(i)+WINDOW<NB_SAMPLES
            if SIGN>0
                [~,indm] = max(ecg(MQRS(i)-WINDOW:MQRS(i)+WINDOW));
                cMQRS(i) = MQRS(i)+indm-WINDOW-1;
            else
                [~,indm] = min(ecg(MQRS(i)-WINDOW:MQRS(i)+WINDOW));
                cMQRS(i) = MQRS(i)+indm-WINDOW-1;
            end
        elseif MQRS(i)<WINDOW
            % managing left boder
             if SIGN>0
                [~,indm] = max(ecg(1:MQRS(i)+WINDOW));
                cMQRS(i) = MQRS(i)+indm-1;
            else
                [~,indm] = min(ecg(1:MQRS(i)+WINDOW));
                cMQRS(i) = MQRS(i)+indm-1;
            end           
        elseif MQRS(i)+WINDOW>NB_SAMPLES
            % managing right border
             if SIGN>0
                [~,indm] = max(ecg(MQRS(i)-WINDOW:end));
                cMQRS(i) = MQRS(i)+indm-WINDOW-1;
            else
                [~,indm] = min(ecg(MQRS(i)-WINDOW:end));
                cMQRS(i) = MQRS(i)+indm-WINDOW-1;
            end           
        else
            cMQRS(i) = MQRS(i);
        end
    end
    
catch ME
    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
    cMQRS = MQRS;
end

% == plots
if debug || isempty(nargout)
    plot(ecg); hold on, plot(MQRS, ecg(MQRS),'ok');
    hold on, plot(cMQRS, ecg(cMQRS),'or');
    legend('ecg','MQRS initial','MQRS adjusted');
    set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold');
end


end





