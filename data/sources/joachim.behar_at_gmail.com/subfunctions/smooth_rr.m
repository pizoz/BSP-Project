function FQRSc = smooth_rr(FQRS,rMQRS,ecg,fs,debug)
% look for smoothing RR time series to remove additional FQRS and fill 
% missing FQRS. This showed to improve the PCinC scoring results in term 
% of F1 measure by >1% on the training set-a. 
% This function takes into account a certain number of cases
% where there might be missing or extra QRS and given context and prior 
% on physiology decides where to add/remove QRS.
%
% inputs
%   FQRS:   FQRS time series (required)
%   rMQRS:  reference MQRS (optional). If rMQRS is inputted then the
%           algorithm looks for FQRS that are overlapping with the MQRS. If
%           that is the case then the FQRS fiducial is moved to the centre
%           between the previous and following FQRS. This is thought to
%           be more accurate.
%   ecg:    residual signal for context (optional). Instead of guessing the
%           FQRS location out of complete context, we look at the sign of
%           the known FQRS and look for a local max/min around the expected
%           position where we should have a FQRS.
%   fs:     sampling frequency (optional, default: 1kHz)
%
% output
%   FQRS:   outputed smoothed FQRS time series. 
%   (empty nargout: in the case no output is specified the results
%   are plotted)
%
% IMPORTANT NOTE: this code is designed for foetal ecg which means that the
% constants used to assess whether a FQRS position is missing or extra are
% based on known typical foetal heart rate. If you want to use this code for
% adult ecg then you should change the constants that are listed in the
% code under 'constants' to reflect adult ECG physiology.
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
if nargin<1; error('smooth_rr: wrong number of input arguments \n'); end;
if nargin<2; rMQRS=[]; end;
if nargin<3; ecg=[]; end;
if nargin<4; fs=1000; end;
if nargin<5; debug=0; end;

% == constants
WINDOW  = 0.020*fs; % window in which looking for local max/min in case of missing QRS
MIN_FRR = 0.35*fs; % minimal foetal RR interval corresponds to 171.4bpm
MAX_FRR = 0.5*fs; % maximal foetal RR interval corresponds to 120bpm
EXTRA_BEAT_THRES  = 0.7; % at what threshold do we consider the QRS to be a FP
MISSED_BEAT_THRES = 1.75; % at what threshold do we consider the QRS to be a FN
if ~isempty(ecg); SIGN = sign(median(ecg(FQRS))); end; % sign of the peaks
MAX_NB_ITER = 5000; % this is to avoid infine loops
MQRS_FQRS_TOL = 0.050; % tolerance for assuming MQRS and MQRS overlaps (in sec)

% == core function
try
    FQRSc = FQRS;
    qrsnb = 2;
    compt = 0;
    while qrsnb<length(FQRSc)-1 && compt<MAX_NB_ITER
        % med RR interval over the past 5 beats
        if qrsnb>6
            med = median(diff(FQRSc(qrsnb-5:qrsnb)));
        else
            med = median(diff(FQRSc(qrsnb:qrsnb+5)));
        end
        % check that the med computed makes sense (i.e. should be between
        % 350 and 500ms i.e. 120 and 172bpm) otherwise do nothing 
        % because prediction might be rubbish
        if med>MIN_FRR && med<MAX_FRR
            dTplus = FQRSc(qrsnb+1)-FQRSc(qrsnb);  % RR interval forward
            dTminus = FQRSc(qrsnb)-FQRSc(qrsnb-1);  % RR interval backward

            if dTplus<EXTRA_BEAT_THRES*med && dTminus<1.2*med
                % == case 1: extra beat
                FQRSc(qrsnb+1) = [];                     % remove extra beat
            elseif dTplus>MISSED_BEAT_THRES*med && dTminus>0.7*med
                % == case 2: missed beat
                if ~isempty(ecg)
                    % if we have context then look for a local max/min
                    if SIGN>0
                        [~,mind] = max(ecg(FQRSc(qrsnb)+med-WINDOW:FQRSc(qrsnb)+med+WINDOW));
                    else
                        [~,mind] = min(ecg(FQRSc(qrsnb)+med-WINDOW:FQRSc(qrsnb)+med+WINDOW));
                    end
                    MissedFQRS = (mind-WINDOW-1)+FQRSc(qrsnb)+med;
                else
                    % otherwise guess location
                    if dTplus<3*med && dTplus>1.5*med
                        % if the dTplus makes sense then better to guess in
                        % the middle of current and next peak
                        MissedFQRS = round(FQRSc(qrsnb)+(dTplus)/2);
                    else
                        % otherwise use the med to predict (less precise)
                        MissedFQRS = round(FQRSc(qrsnb)+med);
                    end
                end
                % insert missed FQRS where it belongs
                FQRSc = [FQRSc(1:qrsnb) MissedFQRS FQRSc(qrsnb+1:end)];
                qrsnb = qrsnb+1;
            else
                % == case 3: normal detection
                qrsnb = qrsnb+1;
            end
        else
            qrsnb = qrsnb+1;
        end
        compt = compt+1;
    end

    % == managing MQRS/FQRS overlap
    if ~isempty(rMQRS)
        FQRStemp = FQRSc;
        [~,d] = dsearchn(rMQRS',FQRStemp');
        ind2cor = find(d<MQRS_FQRS_TOL*fs); % indices to correct
        for cc=1:length(ind2cor)
            % replace all overlapping FQRS/MQRS by the midpoint between
            % previous and next FQRS position which is assumed to be a better
            % estimate
            if ind2cor~=1 & ind2cor~=length(d)
                FQRStemp(ind2cor(cc)) = round((FQRStemp(ind2cor(cc)-1)+FQRStemp(ind2cor(cc)+1))/2);
            end
        end
        FQRSc = FQRStemp; 
    end

catch ME
    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
    FQRSc = FQRS;
end

% == plots
if debug || isempty(nargout)
   if isempty(ecg)
        ax(1)=subplot(211); plot(FQRS/fs,'o','LineWidth',2); hold on, plot(FQRSc/fs,'+r','LineWidth',2); 
        legend('initial','corrected'); title('FQRS correction'); xlabel('Time [sec]');
   else
        tm = 1/fs:1/fs:length(ecg)/fs; 
        ax(1)=subplot(211); plot(tm,ecg);
        hold on; plot(tm(FQRS),ecg(FQRS),'o','LineWidth',2); hold on, plot(tm(FQRSc),ecg(FQRSc),'+r','LineWidth',2);
        legend('ecg residual','initial','corrected'); title('FQRS correction'); xlabel('Time [sec]');     
   end
   hr = 60./(diff(FQRS)/fs);
   hrc = 60./(diff(FQRSc)/fs);
   ax(2)=subplot(212); plot(FQRS(1:end-1)/fs,hr,'--o'); hold on, plot(FQRSc(1:end-1)/fs,hrc,'--r+');
   legend('HR','HRc'); ylabel('HR [bpm]'); xlabel('Time [sec]');
   linkaxes(ax,'x');
   set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold');
end

end




