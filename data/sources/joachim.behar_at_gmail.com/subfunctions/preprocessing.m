function [ecgf] = preprocessing(ecg,paramStruct,varargin)
% this function preprocess the ecg prior applying any source separation 
% technique.
%
% inputs  
%   ecg: 4x60000 (4 channels and 1min of signal at 1000Hz) matrix of
%   abdominal ecg channels.
%   paramStruct: the following parameters needs to be specified
%       paramStruct.fs
%       paramStruct.NotchCheck
%       paramStruct.fhigh
%       paramStruct.fbas
%       paramStruct.fhighNbCoeff
%       paramStruct.fbasNbCoeff
%   varargin
%       createFilters: recreate the filter coefficients?
%       path2filtersCoeff: path to where to save the filter 
%                          coefficients or load them from.
%       debug: enter debug mode
%
% outputs
%   ecgf: normalized and filtered ecg. Filtering
%                     includes Notch, high frequency, baseline wander removal
%                     and notmalisation in range [-1 1]. A non linear step is
%                     applied with the tanh so this might note be suited for
%                     morphological analysis but is usefull for FQRS detection.
%
% IMPORTANT NOTE: too many coefficients  for the IIR filters resulted in 
% ringing. So need to use as few as possible and avoid any useless filtering 
% step which is why I am looking whether there is a peak at 50/60Hz i.e. 
% whether or not to perform a Notch filtering step.
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

createFilters = 0;
path2filtersCoeff = 'subfunctions/FiltersCoeff/';
debug = 0;
% == manage inputs
switch nargin
    case 2
    case 3
        createFilters=varargin{1};       
    case 4
        createFilters=varargin{1};
        path2filtersCoeff=varargin{2};
    case 5
        createFilters=varargin{1};
        path2filtersCoeff=varargin{2};
        debug=varargin{3};
    otherwise
        error('preprocessing: wrong number of input arguments \n');              
end
if size(ecg,2)>size(ecg,1); ecg = ecg'; end;
if ~exist(path2filtersCoeff,'dir'); mkdir(path2filtersCoeff); end;
currFold = pwd;

% == core function
try 
    % some AECG (the artificial ones) saturate sometimes. Because they were
    % rescaled (after being saved in WFDB format) then this leads to weird
    % disptoportionated saturation. The lines below are aimed at solving
    % this problem. (If I open the records in physical units and run the
    % filters, I am gettnig some NaN..)
    % - manage saturation for wfdb records
    % NOTE: this is wrong but works by luck somehow- will need to be fixed
    % at some point!
    [ind_t,ind_ch] = find(ecg>32750);
    ecg(ind_t,ind_ch) = median(ecg(ind_t-5:ind_t-5,ind_ch));
    [ind_t,ind_ch] = find(ecg<-32750);
    ecg(ind_t,ind_ch) = median(ecg(ind_t-5:ind_t-5,ind_ch)); % should be -32767 but they have rescaled the artificial signals

    % - manage NaN for .cvs files
    ind_nan = isnan(ecg);
    for ind_ch_nan=1:4
        if sum(ind_nan(:,ind_ch_nan))>0
            nan_ch_pos = find(ind_nan(:,ind_ch_nan)); % position of the NaN
            if ~isempty(nan_ch_pos)
                for jj=1:length(nan_ch_pos)
                    
                    start = nan_ch_pos(jj)-5;
                    stop = nan_ch_pos(jj)+5;
                    if start<0; start=0; end;
                    if stop>length(ecg); stop=length(ecg); end;
                    
                    vecTrans = ecg(start:stop,ind_ch_nan);
                    vecTrans(isnan(vecTrans)) = [];
                    if ~isempty(vecTrans)
                        ecg(nan_ch_pos(jj),ind_ch_nan) = median(vecTrans);
                    end
                end 
            end
        end
    end
    
    % == Filtering
    if createFilters
        % == 60Hz Notch
        wo = 60/(paramStruct.fs/2);  bw = wo/35;
        [bn_60Hz,an_60Hz] = iirnotch(wo,bw);
        dlmwrite(strcat(path2filtersCoeff,'bn_60Hz.txt'),bn_60Hz,'precision',20);
        dlmwrite(strcat(path2filtersCoeff,'an_60Hz.txt'),an_60Hz,'precision',20);    
        % == 50Hz Notch
        wo = 50/(paramStruct.fs/2);  bw = wo/35;
        [bn_50Hz,an_50Hz] = iirnotch(wo,bw);
        dlmwrite(strcat(path2filtersCoeff,'bn_50Hz.txt'),bn_50Hz,'precision',20);
        dlmwrite(strcat(path2filtersCoeff,'an_50Hz.txt'),an_50Hz,'precision',20);
        % == high freq
        [b_lp,a_lp] = butter(paramStruct.fhighNbCoeff,paramStruct.fhigh/(paramStruct.fs/2),'high');
        dlmwrite(strcat(path2filtersCoeff,'b_lp.txt'),b_lp,'precision',20);
        dlmwrite(strcat(path2filtersCoeff,'a_lp.txt'),a_lp,'precision',20);
        % == baseline
        [b_bas,a_bas] = butter(paramStruct.fbasNbCoeff,paramStruct.fbas/(paramStruct.fs/2),'high');
        dlmwrite(strcat(path2filtersCoeff,'b_bas.txt'),b_bas,'precision',20);
        dlmwrite(strcat(path2filtersCoeff,'a_bas.txt'),a_bas,'precision',20);
        %fvtool(b,a);     
    else
        cd(path2filtersCoeff);
            an_60Hz = importdata('an_60Hz.txt');
            bn_60Hz = importdata('bn_60Hz.txt');
            an_50Hz = importdata('an_50Hz.txt');
            bn_50Hz = importdata('bn_50Hz.txt');
            b_lp = importdata('b_lp.txt');
            a_lp = importdata('a_lp.txt');
            b_bas = importdata('b_bas.txt');
            a_bas = importdata('a_bas.txt');
        cd(currFold)
    end

    % == Notch fileting
    if paramStruct.NotchCheck
        % == remove 60Hz?
        NeedFiltering = check_notch(ecg(:,1),60);
        if NeedFiltering
            ecgn = filtfilt(bn_60Hz,an_60Hz,ecg);
        else
            ecgn = ecg;
        end
        % == remove 50Hz?
        NeedFiltering = check_notch(ecg(:,1),50);
        if NeedFiltering
            ecgn = filtfilt(bn_50Hz,an_50Hz,ecgn);
        end
    else
        % == Notch fileter at 50Hz AND 60Hz
        ecgn = filtfilt(bn_60Hz,an_60Hz,ecg);
        ecgn = filtfilt(bn_50Hz,an_50Hz,ecgn); 
    end

    % == remove high freq and baseline
    ecgnb  = ecgn-filtfilt(b_lp,a_lp,ecgn);
    ecgnbh = filtfilt(b_bas,a_bas,ecgnb);

    % == normalize the data
    [ecgf,~] = normalise_ecg(ecgnbh,1*paramStruct.fs,5*paramStruct.fs);
    % reason for starting at 1sec is because of filtering initialisation at the
    % border which takes some datapoints

    ecgf(isnan(ecgf))=0; %FIXME
    
catch ME
    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
    ecgf = ecg;
end

% == plots
if debug
    % plot first channel before and after preprocessing
    tm = 1/paramStruct.fs:1/paramStruct.fs:length(ecg)/paramStruct.fs; figure(1);
    ax(1) = subplot(211); plot(tm,ecg(:,1)); legend('before filtering');
    ax(2) = subplot(212); plot(tm,ecgf(:,1),'r','LineWidth',2); 
    legend('after filtering'); linkaxes(ax,'x'); xlabel('Time [sec]'); ylabel('Amplitude (NU)');
    set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold');
end
end


%%Plot PSD       
%     plot_psd(ecg(:,1),ecgnbh(:,1),paramStruct.fs,'welch');
%     global NameCurrent
%     saveas(gcf,NameCurrent(1:3));







