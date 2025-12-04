function ecgout = bss_step(ecg,icaecgf,rMQRS,paramStruct,debug)
% Blind Source Separation (BSS) step. In the case of the DEFLATION method
% BSS is applied to the ABD signals and MECG cancellation takes place in
% the source domain. In the other cases then BSS is applied on the ecg
% matrix containing the residuals of the MECG cancellation perfomed in the
% time domain.
%
% inputs
%   ecg:         ecg matrix onto which to perform BSS
%   icaecgf:     ica applied on filtered ecg stack
%   rMQRS:       reference MQRS location in samples
%   paramStruct: important algo parameters
%
% output
%   ecgout:      output ecg that came through the BSS step
%                case DEFLATION: [ICA-TS], [ICA-TS,ICA-TS-ICA]
%                case TS: [TS-ICA]
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
if nargin<4; error('bss_step: wrong number of input arguments \n'); end;
if size(ecg,2)<size(ecg,1); ecg = ecg'; end;
if size(icaecgf,2)>size(icaecgf,1); icaecgf = icaecgf'; end;
if nargin<5; debug=0; end;

% == general
NB_ABD_CHANNELS = size(icaecgf,2);
MQRS = cell(NB_ABD_CHANNELS,1);

% == core function
try
    if isempty(find(strcmp(paramStruct.method,{'DEFLATION','FUSE'}),1))
        % == ICA on residual in time domain
        [~,icaecgr] = jade(ecg);
        ecgout = icaecgr;
    end

    if find(strcmp(paramStruct.method,{'DEFLATION','FUSE'}))
        % == ICA on residual in source space
        recgs = zeros(size(icaecgf')); % residual ecg in source domain
        for i=1:NB_ABD_CHANNELS
            MQRS{i} = adjust_mqrs_location(icaecgf(:,i),rMQRS,paramStruct.fs,0);
            recgs(i,:) = mecg_cancellation(MQRS{i},icaecgf(:,i)','TS',20);
            %[recgs(i,:),~,~] = run_kalman(icaecgf(:,i),MQRS{i},paramStruct);
        end
        [~,icaecgr] = jade(recgs);
        if strcmp(paramStruct.methodComplexity,'COMPLEX')
            ecgout = [recgs;icaecgr]; 
        else
            ecgout = recgs;
        end
    end
    
    if size(ecg,1)>size(ecg,2)
        ecgout = ecgout';
    end
    
catch ME
    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
    ecgout = ecg;
end

% == plots
if debug
   ax(1) = subplot(211); plot(ecg); 
   title('ecgin');
   ax(2) = subplot(212); plot(ecgout);
   title('ecgout');
   linkaxes(ax,'x');
end

end

% SOME NOTES ON THIS FUNCTION:
% == STEP 1: direct BSS on filtered ecg (and not residuals)
    % Adding another step where BSS is performed on the filtered ABD signal was usefull 
    % for the EKF because sometimes the EKF kills the FECG and applying BSS after removing 
    % the MECG does not work really well.
% == STEP 2: ICA on residual in time domain
    % NOTE: jade performed better than PCA and FastICA (default options) when applied to the residual
    % signal but PCA performed better than ICA when applied to the initial
    % prefiltered mixture. Difference in performance was >1% for F1
    % measure but this is due to one specific record that is not separated
    % by ICA but is by ICA. Removing this specific record results in an
    % increase in F1 when using ICA rather than PCA on the ecgf
    % matrix. Having more data in the training set would have been usefull
    % to confirm this observation.
% == STEP 3: ICA on residual in source space    
    % ICA to move to the space domain, performing MECG cancellation in that
    % domain then take ICA of mother-free source signals. 
    % So it is sort of ICA of ICA...   

% == NOTE: in theory ICA would need centering that is subtracting the mean of the
% signal for each channel before applying ICA. Performing this operation
% lowered the F1 measure slightly!
% [~,icaecgr] = jade(bsxfun(@minus,recgs,mean(recgs,2)));
    
    
