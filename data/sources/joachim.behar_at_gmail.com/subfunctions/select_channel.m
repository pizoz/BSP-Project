function [sQRS,chNb] = select_channel(ecg,ecgType,varargin)
% this function automates the channel selection onto
% which to perform sQRS detection for whether to get MQRS or FQRS. 
% In the case of FQRS extraction, the key ideas are to to use the regularity 
% of the extracted HR to decide what FECG channel to choose and to use bxb 
% (proportion of matching sQRS btw two sQRS time series) to check that we 
% are not picking up the MECG.
%
% inputs
%   ecg:                the preprocessed FECG matrix
%   ecgType:            'FECG' or 'MECG'
%   varargin
%       MQRS:           MQRS position (in ms)
%       paramStruct:    structure containing challenge parameters
% outputs
%   sQRS:               selected QRS location (in ms)
%   chNb:               selected ABD channel
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
MQRS = [];
paramStruct.detThres = 0.5;
paramStruct.fs = 1000;
paramStruct.qrsdet_wind = 15;
debug = 0;

switch nargin
    case 3
        MQRS = varargin{1};
    case 4
        MQRS = varargin{1};
        paramStruct = varargin{2};
    case 5
        MQRS = varargin{1};
        paramStruct = varargin{2};
        debug = varargin{3};
    otherwise
        error('select_channel: wrong number of input arguments \n');
end

if size(ecg,2)>size(ecg,1)
    ecg = ecg';
end

NB_AECG = size(ecg,2); % nb of ABD ecg
QRSce = cell(NB_AECG,1); % sQRS cell
sQRS = cell(1,1);
SMI = zeros(NB_AECG,1); % smoothness indices
bxb = zeros(NB_AECG,1); % matching btw MQRS and another QRS time series

try
    % == core function
    if strcmp(ecgType,'FECG')
        % == FECG channel selection
        for ch=1:NB_AECG
            % Normalize residuals to avoid big bumps that crashed the FQRS detector (this showed improvement)
            [ecgn,~] = normalise_ecg(ecg(:,ch),1,5*paramStruct.fs);

            % sQRS detection
            QRSce{ch} = run_qrsdet_by_seg(ecgn,paramStruct.fs,paramStruct.qrsdet_wind,paramStruct.detThres,ecgType);
            
            % smooth the time series (physiological)
            if paramStruct.smoothRR
                QRSce{ch} = smooth_rr(QRSce{ch});
            end
            
            % compute smoothness indice
            [SMI(ch),~] = assess_regularity(QRSce{ch},1);

            % check that the FQRSs are not matching the MQRSs
            bxb(ch) = bsqi(MQRS/paramStruct.fs,QRSce{ch}/paramStruct.fs);
        end

        % if time series too similar to the one from the mother then we are
        % probably picking up the mother residuals. Then put SMI to inf to
        % not select this channel
        SMI(bxb>0.4) = Inf;

        % in order to select the ABD channel onto which the FQRS will be
        % identified we look at SMI (lower SMI -> most likely to be a good detection)
        [~,ind_std_trans] = sort(SMI);
        chNb = ind_std_trans(1);
        sQRS = QRSce{chNb};

        % == plots
        if debug || isempty(nargout)
            tm = 1/paramStruct.fs:1/fs:length(ecg)/paramStruct.fs;
            [m,n] = size(ecg);
            ax = zeros(m,1);
            for pp=1:n
                ax(pp) = subplot(ceil(n/2),2,pp); plot(tm,ecg(:,pp)); 
                legend(strcat('SMI:',num2str(SMI(pp)),' bxb:',num2str(bxb(pp))));
                linkaxes(ax,'x');
            end
        end

        
        
        elseif strcmp(ecgType,'MECG')
             % == MECG channel selection
             % this step is performed to get more robust MQRS detection.
             % This is aimed at finding the MQRS with the best accuracy but
             % ALSO to distinguish between MQRS and FQRS time series since
             % in some instances the FECG is of higher or similar amplitude 
             % that the MECG. In practice a chest electrode should be used
             % to avoid this to happen as the chest signal will be of
             % better quality and will ensure that there is no confusion
             % between MQRS and FQRS time series.
             [~,icaECG] = jade(ecg); 
             ecg = [ecg';icaECG]';
             NB_AECG = 2*NB_AECG; 
             lengthQRS  = zeros(4,1);

             for i=1:NB_AECG
                % run sQRS detector
                sQRS{i} = run_qrsdet_by_seg(ecg(:,i),paramStruct.fs,10,0.6,ecgType);
                lengthQRS(i) = length(sQRS{i});
                
                % compute SMI of MHR at 96% CI
                [SMI(i),~] = assess_regularity(sQRS{i});
             end

             if isempty(find(SMI~=inf,1))
                [~,chNb] = min(lengthQRS);
                sQRS = sQRS{chNb};
             else
                % this set of rules is specific to the PCinC2013
                % challenge and should not be kept for practical
                % applications
                indLowSMI = find(SMI<10);
                if length(indLowSMI)>1
                    [mi,indMin] = min(lengthQRS(indLowSMI));
                    [ma,indMax] = max(lengthQRS(indLowSMI));
                    if ma-mi>7
                        % if there is no big gap in length btw the 'good' 
                        % time series then take the most regular
                        chNb = indLowSMI(indMin);
                    else
                        % if there is the gap then we might be detecting
                        % both the MQRS and FQRS in which case it is safer
                        % to assume that the MQRS is the shortest time
                        % series so take the longuest for the foetus
                        chNb = indLowSMI(indMax);
                    end
                else
                    [~,chNb] = min(SMI);
                end
                sQRS = sQRS{chNb};
             end  
    else          
            disp('Wrong entry \n');
    end

catch ME
    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
    sQRS = [20 400 600 1000]; chNb=1;
end

end


% alternative to bSQI is to use the Physionet bxb function but this requires converting 
% .mat annotations to WFDB format. This is unfortunately taking a lot of computational time!!
% (0.116627sec)
%
%         mat2wfdbanns('qrs','fqrs',QRSce{ch},paramStruct.fs);
%         mat2wfdbanns('qrs','mqrs',MQRS,paramStruct.fs);
%         bxb('qrs','mqrs','fqrs','bxbReport.txt','00:00:00','00:01:00');
%         testFil = importdata('bxbReport.txt');
%         delete('bxbReport.txt');
%         bxb(ch) = testFil.data(1,1)/(testFil.data(1,1)+testFil.data(6,1));

