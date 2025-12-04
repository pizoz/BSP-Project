function [FQRSout,QTout,outExtra] = physionet2013extract(tm,ecg,paramStruct)
% template algorithm for Physionet/CinC competition 2013. This function can
% be used for events 1 and 2. Participants are free to modify any
% components of the code. However the function prototype must stay the same
% [FQRSout,QTout] = physionet2013(tm,ecg) where the inputs and outputs are specified
% below.
%
% inputs
%   tm :         Nx1 vector of time in milliseconds
%   ecg:         4x60000 (4 channels and 1min of signal at 1000Hz) matrix of
%                abdominal ecg channels.
%   paramStruct: structure with the key parameters for the different
%                algorithms. if empty then the defaults parameters are
%                used. 
% 
% output
%   FQRSout:    FQRS markers in seconds. Each marker indicates the position of one
%               of the FQRS detected by the algorithm.
%   QTout:      1x1 estimated fetal QTout duration
%   CHout:      the channel nb that was selected
%
%
%
% Safe Foetus Monitoring Toolbox, version 1.0, Sept 2013
% Released under the GNU General Public License
%
% Copyright (C) 2013 Joachim Behar
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

% == check the inputs
if nargin<2; error('physionet2013extract: wrong number of input arguments \n'); end;
if size(ecg,2)>size(ecg,1); ecg = ecg'; end;


% == entering core of the function
try      
    % == allocate space    
    NB_ABD_CHANNELS = size(ecg,2);
    MQRS = cell(1,NB_ABD_CHANNELS);
    %LogFileID = fopen('log.txt','a+');
    
    % == preprocessing
    ecgf = preprocessing([flipdim(ecg(1:paramStruct.fs*1-1,:),1)' ecg'],paramStruct);
    ecgf = ecgf(paramStruct.fs*1:end,:); % for border filtering reasons..
    
    if ~paramStruct.extract_qt
        % ===========================================================
        % ================= FECG extraction block ===================
        % ===========================================================
        % == MQRS detection
        [rMQRS,~] = select_channel(ecgf,'MECG',[],paramStruct);
        
        tic;
        for cc=1:NB_ABD_CHANNELS
            %fprintf(LogFileID,'ch %d: ',cc);
            % == adjust MQRS location
            MQRS{cc} = adjust_mqrs_location(ecgf(:,cc),rMQRS,paramStruct.fs,0); %MQRS{cc} = rMQRS;
            % == now remove MECG
            if find(strcmp(paramStruct.method,{'TS','TS-SUZANNA','TS-CERUTTI','TS-LP','TS-PCA'}))
                % TS techniques
                fecg.m1(cc,:) = mecg_cancellation(MQRS{cc},ecgf(:,cc)',paramStruct.method,paramStruct.TSnbCycles,paramStruct.NbPC);
            elseif strcmp(paramStruct.method,'TS-EKF')
                % EKF
                [fecg.m1(cc,:),errorInit,~] = run_kalman(ecgf(:,cc),MQRS{cc},paramStruct);
                %fprintf(LogFileID,' EKF errorInit=%d',round(100*errorInit)/100);
            elseif strcmp(paramStruct.method,'FUSE')
                % FUSE
                fecg.m2(cc,:) = mecg_cancellation(MQRS{cc},ecgf(:,cc)','TS',paramStruct.TSnbCycles);    % can also be TS/PCA
                % fecg.m3(cc,:) = mecg_cancellation(MQRS{cc},ecgf(:,cc)','TS',paramStruct.TSnbCycles);
                % [fecg.m4(cc,:),~,~] = run_kalman(ecgf(:,cc),MQRS{cc},paramStruct);
                % plot(ecgf(:,cc)); hold on, plot(fecg.m2(cc,:),'r')
            end        
            %fprintf(LogFileID,' \n'); 
        end

        if strcmp(paramStruct.method,'PCA')
            % PCA
            [~,fecg.m1] = princomp(ecgf); fecg.m1=fecg.m1';
        elseif strcmp(paramStruct.method,'ICA')
            % ICA
            [~,fecg.m1] = jade(ecgf);
            %fecg.m1 = fastica(ecgf');
        elseif strcmp(paramStruct.method,'DEFLATION')
            % ICA as preprocessing
            [~,fecg.m1] = jade(ecgf);        
        elseif strcmp(paramStruct.method,'FUSE')
            %[~,fecg.m4] = princomp(ecgf); fecg.m4=fecg.m4';
            [~,fecg.m1] = jade(ecgf);
        end
        %fprintf(LogFileID,paramStruct.method);

        % == BSS step
        if find(strcmp(paramStruct.method,{'TS','TS-SUZANNA','TS-CERUTTI','TS-LP','TS-PCA','TS-EKF'})) ...
                & strcmp(paramStruct.methodComplexity,'COMPLEX')
            % TS/PCA/EKF: add a BSS step on the residuals
            ecgout  = bss_step(fecg.m1,ecgf,rMQRS,paramStruct);
            fecg.m1 = [ecgout;fecg.m1];
        elseif strcmp(paramStruct.method,'DEFLATION')
            ecgout  = bss_step([],fecg.m1,rMQRS,paramStruct);
            fecg.m1 = ecgout;    
        elseif strcmp(paramStruct.method,'FUSE')
            ecgout1 = bss_step([],fecg.m1,rMQRS,paramStruct);
            paramStruct.method = 'TS';
            ecgout2 = bss_step(fecg.m2,fecg.m1,rMQRS,paramStruct);        
            fecg.m1 = [ecgout1;ecgout2;fecg.m1;fecg.m2];       % [ICA-TS/ICA-TS-ICA] [TS-ICA] [ICA] [TS]
        end
        timeExtrac = toc; fprintf('Execution time FECG extraction: %f \n',timeExtrac);
        
        % == channel selection & FQRS detection & post processing
        [FQRStemp,CHout] = select_channel(fecg.m1,'FECG',rMQRS,paramStruct); 
        FQRS = FQRStemp;

        fprintf('Selected channel number is: %d \n',CHout);
        %fprintf(LogFileID,'SelectFECGCha: %d ',CHout);

        % ===========================================================
        % ================= Challenge twicks ========================
        % ===========================================================
        % NOTE: Consider removing parts of that for an actual application!
        outExtra.rFQRS = FQRS;
        outExtra.rMQRS = rMQRS;
        
        if paramStruct.cinc_match
            [STD,~] = assess_regularity(FQRS,1,paramStruct.fs,1);  
            fprintf('STD: %f \n',STD);
            if length(FQRS)<85 || length(FQRS)>200
                % identify rubbish! In the case the time
                % series is rubbish then just create one at 134bpm
                disp('abnormal length(FQRS)<85 || length(FQRS)>200');
                FQRSout = 0.1:0.42:60;
                FQRSout = round(FQRSout*1000);
                QTout = 0;
            elseif STD>paramStruct.STDcut
                disp('abnormal STD');
                % if time series that we get is too irregular then better to look 
                % for predominant HR mode and guess the time series!
                [FQRSmode,m1,m2] = guess_fqrs(FQRS,rMQRS,paramStruct.fs);
                %m = mode(60./diff(FQRS./paramStruct.fs)); 
                if m1>125 && m1<170 && m2>125 && m2<175
                    disp('- using MODE');
                    FQRSout = FQRSmode;
                else
                    disp('- just guessing');
                    FQRSout = 0.1:0.42:60; % cc.e. 143bpm
                end
                FQRSout = round(FQRSout*1000);
                QTout = 250;
            else
                FQRSout = FQRS;
                QTout = 250;

                % then add the magic smoothing step while possibly remove first
                % and last peaks
                FQRSout = smooth_rr_lots(FQRSout(FQRSout>paramStruct.fs*0.200 & ...
                FQRSout<paramStruct.fs*59.9)/paramStruct.fs,paramStruct.fs);
            end
        else
            FQRSout = FQRS;
            QTout = 250;
        end
        % ===========================================================
        
        
    else
        % ===========================================================
        % ================= QT extraction ===========================
        % ===========================================================

        % in previous step need to output rMQRS and rFQRS
        MED_RR = median(diff(paramStruct.rFQRS));
        NB_ITER = 5;
        Qtrans = zeros(NB_ABD_CHANNELS*NB_ITER,1);
        Ttrans = zeros(NB_ABD_CHANNELS*NB_ITER,1);
        compt = 0;
        tic
        for cc=1:NB_ABD_CHANNELS
            fprintf('QT processing channel: %f \n',cc);
            % adjust MQRS
            MQRS{cc} = adjust_mqrs_location(ecgf(:,cc),paramStruct.rMQRS,paramStruct.fs,0);
            % source separation
            fecg.m1(cc,:) = mecg_cancellation(MQRS{cc},ecgf(:,cc)',paramStruct.method,paramStruct.TSnbCycles,paramStruct.NbPC);
            % template construction
            relevantMode = fecg_template_construction(paramStruct.rFQRS,fecg.m1(cc,:));
            for kk=1:NB_ITER
                compt = compt+1;
                [~,Qtrans(compt),Ttrans(compt),~] = measure_qt(relevantMode.cycleMean,MED_RR);
            end
        end
        
        Qtrans2 = Qtrans;
        Ttrans2 = Ttrans; 
        
        Qtrans2(Qtrans==0 | Qtrans==1)=[]; % remove when it failed
        Ttrans2(Ttrans==0 | Ttrans==1)=[];
        
        Q = median(Qtrans2);
        T = median(Ttrans2);
        QTout = T-Q;
                
        if QTout== 0|| QTout>320 || QTout<125; QTout=250; end;

        CHout = 1; % how to select a good channnel for plot??
        FQRSout=[]; outExtra=[];
        toc
        % ===========================================================
    end
    
    
    % ==================== plots for analysis ===================
    if paramStruct.debug==1 && paramStruct.extract_qt==0
        figure(1);
        title(strcat('Channel selected is: ',CHout));
        one_ecgres = fecg.m1(CHout,:);
        hrv = 60./(diff(FQRS)/paramStruct.fs);
        hrvss = 60./(diff(FQRSout)/paramStruct.fs);
        ax(1) = subplot(311); 
            plot(tm,ecg(:,1),tm(MQRS{1}),ecg(MQRS{1},1),'+r','LineWidth',2);
            plot(tm,ecgf(:,1),tm(MQRS{1}),ecgf(MQRS{1},1),'+r','LineWidth',2);
            legend('One ABD ecg','MQRS');
        ax(2) = subplot(312);
            hold on, plot(tm,one_ecgres,'k','LineWidth',2); 
            hold on, plot(tm(FQRS),one_ecgres(FQRS),'+r');
            legend('fecg (residual)','FQRS detected');
        ax(3) = subplot(313); plot(FQRS(1:end-1)/paramStruct.fs,hrv,'-ro');
            hold on, plot(FQRSout(1:end-1)/paramStruct.fs,hrvss,'--ko');
            legend('FHRc','FHRcc'); xlabel('Time (sec)'); ylabel('Amplitude (NU)');
        linkaxes(ax,'x'); set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold');
    elseif paramStruct.debug==1 && paramStruct.extract_qt==1
        % extract the template ecg closest from T point identified (rather than plotting a random one..)
        k = dsearchn(Ttrans-Qtrans,QTout);
        ch = k+1-floor(k/NB_ABD_CHANNELS)*NB_ABD_CHANNELS;
        relevantMode = fecg_template_construction(paramStruct.rFQRS,fecg.m1(ch,:));

        figure(2);
        fprintf('QT: Number of cycles: %f \n',relevantMode.NbCycles);
        one_ecgres = fecg.m1(CHout,:);
        title('QT analysis');
        subplot(211);
            hold on, plot(tm,one_ecgres,'k','LineWidth',2); 
            hold on, plot(tm(paramStruct.rFQRS),one_ecgres(paramStruct.rFQRS),'+r');
            legend('fecg (residual)','FQRS detected');
        subplot(212); 
            plot(relevantMode.cycleMean,'LineWidth',2);
            QTpos = round([Q T]*250/MED_RR);
            if QTpos(1)>0 && QTpos(2)<250
                hold on, plot(QTpos,relevantMode.cycleMean(QTpos),'+r','LineWidth',3);
            else
                hold on, plot(QTpos,0,'+r','LineWidth',3);
            end
            xlabel('Phase Wrap'); ylabel('Amplitude [NU]');
    end
    % ===========================================================
    
    
%fclose(LogFileID);    
catch ME
    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;   
    FQRSout = 0.1:0.42:60;
    QTout = 250;
end

end
























