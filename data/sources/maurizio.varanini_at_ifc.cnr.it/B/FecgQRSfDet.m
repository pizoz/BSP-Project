function qrsFes=FecgQRSfDet(X,fs,cName,qrsM,graph,dbFlag,saveFig,saveFigRRf,qrsAfs)
%---------------------------------------------------------------------------------------------
%   Fecg: Fetal QRS detection
%   - For each ECG channel:
%     - Fetal QRS pre-detection
%     - Identification of a good interval
%     - Fetal QRS detection in forward direction starting from the beginning
%       and in backward starting from the end of such a interval
%   - Best fetal QRS detection selection based on fetal RR variability and
%     fetal-mother QRS position coincidence
% --------------------------------------------------------------------------------------------
% Author: Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
% For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% --------------------------------------------------------------------------------------------
% This program is free software; you can redistribute it and/or modify it under the terms
% of the GNU General Public License as published by the Free Software Foundation; either
% version 2 of the License, or (at your option) any later version.
%
% This program is distributed "as is" and "as available" in the hope that it will be useful,
% but WITHOUT ANY WARRANTY of any kind; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% --------------------------------------------------------------------------------------------
if(nargin<3), cName=''; end
if(nargin<4), qrsM=[]; end
if(nargin<5), graph=1; end
if(nargin<6), dbFlag=0; end
if(nargin<7), saveFig=0; end
if(nargin<8), saveFigRRf=0; end
if(nargin<9), qrsAfs=0; end

fprintf('\n --------------------------------------------------------- \n');
[progpath, progname] = fileparts(which(mfilename));
fprintf('Program: %s,   record name: %s\n', progname, cName);

if(saveFig), figFmt='png';
    figPath=fullfile('../Figure/',progname);
    if(~exist(figPath,'dir')), mkdir(figPath); end
end

if(~isempty(qrsAfs));
    learning=1;   % learning set, annotations are available
else
    learning=0;  %test set
end

%-------------------------------------------------------------

% recording duration
[ndt, ns]=size(X);

vtime= [1:ndt]/fs;
if(dbFlag && learning), PlotSgnMrkNc(X, qrsAfs*fs, fs, [cName,': faQRS']); end

% -----
if(dbFlag), PlotSgnMrkNc(X, qrsM, fs, [cName,': mQRS']); end
qrsMs=qrsM/fs;
%--------------------------------------------------------------------------------------------
for is=1:ns,
    X(:,is)= (X(:,is)-mean(X(:,is)))/std(X(:,is));
end
% X(:,ns+1)= mean(X')';  % if ICA fails S/N improvment could be obtained
% ns=ns+1;

if(graph)
    figure, set(gcf,'Color','white');
    for is=1:ns,
        subplot(ns,1,is), plot(vtime,X(:,is));
        wgmi1= min(X(:,is)) -2;
        wgma1= max(X(:,is)) +2;
        ylim([wgmi1, wgma1]);
        set(gca,'YTick',[-5 0 5])
        % if(is~=ns), set(gca,'XTickLabel',''); end
        if(is==1), title([cName,': ICA fetal ECG']); end
    end
end

ecg=X;

% raw derivative filter coefficients
% 5ms before and after, 3ms in the middle
nu=ceil(0.005 * fs); nz=floor(0.0030*fs /2)*2 +1;  % nz=nearest odd value
B=[ones(nu,1);zeros(nz,1);-ones(nu,1)];
delay=floor(length(B)/2);
% ---   compute the derivative signal
ecgfx=[repmat(ecg(1,:),delay,1);ecg;repmat(ecg(end,:),delay,1)];
decgr=filter(B,1,ecgfx);   decgr= decgr(2*delay+1:end,:);
if(dbFlag)
    figure; freqz(B,1,1000,fs);
end
% ----- Butterworth forward and backward bandpass filtered (1-6Hz)
fmind=0.7; fmaxd=8;
Wn = [fmind, fmaxd]/(fs/2);  % normalized cut-off frequency (0,1)
[b,a]= butter(1,Wn);
decgBp=filtfilt(b,a,decgr);
% -----   absolute derivative ------
adecg=abs(decgr);
if(graph)
    figure, set(gcf,'Color','white');
    wgmima= mimaxscG(adecg,0,0,0.1);
    for is=1:ns,
        subplot(ns,1,is), plot(vtime,adecg(:,is));
        ylim(wgmima);
        % if(is~=ns), set(gca,'XTickLabel',''); end
        if(is==1), title([cName,': ECG absolute derivative']); end
    end
end
% ----- Butterworth forward and backward bandpass filtered (1-6Hz)
fmind=0.7; fmaxd=8;
Wn = [fmind, fmaxd]/(fs/2);  % normalized cut-off frequency (0,1)
[b,a]= butter(1,Wn);
adecgBp=filtfilt(b,a,adecg);

if(graph)
    figure, set(gcf,'Color','white');
    wgmima= mimaxscG(adecgBp(fs:end-fs),0,0,0.1);
    for is=1:ns,
        subplot(ns,1,is), plot(vtime,adecgBp(:,is));
        ylim(wgmima);
        % if(is~=ns), set(gca,'XTickLabel',''); end
        if(is==1), title([cName,': ECG filtered absolute derivative']); end
    end
end

% ----
if(learning),  nQRSaf=length(qrsAfs);
    fprintf('Number of annotated fetal QRSs= %d\n', nQRSaf);
    RRafs= diff(qrsAfs);          % annotated fetal RR
end
RRmc= diff(qrsM); RRms= RRmc/fs;  % detected mother RR
%
% -------------------------------------------------------------------------------
qrsFs=cell(1,ns);
RRfs=cell(1,ns);

for ic=1:ns,
    
    % ----
    xic=X(:,ic);
    ecgic=ecg(:,ic);
    decgic=decgr(:,ic);
    decgBpic=decgBp(:,ic);
    adecgic=adecg(:,ic);
    adecgBpic=adecgBp(:,ic);
    
    if(dbFlag && learning), PlotSgnMrkNc(xic, qrsAfs*fs, fs, [cName,' - eS']); end
    
    % ----  QRS detection   (first step) ---------------------------------------------
    %        qrsM=QRSdetectorQeA2T(vadx,vdx,fs,pth,RRts,pmQT)
    pth=0.35; qrsFp1=QRSdetectorF1(adecgBpic,decgic,fs,pth,0.45,0.8);
    %
    if(dbFlag), fprintf('Number of pre-detected fetal QRSs= %d\n', length(qrsFp1)); end
    RRfc1= diff(qrsFp1); RRfs1= RRfc1/fs;  % predetected fetal RR
    
    if(graph)
        figure, set(gcf,'Color','white'); hold on;
        plot(qrsMs, [RRms(1);RRms],'m.-');
        if(learning),   plot(qrsAfs, [RRafs(1);RRafs],'b.-');  end
        plot(qrsFp1/fs, [RRfs1(1);RRfs1],'r.-'); ylim([0.1, 1]);
        if(learning), title([cName,': detected mother(m), annotated (b) & pre-detected(r) fetal RR series']);
        else title([cName,': detected mother (m) & pre-detected fetal(r) RR series']); end
        figResize(0, 1, 1, .35);
    end
    
    %  ---- Identification of a "good" RR interval  ----

    % [nb,vp]=hist(RRfs,[0.25:0.05:0.9]);
    RRfsc1=sort(RRfs1); RRfsc1=RRfsc1(round(0.05*length(RRfsc1)):round(0.95*length(RRfsc1)));
    [nb,vp]=hist(RRfsc1);
    if(graph)
        figure, bar(vp,nb,1); title([cName,': pre-detected RR hist']);
    end
    [nbm,ibm]=max(nb);  RRst=vp(ibm);   % most frequent value
    
    adRR=abs(diff(RRfs1));
    adRRx=[repmat(adRR(1),2,1);adRR;repmat(adRR(end),2,1)];
    %adRRx=medfilt1(adRRx,5);
    mvadRR=filter(ones(1,5),1,adRRx); mvadRR= mvadRR(2*2+1:end,:);
    if(graph), figure, plot(qrsFp1(3:end)/fs,mvadRR); end
    thad=0.4;
    RRoklen=0;
    while(~RRoklen)
        for j=1:4
            isRRok=[]; ieRRok=[];
            jint=1;  RRokFlag=0;
            for i=1:length(mvadRR)
                if(RRst-0.1<RRfs1(i) && RRfs1(i)<RRst+0.1 && mvadRR(i)<thad)
                    if(~RRokFlag), RRokFlag=1; isRRok(jint)=i; else ieRRok(jint)=i; end
                else
                    if(RRokFlag), ieRRok(jint)=i; jint=jint+1; RRokFlag=0; end
                end
            end
            if(RRokFlag), ieRRok(jint)=i+1; RRokFlag=0; end
            if(isempty(isRRok)), RRoklen=0;
            else [RRoklen,im]=max(ieRRok-isRRok);
            end
            if(RRoklen>5), break; end
            thad=thad+0.1;
        end
    end
    RRsmi=mean(RRfs1(isRRok(im):ieRRok(im)));
    iqis=qrsFp1(isRRok(im)); iqie=qrsFp1(ieRRok(im));
    sgnOklens=(iqie-iqis)/fs;
    
%   --- Fetal QRS detection in backward direction starting from the end 
%       of the identified "good" RR interval
    pth=0.2; qrsFb=QRSdetectorF2(adecgBpic(iqie:-1:1),decgic(iqie:-1:1),fs,pth,RRst,RRsmi,sgnOklens,1);  % backward
    qrsFb=iqie-qrsFb(end:-1:1);
    
%   --- Fetal QRS detection in forward direction starting from the beginning 
%       of the identified "good"  RR interval
    pth=0.2; qrsFf=QRSdetectorF2(adecgBpic(iqis:end),decgic(iqis:end),fs,pth,RRst,RRsmi,sgnOklens,0);    % forward
    qrsFf=iqis+qrsFf;
    
    if(isempty(qrsFb))
        qrsFic=qrsFf;
    elseif(isempty(qrsFf))
        qrsFic=qrsFb;
    else
        ibe=min(find(qrsFb>qrsFf(1)-0.5*RRst*fs));
        ifs=min(find(qrsFf>qrsFb(ibe)+0.5*RRst*fs));
        qrsFic=[qrsFb(1:ibe); qrsFf(ifs:end)];
    end
    
    nQRSdf=length(qrsFic);
    if(dbFlag), fprintf('Channel=%d, number of detected fetal QRSs= %d\n', ic, nQRSdf); end
    qrsFsic=qrsFic/fs;
    RRfsic= diff(qrsFsic);   % detected fetal RR
    [nb,vp]=hist(RRfsic,[0.25:0.05:0.9]);
    if(graph)
        figure, bar(vp,nb,1);  title([cName,': detected RR hist']);
    end
    [nbm,ibm]=max(nb);  RRst=vp(ibm);
    
    RRfmean=meansc(RRfsic,4,4);
    %RfRmean= mean(RRfsic);
    RRfstd= std(RRfsic);
    if(dbFlag), fprintf('fetal RR mean= %d,  stdev=%d\n', RRfmean, RRfstd); end
    debugQRS=1;
    if(dbFlag && debugQRS)
        figure
        subplot(2,1,1),plot(adecgic,'b'); title([cName,': filter absolute derivative']);
        DrawVertMarker(qrsFic,'r',':','none');
        %     DrawVertMarker(qrsF(2,:)','b',':','none');
        %     DrawVertMarker(qrsF(end,:)','r',':','none');
        subplot(2,1,2),plot(ecgic(:,1),'b'); title('ecg signal');
        DrawVertMarker(qrsFic,'r',':','none');
        %     DrawVertMarker(qrsF(2,:)','b',':','none');
        %     DrawVertMarker(qrsF(end,:)','r',':','none');
    end
    
    if(graph || saveFig)
        if(graph), figure; else figure('visible','off'); end
        set(gcf,'Color','white'); hold on;
        plot(qrsMs, [RRms(1);RRms],'m.-');
        if(learning),  plot(qrsAfs, [RRafs(1);RRafs],'b.-');  end
        plot(qrsFsic, [RRfsic(1);RRfsic],'r.-');
        ylim([0.1, 1]);
        if(learning), title([cName,': ic=',num2str(ic),', detected mother(m), annotated (b) & detected(r) fetal RR series']);
        else title([cName,': ic=',num2str(ic),', detected mother (m) & fetal(r) RR series']); end
        figResize(0, 1, 1, .35);
        if(saveFigRRf), figFmtRR='png';
            figPathRR='../Figure/';
            if(~exist(figPathRR,'dir')), mkdir(figPathRR); end
            figName=fullfile(figPathRR,[cName,'_',num2str(ic),'_RRfetal']);
            print(gcf, ['-d',figFmtRR],figName);
        end
    end
    
    qrsFs{ic}=qrsFsic;
    RRfs{ic}=RRfsic;
    
end

% -----------------------------------------------------------------------------
% --- Choosing the best fetal RR series and QRS detection
%     based on fetal RR variability and fetal-mother QRS position coincidence
cfact=zeros(ns,1);  meRR=zeros(ns,1);
madRR = zeros(ns,1);  maddRR=zeros(ns,1);  % rmssd = zeros(ns,1);
cFMsim=zeros(ns,1);   chanIqf=zeros(ns,1);
for ic=1:ns,
    meRR(ic) = median(RRfs{ic});
    if(meRR(ic)<0.35), cfact(ic)=0.3*(0.35-meRR(ic));
    elseif(meRR(ic)>0.5), cfact(ic)=0.05*(meRR(ic)-0.5);
    else cfact(ic)=0;
    end
    dRRfs=diff(RRfs{ic}); ddRRfs=diff(dRRfs);
    madRR(ic) = meansc(abs(dRRfs),0,5);
    maddRR(ic) = meansc(abs(ddRRfs),0,5);
    %rmssd(ic) = std(diff(RRfs{ic}));
    %dRRfsz=dRRfs-mean(dRRfs); rmssd(ic) = sqrt(meansc(dRRfsz.*dRRfsz,0,8));
    cmpRes=QRSdet_Ann_cmp(qrsFs{ic},qrsMs, 0.1, 0);
    %FMsim(ic)=cmpRes.nTP/(cmpRes.nFP+cmpRes.nFN);
    cFMsim(ic)=0.05*cmpRes.nTP/(0.05*cmpRes.nTP+cmpRes.nFP+cmpRes.nFN);
    fprintf('channel=%1d,  meRR=%7.4f,  madRR=%7.4f,  maddRR=%7.4f, cfact=%7.4f, cFMsim=%7.4f\n', ...
        ic, meRR(ic), madRR(ic),maddRR(ic),cfact(ic),cFMsim(ic));
    chanIqf(ic)=madRR(ic)+maddRR(ic)+cfact(ic)+cFMsim(ic);
end
for ic=1:ns,
    fprintf('channel=%1d,  chanIqf=%7.4f\n', ic, chanIqf(ic));
end
[dummy,icsV]=sort(chanIqf);  % select the best channel
ics=icsV(1);
fprintf('Selected channel=%1d\n', ics);

qrsFes=qrsFs{ics};
RRfes=RRfs{ics};

% xs=X(:,ics);
ecgs=ecg(:,ics);
% decgs=decgr(:,ics);
% adecgs=adecg(:,ics);

if(graph)
    figure, set(gcf,'Color','white');
    for is=1:ns,
        subplot(ns,1,is), plot(vtime, X(:,is));
        wgmi1= min(X(:,is)) -2;
        wgma1= max(X(:,is)) +2;
        ylim([wgmi1, wgma1]);
        set(gca,'YTick',[-5 0 5])
        % if(is~=ns), set(gca,'XTickLabel',''); end
        if(is==1), title([cName,': sorted fetal ECG']); end
        shg
    end
end
% -----------------------------------------------------------------------------
% if(dbFlag) PlotSgnMrkNc(ecgs, qrsFes, fs, cName);  end

if(learning && graph), PlotSgnMrkNc(ecgs, {qrsAfs*fs,qrsFes}, fs, [cName,': fecg, ann.(m) & det.(b) qrs']); end
if(saveFig)
    figName=fullfile(figPath,[cName,'_ECGfetal']);
    print(gcf, ['-d',figFmt],figName);
end

if(graph || saveFig)
    figure, set(gcf,'Color','white');
    for is=1:ns,
        subplot(ns,1,is), plot(vtime, X(:,is));
        wgmi1= min(X(:,is)) -2;
        wgma1= max(X(:,is)) +2;
        ylim([wgmi1, wgma1]);
        %   DrawVertLines(vtime', qrsF(1,:)/fs, 'r', ylim);
        if(learning), DrawVertMarker(qrsAfs,'k',':','+'); end
        DrawVertMarker(qrsFes,'r',':','none');
        set(gca,'YTick',[-5 0 5])
        % if(is~=ns), set(gca,'XTickLabel',''); end
        if(is==1)
            if(learning), title([cName,': annotated (k) detected (r) fetal qrs']);
            else title([cName,': detected (r) fetal qrs']); end
        end
    end
    shg
    if(saveFig)
        figName=fullfile(figPath,[cName,'_fetal']);
        print(gcf, ['-d',figFmt],figName);
    end
end
%-------------------------------------------------------------

if(graph || saveFigRRf)
    if(graph), figure; else figure('visible','off'); end
    set(gcf,'Color','white'); hold on;
    plot(qrsMs, [RRms(1);RRms],'m.-');
    if(learning),  plot(qrsAfs, [RRafs(1);RRafs],'b.-');  end
    plot(qrsFes, [RRfes(1);RRfes],'r.-');
    ylim([0.1, 1]);
    if(learning), title([cName,': detected mother(m), annotated (b) & detected(r) fetal RR series']);
    else title([cName,': detected mother (m) & fetal(r) RR series']); end
    figResize(0, 1, 1, .35);
    if(saveFigRRf), figFmtRR='png';
        figPathRR='../FigureRRf/';
        if(~exist(figPathRR,'dir')), mkdir(figPathRR); end
        figName=fullfile(figPathRR,[cName,'_RRfetal']);
        print(gcf, ['-d',figFmtRR],figName);
    end
end

end %== function ================================================================

