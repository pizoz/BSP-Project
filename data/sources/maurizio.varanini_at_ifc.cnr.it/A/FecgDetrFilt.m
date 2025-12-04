function Xd=FecgDetrFilt(X,fs,cName,graph,dbFlag)
% -------------------------------------------------------------------------------------------
%   ECG detrending
%
%  A baseline signal is estimated applying a low pass Butterworth filter 
%  in forward and backward direction.  The filter cut frequency is 3.17Hz. 
%  Each detrended signal is obtained as difference between the original signal 
%  and the estimated baseline. 
%  In case of residual artifacts due to fast baseline movements, median filtering is used.
%
%  function Xd=FecgDetrFilt(X,fs,cName,graph,dbFlag)
%
%  X      : signal matrix (one signal per column)
%  fs     : sampling frequency
%  cName  : record name
%  graf   : flag enabling figure drawing
%  dbFlag : flag enabling figure drawing for tuning and debugging
%
% -------------------------------------------------------------------------------------------
%   Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
%   For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% -------------------------------------------------------------------------------------------
% This program is free software; you can redistribute it and/or modify it under the terms
% of the GNU General Public License as published by the Free Software Foundation; either
% version 2 of the License, or (at your option) any later version.
%
% This program is distributed "as is" and "as available" in the hope that it will be useful,
% but WITHOUT ANY WARRANTY of any kind; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% -------------------------------------------------------------------------------------------

if(nargin<3), cName=''; end
if(nargin<4), graph=0; end
if(nargin<5), dbFlag=0; end
graphD= graph &dbFlag ;
graphSpt= graph &dbFlag ;

fprintf('\n --------------------------------------------------------- \n');
[progpath, progname] = fileparts(which(mfilename));
fprintf('Program: %s,   record name: %s\n', progname, cName);
%-------------------------------------------------------------
% recording duration
[ndt, ns]=size(X);
vtime= [1:ndt]/fs;
%-------------------------------------------------------------------------------
fmaxd=5;        % --> 3.17 Hz = filtfilt 3dB cut frequency
fmaxn = fmaxd/(fs/2);
[b,a]= butter(1,fmaxn,'low');
Xb=filtfilt(b,a,X);  %  estimated baseline

%  remove baseline from original ecg
Xd=X-Xb;
%-------------------------------------------------------------------------------
RRpmean=1.;
ww=fix(RRpmean*1.5*fs);   % wide windows  (many contain at least one QRS)
for is=1:ns
    x=Xd(:,is); n=length(x);
    [xmeaMi,xmeaMa] = meanMiMaSc(x,ww,8,8);
    thmima=1.3;
    xthmi=xmeaMi*thmima;    xthma=xmeaMa*thmima;
    if(graphD)
        figure, set(gcf,'Color','white');  hold on;
        plot(x,'b.-');
        plot(xmeaMi*ones(n,1),'m');
        plot(xmeaMa*ones(n,1),'m');
        plot(xthmi*ones(n,1),'r');
        plot(xthma*ones(n,1),'r');
        title(['EcgBaseline_atest:',' x(b), mi&ma (m), th (r)']);
    end
    
    iibadDetr=find(x<xthmi | xthma<x);
    % iibadDetr(find(diff(iibadDetr)>1*fs));
    if(length(iibadDetr)>1.3*fs)
        fprintf('Switching to median filtering,rec=%s, is=%d, bads=%d\n',cName, is, length(iibadDetr));
        wm=fix(0.26*fs); wm=2*floor(wm/2)+1;    %  window width
        Xb(:,is) = medfilt1mit(X(:,is),wm,1);
        x=X(:,is)-Xb(:,is);         % remove baseline from original ECG
        [xmeaMi,xmeaMa] = meanMiMaSc(x,ww,8,8);
        thmima=1.3;
        xthmi=xmeaMi*thmima;    xthma=xmeaMa*thmima;
        iinf=find(x<xthmi); x(iinf)=xthmi;
        isup=find(x>xthma); x(isup)=xthma;
        fmaxd=150;    fmaxn=fmaxd/(fs/2);  % Normalized cut-off frequency
        [b,a]= butter(1,fmaxn,'low');
        x=filtfilt(b,a,x);     %  low-pass filtered signal
        Xd(:,is)=x;
        if(graphSpt)
            figure; freqz(conv(b,b(end:-1:1)),conv(a,a(end:-1:1)),[0:0.005:2],fs);
        end
    end
end

if(graphSpt)
    for is=1:ns,
        figure, set(gcf,'Color','white');
        pwelch(X(:,is),[],[],[],fs);
        title('Original ECG Welch spectrum');
        shg
    end
end
if(graphSpt)
    for is=1:ns,
        figure, set(gcf,'Color','white');
        pwelch(Xd(:,is),[],[],[],fs);
        title('Detrended ECG Welch spectrum');
        shg
    end
end

if(graphD)
    figure, set(gcf,'Color','white');
    for is=1:ns,
        subplot(ns,1,is), hold on, plot(vtime,X(:,is),'r'), plot(vtime,Xb(:,is),'b');
        wgmimaV=mimaxscG(X(:,is),0,0,.1);
        ylim(wgmimaV);
        if(is~=ns), set(gca,'XTickLabel',''); end
        if(is==1), title('Original ECG (r) & estimated baseline(b)'); end
    end
    shg
end
if(graph)
    figure, set(gcf,'Color','white');
    for is=1:ns,
        subplot(ns,1,is), hold on, plot(vtime,X(:,is),'r'), plot(vtime,Xd(:,is),'b');
        wgmimaV=mimaxscG(Xd(:,is),0,0,.4);
        ylim(wgmimaV);
        if(is~=ns), set(gca,'XTickLabel',''); end
        if(is==1), title([cName,': original(r) & detrended(b) ECG']); end
    end
end
%if(dbFlag & learning) PlotSgnMrkNc(X', QRSa, fs, cName); end

end %== function ================================================================
%

