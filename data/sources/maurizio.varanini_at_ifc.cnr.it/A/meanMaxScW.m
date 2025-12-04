function vme = meanMaxScW(v, wl, wm, stepwl, percmi,percma)
% ---------------------------------------------------------------------------------------------
% meanMaxScW.m: Compute the average value of the maxima on data windows of the input vector.
% ---------------------------------------------------------------------------------------------
% Compute the mean values of the maxima of signal values in subwindows of a sliding window
% The distribution tails can be excluded.
%   vme = meaMaxScW(v,wl,wm,stepwl,perci,percf)
%   v = signal vector
%   wl = width of the sliding window
%   wm = width of the subwindows
%   stepwl = sliding window step
%   perci = % of min values to be excluded
%   percf = % of max values to be excluded
%
% ---------------------------------------------------------------------------------------------
%   Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
%   For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% ---------------------------------------------------------------------------------------------

if(nargin<4), stepwl=1; end
if(nargin<5), percmi=5; end
if(nargin<6), percma=percmi; end
if nargin==0,  meanMaxScW_test; return; end

if(percmi<0 || percma<0 || (percmi+percma)>=100), fprintf('\n"meanMaxScW": error in input parameters\n'); pause; end
if(stepwl<1), stepwl=1; end
if(isempty(v)), vme=[]; return; end
[nt,dummy]=size(v);
if(nt<dummy), v=v'; nt=dummy; end

wlm=fix(wl/wm);
ii=1+floor(wlm*percmi/100);
fi=wlm-floor(wlm*percma/100);

wl=wlm*wm;
nwl=fix((nt-wl)/stepwl)+1;
vme=zeros(nwl,1);
iwl=1;
for is=1:nwl
    vmawm= max(reshape(v(iwl:iwl+wl-1),wm,wlm))';
    omawl= sort(vmawm);
    vme(is)=mean(omawl(ii:fi));
    iwl=iwl+stepwl;
end

return
end % = function ===============================================================

% -------------------------------------------------------------------------------------------------
% Test program (simulating ecg channel selection for QRS detection)
function meanMaxScW_test
close all, clear variables
Fs=1000; vlen=60*Fs;
RRpmean=0.4;
ww=fix(RRpmean*1.5*Fs);      % wide windows (each window contains at least one QRS)
wn=fix(RRpmean/5*Fs);      % short windows (many windows do not contain QRSs)
wl=20*ww;   %10*ww;
stepwl=0.2*wl;   %1; %0.5*wl; %wl; %fix(0.01*Fs);
fprintf('ww=%d,  wn=%d,  wl=%d, stepwl=%d\n', ww,wn,wl,stepwl);
% building of test signal
sigType=3;    % <== Test signal type
ns=2;
t=(1:vlen)'/Fs;
switch(sigType)
    case 1
        xb1 = abs(0*sin(t) + 3*randn(size(t)));
        xb2 = abs(0*sin(t) + 3*randn(size(t)));
        pp = floor(RRpmean*Fs);
        tp = 1:pp:length(xb1);
        xp=zeros(size(t)); xp(tp) = 3; %+1*abs(randn(size(tp)));
        xs(:,1)=xb1;
        xs(:,2)=0.3*xb2+xp;
    case 2
        triaNoise=ones(size(t)); nptr=1001; npCyc=1750; np1=npCyc-nptr;
        nptrt=floor(size(t,1)/(npCyc))*(npCyc);
        triaNoise(1:nptrt)=repmat([2+triang(nptr);ones(np1,1)],floor(size(t,1)/npCyc),1);
        xb1 = abs(.5*abs(sin(t))+3*triaNoise .*randn(size(t)));
        aecgd=ones(size(t)); npq=fix(0.1*Fs); npRR=floor(RRpmean*Fs); npic=npRR-npq;
        npecg=floor(size(t,1)/(npRR))*(npRR);
        aecgd(1:npecg)=repmat([2+triang(npq);ones(npic,1)],floor(size(t,1)/npRR),1);
        xb2 = abs(.5*abs(sin(t))+2*aecgd .*randn(size(t)));
        xs(:,1)=xb1;
        xs(:,2)=xb2;
    case 3
        triaNoise=ones(size(t)); nptr=350; npCyc=1350; np1=npCyc-nptr;
        nptrt=floor(size(t,1)/(npCyc))*(npCyc);
        triaNoise(1:nptrt)=repmat([2+triang(nptr);ones(np1,1)],floor(size(t,1)/npCyc),1);
        xb1 = abs(.5*abs(sin(t))+3*triaNoise .*randn(size(t)));
        aecgd=ones(size(t)); npq=fix(0.1*Fs); npRR=floor(RRpmean*Fs); npic=npRR-npq;
        npecg=floor(size(t,1)/(npRR))*(npRR);
        aecgd(1:npecg)=repmat([2+triang(npq);ones(npic,1)],floor(size(t,1)/npRR),1);
        xb2 = abs(.5*abs(sin(t))+2*aecgd .*randn(size(t)));
        xs(:,1)=xb1;
        xs(:,2)=xb2;
end
percmi=0;
percma=5;
for is=1:ns,
    vmew(:,is)=meanMaxScW(xs(:,is),wl,ww,stepwl,percmi,percma);
    vmen(:,is)=meanMaxScW(xs(:,is),wl,wn,stepwl,percmi,percma);
end

figure, set(gcf,'Color','white');
wgmima= mimaxscG(xs(Fs:end-Fs,1),0,0,0.1);
for is=1:ns,
    subplot(ns,1,is), plot(t,xs(:,is));
    ylim(wgmima);
    if(is~=ns), set(gca,'XTickLabel',''); end
    if(is==1), title('Signals'); end
end

vmtime=1/2*wl/Fs+(1:size(vmew,1))*stepwl/Fs; vmtime=vmtime';
figure, set(gcf,'Color','white');
wgmima= mimaxscG(vmew,0,0,0.1);
for is=1:ns,
    subplot(ns,1,is), plot(vmtime,vmew(:,is));
    ylim(wgmima);
    if(is~=ns), set(gca,'XTickLabel',''); end
    if(is==1), title('vmew'); end
end
figure, set(gcf,'Color','white');
wgmima= mimaxscG(vmen,0,0,0.1);
for is=1:ns,
    subplot(ns,1,is), plot(vmtime,vmen(:,is));
    ylim(wgmima);
    if(is~=ns), set(gca,'XTickLabel',''); end
    if(is==1), title(' vmen'); end
end
vqf=vmew./vmen;
%vqf=vmew-1.0*vmen;
figure, set(gcf,'Color','white');
wgmima= mimaxscG(vqf,0,0,0.1);
for is=1:ns,
    subplot(ns,1,is), plot(vmtime,vqf(:,is));
    ylim(wgmima);
    if(is~=ns), set(gca,'XTickLabel',''); end
    if(is==1), title('vqf'); end
end
figure, set(gcf,'Color','white');
wgmima= mimaxscG(xs,0,0,0.1);
for is=1:ns,
    subplot(ns+1,1,is), plot(t,xs(:,is));
    ylim(wgmima);
    set(gca,'XTickLabel','');
end
subplot(ns+1,1,ns+1), hold on;
plot(vmtime,vqf(:,2)-vqf(:,1));
plot(vmtime,zeros(size(vmtime)),'r');

end %== function ================================================================
%

