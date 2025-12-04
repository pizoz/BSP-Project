function xc=ImpArtElimS(x,thE,wm,pvsc,plotFlag)
% -------------------------------------------------------------------------------------------
%   Impulsive Artifact Canceling from a signal.
%     xc=ImpArtElimS(x,thE,wm,pvsc,plotFlag)
%   Input parameters:
%     x        : array of data to be cleaned from wild points
%     thE      : threshold on absolute error
%     wm       : window length for "deviation" estimation
%     pvsc     : % of value to discard in deviation estimation
%     plotFlag :
%
%   Ouput parameters:
%     xc       : corrected array of data
%
% -------------------------------------------------------------------------------------------
% Author: Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
% For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% -------------------------------------------------------------------------------------------
% This program is free software; you can redistribute it and/or modify it under the terms
% of the GNU General Public License as published by the Free Software Foundation; either
% version 2 of the License, or (at your option) any later version.
%
% This program is distributed "as is" and "as available" in the hope that it will be useful,
% but WITHOUT ANY WARRANTY of any kind; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% --------------------------------------------------------------------------------------------

if (nargin < 1), [x, xv, thE, wm,pvsc,plotFlag]=ImpArtElim_atest; testFlag=1; 
else
    if (nargin < 2 || isempty(thE)), thE=4;  end
    if (nargin < 3 || isempty(wm)),  wm=min(33, length(x)-1);  end  % default values for RR interval series
    if (nargin < 4 ),  pvsc=10;  end
    if (nargin < 5 ),  plotFlag=0;  end
    testFlag=0; 
end
[nr,nc]=size(x);
if(nc>nr), n=nc; x=x'; else n=nr; end

debugFlag=0;

if(wm>n), fprintf('Error: parameter wm must be smaller than length(x)\n'); xc=[]; return; end
% --------
xc=x;

if(debugFlag)
    wLd=60*5000;
    xin=x(1:min(wLd,length(x)));
    xmed= median(xin);    % median value
    xad=abs(xin-xmed);    % vector of absolute error respect to median
    admed=median(xad);    % median of absolute error (deviation)
    admedsd=1.483 *admed;  % robust estimate of standard deviation
    % for a Gaussian signal
    figure; plot(xin);   axis tight;
    figure; plot(xad);   axis tight;
    hold on; plot([1; n],admedsd*[thE; thE],'r');
end

wm=2*floor(wm/2)+1;

xmed = medfilt1mit(xc,wm,1);
xad=abs(xc-xmed);
xadm=maxsc(xad(find(xad>0)),pvsc);

if(plotFlag)
    figure;  hold on;
    plot(x,'r.-'); plot(xmed,'b.-');
    ht=title(['ImpArtElim_atest:',' x (r) - xmed (b)']);
    set(ht,'Interpreter','none');
    shg;   % it need in Linux
    figResize(0, 1, 1, .35); axis tight;
    figure;  hold on;
    plot(xad,'b.-');
    plot(xadm*ones(n,1),'m');
    plot(thE*xadm*ones(n,1),'r');
    ht=title(['ImpArtElim_atest:',' xad (b), xadm (m), th (r)']);
    set(ht,'Interpreter','none');
    shg;   % it need in Linux
    figResize(0, 1, 1, .35); axis tight;
end

kv=find(xad - thE*xadm > 0);
if(~isempty(kv))
    i=1;
    while i<=length(kv)
        iivx=kv(i);
        while (i<length(kv) && kv(i+1)==kv(i)+1), i=i+1; end
        ifvx=kv(i);
        xck=(xc(max(iivx-1,1))+xc(min(ifvx+1,n)))/2;
        xc(iivx:ifvx)= xck;
        i=i+1;
    end
end

if(plotFlag)
    figure;  hold on;
    if(testFlag), plot(xv,'b.-');  end
    plot(x,'r.-'); plot(xc,'g.-');
    ht=title(['ImpArtElim_atest:',' xv (b) - x (r) - xc (g)']);
    set(ht,'Interpreter','none');
    shg;   % it need in Linux
    figResize(0, 1, 1, .35); axis tight;
end

if(nc>nr), xc=xc'; end
if(testFlag),  xc=[];  end
%
end %== function ================================================================

% -------------------------------------------------------------------------------------------------
function [x, xv, thE, wm, pvsc, plotFlag]=ImpArtElim_atest
close all
plotFlag=1;
thE=1.4;  % thE=2;
wm=33;  % wm=33;
pvsc=10;
% --   artificial RR series
% fc_bpm=75bpm => fc=fc_bpm/60= 1.25Hz  => Tc= 0.8s
% fHF=0.25Hz => fHF=fHF/fc
fc=1.25; Tc=1/fc;
fHF=0.25; aHF=0.03*Tc;
fLF=0.1;   aLF= 0.02*Tc;
fHFn=fHF/fc; fLFn=fLF/fc;  % normalized frequencies (cycles/point)
t = (0:Tc:3*60)';
xv = Tc+ aLF*sin(2*pi*fLF*t)+ aHF*sin(2*pi*fHF*t)+ 0.005*Tc*randn(size(t));
fprintf('--- Test for "ImpArtElim" routine ---\n');
fprintf('Artificial  RR interval series - ');
fprintf(' normalized frequencies: LF=%3.2g, HF=%3.2g\n', fLFn, fHFn);
fprintf('n=%d,  thE=%3.2g, wm=%d\n', length(xv), thE, wm);
% --------------------------
% add artifact to the x series
x=xv;
xmed=median(x); mdaerr=median(abs(x-xmed));
fprintf('median= %5.1f,  medAbsErr= %5.2f\n',xmed, mdaerr);
ia=floor(length(x)/4); x(ia)= x(ia)+14*mdaerr; x(ia+2)= x(ia+2)+14*mdaerr;
ia=floor(length(x)/3); x(ia:ia+5)= x(ia:ia+5)+14*mdaerr;
ia=floor(length(x)/2); x(ia:ia+3)= x(ia:ia+3)+8*mdaerr;
ia=floor(length(x)*2/3); x(ia:ia+2)= x(ia:ia+2)+4*mdaerr;
ia=floor(length(x)*3/4); x(ia)= x(ia)-8*mdaerr;
return
end %== function ================================================================
%

