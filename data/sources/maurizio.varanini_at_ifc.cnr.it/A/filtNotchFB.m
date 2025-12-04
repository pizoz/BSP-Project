function xf=filtNotchFB(x,fnotchn,bwn, graf)
% -------------------------------------------------------------------------------------
% notch filter
%
% xf=filtNotchFB(x,fnotchn,bwn,graf)
% x       = input signal
% fnotchn = normalized notch frequency
% bwn     = normalized bandwidth  (Bw/freq)
%
% -------------------------------------------------------------------------------------
%   Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
%   For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% -------------------------------------------------------------------------------------

if(nargin<3), bwn=0.01; end
if(nargin<4), graf=0; end

ro= 1- bwn*2.166;     % approximated formula for two pass filter
[a,b]=notchCoeff(fnotchn, ro);
xf=filtfilt(b,a,x);   % zero-phase forward and reverse filtering

if(graf)
    figure;
    subplot(2,1,1), plot(x);
    subplot(2,1,2), plot(xf);
    title('Notch (FB)');
end
end %== function ================================================================
%

