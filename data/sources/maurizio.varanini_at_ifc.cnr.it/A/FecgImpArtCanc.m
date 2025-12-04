function Xc=FecgImpArtCanc(X,fs,cName,graph,dbFlag)
% -------------------------------------------------------------------------------------
% Impulsive artifact removal from each ECG channel,
% applying the function "ImpArtElimS" to each channel.
%
% Xc=FecgImpArtCanc(X,fs,cName,graph,dbFlag)
%
% X       : input signal matrix (one signal per column)
% fs      : sampling frequency
% cName   : record name
% graph   : flag enabling figure drawing
% dbFlag  : flag enabling figure drawing for debugging
%
% Xc      : output cleaned signal matrix (one signal per column)
%
%   Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
%   For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% -------------------------------------------------------------------------------------

if(nargin<3), cName=''; end
if(nargin<4), graph=0; end
if(nargin<5), dbFlag=0; end

[ndt, ns]=size(X);
%-------------------------------------------------------------------------------
% Set the first ten samples to the median of the following three
npti=10; nptim=3;
for is=1:ns
    X(1:npti,is)=median(X(npti+1:npti+nptim,is));
end
%-------------------------------------------------------------------------------
thE=4;               %  threshold on absolute error
wm=fix(0.06*fs);   %  window length for "deviation" estimation (60ms)
pvsc=2;
Xc=zeros(size(X));
for is=1:ns
    Xc(:,is)=ImpArtElimS(X(:,is),thE,wm,pvsc,dbFlag);
end

if(graph)
    vtime= [1:ndt]/fs;
    figure, set(gcf,'Color','white');
    for is=1:ns,
        subplot(ns,1,is), hold on, plot(vtime,X(:,is),'r'), plot(vtime,Xc(:,is),'b');
        wgmimaV=mimaxscG(Xc(:,is),0,0,.1);
        ylim(wgmimaV);
        if(is~=ns), set(gca,'XTickLabel',''); end
        if(is==1), title([cName,': original (r) & cleaned (b) ECG']); end
        shg
    end
end
end %== function ================================================================
%
