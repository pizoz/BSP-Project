function [Xi,fsi]=FecgInterp(X,fs,interpFact,cName,graph)
% ----------------------------------------------------------------------------------------------
%   ECG  interpolation
%
% [Xi,fsi]=FecgInterp(X,fs,interpFact,cName,graph)
%
% X          : input signal matrix (one signal per column)
% fs         : sampling frequency
% interpFact : interpolation factor
% cName      : record name
% graph      : flag enabling figure drawing
%
% Xi         : interpolated signal matrix
% fsi        : sampling frequency of interpolated signals
%
% Author: Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
% For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% ----------------------------------------------------------------------------------------------

if(nargin<3), cName=''; end
if(nargin<4), graph=0; end

[progpath, progname] = fileparts(which(mfilename));
fprintf('\n --------------------------------------------------------- \n');
fprintf('Program: %s,  record name: %s\n', progname, cName);

%-------------------------------------------------------------
% recording duration
[ndt, ns]=size(X);

% interpolation by Fourier Transform
Xi=interpft(X,round(interpFact*ndt));
fsi=interpFact*fs;

if(graph)
    figure, set(gcf,'Color','white');
    vtime= [1:ndt]/fs;
    vtimei= [1:size(Xi,1)]/(fsi);
    for is=1:ns,
        subplot(ns,1,is), hold on, plot(vtime,X(:,is),'r'), plot(vtimei,Xi(:,is),'b');
        wgmimaV=mimaxscG(Xi(:,is),0,0,.1);
        ylim(wgmimaV);
        if(is~=ns), set(gca,'XTickLabel',''); end
        if(is==1), title([cName,': original & interpolated ECG']); end
        shg
    end
end
end %== function ================================================================
%
