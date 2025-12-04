function Ser=FecgICAf(Xe,fs,cName,graph,dbFlag,saveFig)
% ----------------------------------------------------------------------------------
%   Fecg: Fetal ecg enhancement by ICA
%   Fixed point algorithm of Hyvarinen with deflationary ortogonalization is applied.
%   In a first attempt the hyperbolic cosine as contrast function is used because it 
%   produces more robust estimates. In case of failure of convergence the algorithm 
%   was run a second time using the kurtosis.
%
% function Ser=FecgICAf(Xe,fs,cName,graph,dbFlag,saveFig)
% Xe      : mixed source signals (one signal per column)
% fs      : sampling frequency
% cName   : record name
% graph   : flag enabling figure drawing
% dbFlag  : flag enabling figure drawing for debugging
% saveFig : flag enabling figure saving
% Ser     : output separated sources (one signal per column)
%
% Author: Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
% For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% -----------------------------------------------------------------------------------

if(nargin<3), cName=''; end
if(nargin<4), graph=0; end
if(nargin<5), dbFlag=0; end
if(nargin<6), saveFig=0; end
% graphD= graph &dbFlag ;

fprintf('\n --------------------------------------------------------- \n');
[progpath, progname] = fileparts(which(mfilename));
fprintf('Program: %s,  record name: %s\n', progname, cName);
%-------------------------------------------------------------
% recording duration
[ndt, ns]=size(Xe);
vtime= [1:ndt]/fs;
%------------------------------------------------------------------------------

for is=1:ns
    Xe(:,is)= (Xe(:,is)-mean(Xe(:,is)))/std(Xe(:,is));
end

if(graph)
    figure, set(gcf,'Color','white');
    for is=1:ns
        subplot(ns,1,is), plot(vtime,Xe(:,is));
        wgmi1= min(Xe(:,is)) -2;
        wgma1= max(Xe(:,is)) +2;
        ylim([wgmi1, wgma1]);
        set(gca,'YTick',[-5 0 5])
        % if(is~=ns), set(gca,'XTickLabel',''); end
        if(is==1), title([cName,': mixed residual signals']); end
    end
end

% ----
rand('state',0);  % to get a reproducible behaviour

% [ica,Ae,cerr]=coshFpDeIca(X,epsilon,maNumIter);
[Se,Ae,cerr]=coshFpDeIca(Xe');  %  X: mixed signals (must be centered)
if(cerr), rand('state',0); [Se,Ae,cerr]=kurtFpDeIca(Xe'); end
% ----
Se=Se';

if(graph || saveFig)
    figure, set(gcf,'Color','white');
    for is=1:ns,
        subplot(ns,1,is), plot(vtime, Se(:,is));
        wgmi1= min(Se(:,is)) -2;
        wgma1= max(Se(:,is)) +2;
        ylim([wgmi1, wgma1]);
        set(gca,'YTick',[-5 0 5])
        % if(is~=ns), set(gca,'XTickLabel',''); end
        if(is==1), title([cName,': residual separated sources']); end
    end
    shg
    if(saveFig), figFmt='png';
        figPath=fullfile('../Figure/',progname);
        if(~exist(figPath,'dir')), mkdir(figPath); end
        figName=fullfile(figPath,[cName,'_ICAf']);
        print(gcf, ['-d',figFmt],figName);
    end
end

Ser=Se;
end %== function ================================================================
%
