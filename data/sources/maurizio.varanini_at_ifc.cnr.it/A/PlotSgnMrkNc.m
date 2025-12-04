function PlotSgnMrkNc(sgn,pqrs,freq, name, nband, nsecr, normModef, swYTickLabel)
% --------------------------------------------------------------------------------------------
% PlotSgnMrkNc.m: Plot signals (on multiple horizontal strips) with vertical markers.
%                --> With signal normalization
%                sgn       = signal matrix
%                pqrs      = cell array (or matrix) of time indexes (markers position)
%                freq      = sampling frequency
%                name      = figure title
%                nbands    = number of strips per page
%                nsecr     = number of seconds per row
%                normMode  = normalization mode
%                       (0 -> on all signals,
%                        1 -> for each channel
%                        2 -> on each page
%                        3 -> on each row)
%                swYTickLabel  = 0 -> no yticklabel; 1 -> yticklabel
%
%  Version allowing marker positions in cell array
% --------------------------------------------------------------------------------------------
%   Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
%   For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% --------------------------------------------------------------------------------------------

% convert column to row if necessary
[m,n] = size(sgn);  if n<m,   sgn = sgn.';   end
if(iscell(pqrs))
    [m,n] = size(pqrs); if n>m,  pqrs = pqrs.';   end
    for ia=1:length(pqrs)
        [m,n] = size(pqrs{ia}); if n<m,  pqrs{ia} = pqrs{ia}.';   end
    end
else
    [m,n] = size(pqrs); if n<m,  pqrs = pqrs.';   end
end
nsgn=size(sgn,1);
nsampt=size(sgn,2);
if(nargin<3 || isempty(freq)), freq=1; end
nsect=nsampt/freq;
if(nargin<4 || isempty(name)), name=''; end
if(nargin<5 || isempty(nband)), nband=6; end
if(nargin<6 || isempty(nsecr)), nsecr=nsect/nband; end
if(nargin<7 || isempty(normModef)), normModef=1; end
if(nargin<8), swYTickLabel=1; end

normMode= floor(normModef);
ymargc = normModef-normMode;

ymint=min(min(sgn));
ymaxt=max(max(sgn));
nsampr=nsecr*freq;   % number of samples per row
isampi=1;
isampf=0;
isampf_i=floor(isampf);
vmrkcol=['m', 'b', 'r', 'g', 'k'];

if(normMode==0)
    yminc=ymint*ones(1,nsgn);
    ymaxc=ymaxt*ones(1,nsgn);
    ymaxc= ymaxc + (ymaxc==yminc);
    yrangec=ymaxc-yminc;
    yminc=yminc-ymargc*yrangec;
    ymaxc=ymaxc+ymargc*yrangec;
elseif(normMode==1)
    yminc=min(sgn,[],2);
    ymaxc=max(sgn,[],2);
    ymaxc= ymaxc + (ymaxc==yminc);
    yrangec=ymaxc-yminc;
    yminc=yminc-ymargc*yrangec;
    ymaxc=ymaxc+ymargc*yrangec;
end
while isampf_i < nsampt
    figure, set(gcf,'Color','white');
    if(normMode==2)
        yminc=min(sgn(:,isampi:isampf_i),[],2);
        ymaxc=max(sgn(:,isampi:isampf_i),[],2);
        ymaxc= ymaxc + (ymaxc==yminc);
        yrangec=ymaxc-yminc;
        yminc=yminc-ymargc*yrangec;
        ymaxc=ymaxc+ymargc*yrangec;
    end
    for i=0:nband-1
        subplot(nband,1,i+1);
        isampf = isampi+nsampr-1;
        isampf_i=ceil(isampf);
        if(isampf_i > nsampt), isampf_i=nsampt; end;
        if(normMode>=3)
            yminc=min(sgn(:,isampi:isampf_i),[],2);
            ymaxc=max(sgn(:,isampi:isampf_i),[],2);
            ymaxc= ymaxc + (ymaxc==yminc);
            yrangec=ymaxc-yminc;
            yminc=yminc-ymargc*yrangec;
            ymaxc=ymaxc+ymargc*yrangec;
        end
        ymin=-1; ymax=1;
        yn=(ymax-ymin)./(ymaxc-yminc);
        yoffset=(1*(ymax-ymin)).*((nsgn:-1:1)-1)';
        clear ysgn;
        ysgn=zeros(nsgn,isampf_i-isampi+1);
        for is=1:nsgn
            ysgn(is,:)= yoffset(is) +ymin + yn(is)*(sgn(is,isampi:isampf_i)-yminc(is));
        end
        yming=ymin;
        ymaxg=ymax+max(yoffset);
        xmin=isampi/freq;
        xmax=(isampi+nsampr)/freq;
        plot([isampi:isampf_i]/freq ,ysgn);
        if(~swYTickLabel), set(gca,'YTickLabel',[]); end
        axis tight;
        if(i==0), ht=title(name); set(ht,'Interpreter','none'); end;
        xlim([xmin, xmax]);
        ylim([yming, ymaxg]);
        if(~isempty(pqrs))
            hold on;
            for im=1:size(pqrs,1);  % loop on marker type
                if(iscell(pqrs))
                    vmrk=pqrs{im};
                else
                    vmrk=pqrs(im,:);
                end
                mrkcol=vmrkcol(1+mod(im-1,size(vmrkcol,2)));
                iqrsint=find((isampi < vmrk)&(isampf_i > vmrk));
                if(~isempty(iqrsint))
                    pqrsint=vmrk(iqrsint);
                    fxpqrs=[pqrsint; pqrsint]/freq;
                    fypqrs = ones(size(fxpqrs)); fypqrs(1,:)=yming; fypqrs(2,:)=ymaxg;
                    plot(fxpqrs,fypqrs, [mrkcol ':+']);
                end
            end
        end
        if(isampf_i == nsampt), break; end;
        isampi = isampf_i+1;
    end
end
end %== function ================================================================
%

