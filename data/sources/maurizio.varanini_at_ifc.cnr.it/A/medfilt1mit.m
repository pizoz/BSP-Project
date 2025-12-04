function xmf = medfilt1mit(x,m,nit)
% --------------------------------------------------------------------------------------------
%  medfilt1mit.m: median filter 
%   median of the first/last m values are assumed to the left and right of v
%   nit - number of iteration of the filter
%   Version 1.00, Date: 07/04/2010
% --------------------------------------------------------------------------------------------
%   Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
%   For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% --------------------------------------------------------------------------------------------

if(nargin<3), nit=1; end
if(nargin<2), m=3; end
if(size(x,2)>size(x,1)), x=x'; colV=0; else colV=1; end
n=length(x);
m2=floor(m/2);
xi=median(x(1:min(n,m)));
xf=median(x(end-min(n,m)+1:end));
xt=[xi+zeros(m2,1); x; xf+zeros(m2,1)];
for it=1:nit
    xx=xt;
    xt=medfilt1(xx,m);
    if(all(~(xt-xx))), break, end
end
if(colV), xmf=xt(m2+1:end-m2); else xmf=xt(m2+1:end-m2)'; end
end %== function ================================================================
%

