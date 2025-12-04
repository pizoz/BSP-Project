function ma = maxsc(v,perc)
% ---------------------------------------------------------------------------------------------
% maxsc.m: Compute the maximum value of a vector excluding the distribution tail.
%   ma = maxsc(v,perc)
%            "perc" = % of max values to be excluded
%
% ---------------------------------------------------------------------------------------------
%   Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
%   For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% ---------------------------------------------------------------------------------------------

vo=sort(v(:));
if(nargin<2), perc=5; end

if(perc<0), fi=length(v)+perc;
else fi=length(v)-floor(length(v)*perc/100); end

ma = vo(fi);
end %== function ================================================================
%

