function [a,b]=notchCoeff(fnotchn,roc)
% --------------------------------------------------------------------------------------------
% Generate notch filter coefficients
% fnotchn = normalized notch frequency
% roc     = module of the pole (roc <1 )
% --------------------------------------------------------------------------------------------
%   Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
%   For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% --------------------------------------------------------------------------------------------

if(nargin<2), roc=0.980; end

ro=1;
phi=2*pi*fnotchn;
% zeri
Zz=[ro*exp(1i*phi); ro*exp(-1i*phi)];
% poli
Zp=[roc*exp(1i*phi); roc*exp(-1i*phi)];

b = poly(Zz);    % MA filter coefficients
a = poly(Zp);    % AR filter coefficients
end %== function ================================================================
%

