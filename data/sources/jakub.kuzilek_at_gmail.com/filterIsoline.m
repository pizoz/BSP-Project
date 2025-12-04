function [out] = filterIsoline(in,fs)
% Isoline filtration of ECG using decimation and interpolation
%   (c) Jakub Kuzilek
% 	http://bio.felk.cvut.cz/~kuziljak/     E-mail: jakub.kuzilek@gmail.com
%   Version: 1.0                    Last update:   12/04/2013.
% 					                (Version: 1.0, 12/04/2013)
%
%======================================================
%
% PURPOSE:    This function removes isoline drift from ECG signal.
%
% MANDATORY INPUT ARGUMENTS
%   in ..... input data Mx1, M - length of data (one lead)
%   fs ..... sampling frequency in Hz
% OPTIONAL INPUT ARGUMENTS
%   none
% OUTPUT ARGUMENTS
%   out .... filtered data Mx1
%=========================================================

% in order to remove filter initial slope extend input
fin=[in(1:fs);in];   
np = length(fin);

xd = decimate(fin,round(fs/20),'fir'); % decimation
lbx = medfilt1(xd,10); % median filtering - get isoline
lb = interp(lbx,round(fs/20)); % interpolation
out = fin-lb(1:np); 
out(1:fs)=[];  % result
