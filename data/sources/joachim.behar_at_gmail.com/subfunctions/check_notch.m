function NeedFiltering = check_notch(ecg,fc,fs,debug)
% this function checks if there is an energy peak in the frequency band
% +/-1 Hz around fc. The reason for this function is to avoid doing too much
% filtering if we can, i.e. decide to Notch filter or not the signal, 
% because this might create ripples that affect the ecg.
% 
% inputs
%   ecg:    a single channel ecg
%   fc:     the central frequency we are studying
%   fs:     sampling frequency
% 
% ouputs
%   NeedFiltering: boolean. 1 if there is an energy peak and 0 otherwise
%
%
% FECG extraction toolbox, version 1.0, Sept 2013
% Released under the GNU General Public License
%
% Copyright (C) 2013  Joachim Behar
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2013
% joachim.behar@eng.ox.ac.uk
%
% Last updated : 02-08-2013
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

% == managing inputs
if nargin<2; error('check_notch: wrong number of input arguments \n'); end;
if nargin==2; fs=1000; debug=0; end;
if nargin==3; debug=0; end;
if size(ecg,2)>size(ecg,1); ecg=ecg'; end;

% == general parameters
NFFT = 10000;        % number of NFFT coefficients
pace = fs/NFFT;      % one data point corresponds to pace Hz
df = 4;              % delta frequency: we look at +/-df around fc
NeedFiltering = 0;

% == core function
try
    % == define spectral band we are looking for a peak into
    specBand = round(fc/pace)-round(df/pace):round(fc/pace)+round(df/pace);

    % == PSD
    [Pxx F1] = pwelch(ecg(:,1),5*fs,fs,NFFT,fs,'onesided');
    PxxdB = 10*log10(Pxx/(mean(Pxx)));

    % == look for local max in specBand
    [~,pos] = max((PxxdB(specBand)));

    % == and the frequency it corresponds to
    CorrespondingFreq = F1(round(fc/pace)-round(df/pace)+pos-1);

    % == is it in +/-1Hz around fc?
    if CorrespondingFreq<fc+1 && CorrespondingFreq>fc-1; NeedFiltering = 1; end;
catch ME
    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
    NeedFiltering = 1;
end

% == debug
if debug
    plot(F1,PxxdB); hold on, plot(F1(specBand),1,'+r');
    hold on, plot(CorrespondingFreq,PxxdB(round(fc/pace)-round(df/pace)+pos-1),'+g');
    set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold');
end

end

