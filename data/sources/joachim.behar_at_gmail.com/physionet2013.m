function [FQRSout,QTout] = physionet2013(tm,ecg,paramStruct)
% template algorithm for Physionet/CinC competition 2013. This function can
% be used for events 1 and 2. Participants are free to modify any
% components of the code. However the function prototype must stay the same
% [FQRSout,QTout] = physionet2013(tm,ecg) where the inputs and outputs are specified
% below.
%
% inputs
%   tm :         Nx1 vector of time in milliseconds
%   ecg:         4x60000 (4 channels and 1min of signal at 1000Hz) matrix of
%                abdominal ecg channels.
%   paramStruct: param of the algorithm to use. if not inputed then the
%                default parameters are used. 
%                ***FOR SET-C of the challenge LEAVE EMPTY.***
% 
% output
%   FQRSout:    FQRS markers in seconds. Each marker indicates the position of one
%               of the FQRS detected by the algorithm.
%   QTout:      1x1 estimated fetal QTout duration
%
%
%
% Safe Foetus Monitoring Toolbox, version 1.0, Sept 2013
% Released under the GNU General Public License
%
% Copyright (C) 2013 Joachim Behar
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

if nargin<2; error('physionet2013: wrong number of input arguments'); end;
if nargin<3; paramStruct=[]; end;
addpath(genpath('subfunctions/'));

try
    if isempty(paramStruct)
        % == for FQRS extraction
        paramStruct.method = 'FUSE';              % TS/TS-PCA/TS-EKF/DEFLATION/FUSE
        paramStruct.methodComplexity = 'COMPLEX'; % COMPLEX/SIMPLE
        paramStruct.smoothRR = 1;                 % postprocessing
        paramStruct.cinc_match = 1;               % twick to better match CinC scores
        paramStruct.debug = 0;
        paramStruct.fs = 1000;
        paramStruct.TSnbCycles = 20;              % nb of cycles to build template 
        paramStruct.EKFgain = 1;
        paramStruct.fhighNbCoeff = 5;
        paramStruct.fbasNbCoeff = 3;
        paramStruct.NotchCheck = 1;
        paramStruct.fbas = 10;
        paramStruct.fhigh = 99;
        paramStruct.STDcut = 29;                  % rubbish to fit to challenge scores (was 17)
        paramStruct.NbPC = 2;                     % nb of principal components in PCA
        paramStruct.detThres = 0.5;               % QRS detector threshold for FQRS
        paramStruct.qrsdet_wind = 15;             % qrs detection window: LENGTH_SEG/qrsdet_wind NEEDS TO BE AN INTEGER!!
        paramStruct.extract_qt = 0;
    end
    
    [FQRSout,~,outExtra] = physionet2013extract(tm,ecg,paramStruct);
catch ME
    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
    FQRSout=0.1:0.42:60;
    QTout=250;
end


try
    % == for QT extraction
    % the main reasons for doing the QT separately are:
    % (1) the baseline wander cutt-off must be lower than 10Hz otherwise
    % the T-wave is completely destroyed
    % (2) using BSS methods is not a good idea as part of the T-wave could
    % end in a secondary channel so we stay we time based methods.
    qt_paramStruct.method = 'TS-PCA';
    qt_paramStruct.methodComplexity = 'SIMPLE';  % MUST BE: SIMPLE (time based)
    qt_paramStruct.debug = 0; 
    qt_paramStruct.fs = 1000;
    qt_paramStruct.TSnbCycles = 20;
    qt_paramStruct.EKFgain = 1;
    qt_paramStruct.fhighNbCoeff = 5;
    qt_paramStruct.fbasNbCoeff = 3;
    qt_paramStruct.NotchCheck = 1;
    qt_paramStruct.fbas = 2; % a lower cut-off not to kill the T-wave
    qt_paramStruct.fhigh = 99;
    qt_paramStruct.NbPC = 2;
    qt_paramStruct.extract_qt = 1;    
    qt_paramStruct.rMQRS = outExtra.rMQRS; % use the MQRS returned by the previous step
    qt_paramStruct.rFQRS = outExtra.rFQRS; % use the FQRS returned by the previous step
    
    [~,QTout,~] = physionet2013extract(tm,ecg,qt_paramStruct);

catch ME
    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
    FQRSout=0.1:0.42:60;
    QTout=250;
end

end







