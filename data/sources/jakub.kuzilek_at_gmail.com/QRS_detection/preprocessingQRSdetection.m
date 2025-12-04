function out = preprocessingQRSdetection(in, fs, mode)
% preprocessing of ecg signal for QRS detector published in [1]
%   (c) Jakub Kuzilek, Lenka Lhotska
% 	http://bio.felk.cvut.cz/~kuziljak/     E-mail: jakub.kuzilek@gmail.com
%   Version: 1.0                    Last update:   07/01/2013.
% 					                (Version: 1.0, 07/01/2013)
%
%======================================================
%
% PURPOSE:    This function preprocess ECG signal in order to obtain
%             complex lead data.
%
% MANDATORY INPUT ARGUMENTS
%   in ..... input data MxN, M - length of data, N - number of leads
%   fs ..... sampling frequency in Hz
%   mode ... 0 - Christov, 1 - Christov + ICA
% OPTIONAL INPUT ARGUMENTS
%   none
% OUTPUT ARGUMENTS
%   out .... complex lead signal
%=========================================================
%
% Related Bibliography:
%
% [1] Kuzilek J, Lhotska L. Electrocardiogram beat detection 
%      enhancement using independent component analysis. Medical
%      Engineering and Physics. 2012.
%
%=========================================================

in = in-ones(size(in,1),1)*mean(in); % remove mean

out = movingAverage(in, fs, 50); % remove grid noise
out = movingAverage(out, fs, 35); % remove muscle noise

if (mode)
    out = performICA(out);

    out = movingAverage(out, fs, 50); % remove grid noise
    out = movingAverage(out, fs, 35); % remove muscle noise
end

out = complexL(out); % transform to complex lead

out = movingAverage(out, fs, 25); % remove transforming noise
