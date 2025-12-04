function out = movingAverage(in, fs, fn)
% moving averaging of input signal using recursive moving average
%   (c) Jakub Kuzilek, Lenka Lhotska
% 	http://bio.felk.cvut.cz/~kuziljak/     E-mail: jakub.kuzilek@gmail.com
%   Version: 1.0                    Last update:   07/01/2013.
% 					                (Version: 1.0, 07/01/2013)
%
%======================================================
%
% PURPOSE:    This function filters data using recursive moving average
%             filter.
%
% MANDATORY INPUT ARGUMENTS
%   in ..... input data MxN, M - length of data, N - number of leads
%   fs ..... sampling frequency in Hz
%   mode ... 0 - Christov, 1 - Christov + ICA
% OPTIONAL INPUT ARGUMENTS
%   none
% OUTPUT ARGUMENTS
%   out .... filtered signal
%=========================================================
%
% Related Bibliography:
%
% [1] Kuzilek J, Lhotska L. Electrocardiogram beat detection 
%      enhancement using independent component analysis. Medical
%      Engineering and Physics. 2012.
%
%=========================================================


N = round(fs/fn); % length of filter = filter order

gdel = grpdelay([1 zeros(1,N-1) -1],[N -N]);

gdel = ceil(gdel(find(gdel ~= Inf,1,'first')));

in = [in; zeros(gdel,size(in,2))];

out = filter([1 zeros(1,N-1) -1],[N -N],in); % filtering

out = out(gdel+1:end,:);
