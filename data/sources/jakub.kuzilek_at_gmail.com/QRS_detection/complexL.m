function out = complexL(in)
% transforms ecg into a complex lead using equation:
%        Y[i] = 1/L * sum(abs(Xj[i+1]-Xj[i-1]))
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
% OPTIONAL INPUT ARGUMENTS
%   none
% OUTPUT ARGUMENTS
%   out .... complex lead signal, two zeros on start compensates time  delay
%=========================================================
%
% Related Bibliography:
%
% [1] Kuzilek J, Lhotska L. Electrocardiogram beat detection 
%      enhancement using independent component analysis. Medical
%      Engineering and Physics. 2012.
%
%=========================================================


out = [0; 0; sum(abs((in(3:end,:)-in(1:end-2,:))'),1)'/size(in,2)]; % sum + /n is better than mean, don't know why
