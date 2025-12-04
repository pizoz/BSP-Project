function out = performICA(in)
% ICA computation and component selection function
%   (c) Jakub Kuzilek, Lenka Lhotska
% 	http://bio.felk.cvut.cz/~kuziljak/     E-mail: jakub.kuzilek@gmail.com
%   Version: 1.0                    Last update:   07/01/2013.
% 					                (Version: 1.0, 07/01/2013)
%
%======================================================
%
% PURPOSE:    This function deploys JADE ICA algorithm and selects 
%             components containing ECG.
%
% MANDATORY INPUT ARGUMENTS
%   in ..... input data MxN, M - length of data, N - number of leads
% OPTIONAL INPUT ARGUMENTS
%   none
% OUTPUT ARGUMENTS
%   out .... selected components MxO, where O<=N
%=========================================================
%
% Related Bibliography:
%
% [1] Kuzilek J, Lhotska L. Electrocardiogram beat detection 
%      enhancement using independent component analysis. Medical
%      Engineering and Physics. 2012.
%
%=========================================================

M = size(in,2); % input ecg data size

% if there is small number of mixtures perform trick 
% to double them (widely used in ICA) - move signals
% by one sample and append them to original data
if M < 4 
    in = [in, [in(2:end,:);zeros(1,M)]];
end

% perform JADE ICA algorithm
in = real((jadeR(in')*in')');

% compute kurtosis of components
x0 = in - ones(size(in,1),1)*mean(in);
s2 = mean(x0.^2); % this is the biased variance estimator
m4 = mean(x0.^4);
kurt = m4 ./ s2.^2;
% kurt = kurtosis(in);

% sort and normalize them
[skurt, order] = sort(kurt, 'descend');

skurt = cumsum(skurt);

skurt = skurt/skurt(end);

% choose componets whose sum of kurt is larger than 75 percents of sum(kurt)
order = order(1:find(skurt >= 0.76));

% output data
out = in(:,order);

