function y = normalizeToStd( x, sd )
%normalizeToStd Zero-means and scales rows of x to variance(-s) sd
% (default is 1).
% by Mantas Lukosevicius

if nargin < 2
	sd = 1;
end	
inSize = size( x, 1 );	
if inSize > 1 && size(sd,1) == 1
	sd = ones( inSize, 1 ) * sd;
end
inputStd = std( x, [], 2 );
inputMean = mean( x, 2 );
y = zeros(size(x));
for i = 1:inSize
	y(i,:) = (x(i,:) - inputMean(i)) * ( sd(i) / inputStd(i) );
end
