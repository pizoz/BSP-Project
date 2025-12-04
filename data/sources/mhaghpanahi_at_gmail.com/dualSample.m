function [xd,dind] = dualSample(x,peaks,phase)

lind = find(peaks,2,'last'); % the indeces of the last two peaks

xx = x(1:lind(1));
I = indFinder(xx,peaks);
k = I(:,1) + floor((phase(1:lind(1))+2*pi*(phase(1:lind(1))<=0))'.*(I(:,2)-I(:,1))./(2*pi));

dind = [1:lind(1)]' + k;
xd = x(dind);
end

function I = indFinder(x,peaks)

% This function finds the indeces of the next two R-peak values
% stored in the "peaks" vector.
% x: input vector 
% peaks: the vector indicating the location of R peaks.
% I: output matrix. i'th row of I contains the indeces of the next two peaks
% from x(i)

n = length(x);
I = zeros(n,2);

for i=1:n
    I(i,:) = find(peaks(i:end),2) - 1;
end

end