function [out,out1] = findMeanHigh(x,limit)
% FHR outliers detection
%   (c) Lukas Hnizdo, Jakub Kuzilek
% 	http://bio.felk.cvut.cz/~kuziljak/     E-mail: jakub.kuzilek@gmail.com
%   Version: 1.0                    Last update:   24/08/2013.
% 					                (Version: 1.0, 23/08/2013)
%
%======================================================
%
% PURPOSE:    This function detects outliers in RR data. Developed for
%             purpose of Physionet Challenge 2013.
%
% MANDATORY INPUT ARGUMENTS
%   RR ...... fetal HRV Mx1, M - length of data (one lead)
%   limit ... detection limit in %
% OPTIONAL INPUT ARGUMENTS
%   none
% OUTPUT ARGUMENTS
%   out .... positions of outliers
%=========================================================
lim = 1+limit/100;
N = length(x);
out = zeros(N,1);

for n = 2:N-1
     %OZNACENI HIGH POINTS
     mean = (x(n-1)+x(n+1))/2;
     if((x(n) > mean*lim))
          out(n) = 1;       
     end
end

out1 = find(out);

n = [];
t = [];
time = cumsum(x);
for i=1:length(out1)
    n = [n; x(out1(i))];
    t = [t; time(out1(i))];
end

% figure
% plot(time,x);
% hold on
% stem(t,n,'r');
% xlabel('time [ms] -->');
% ylabel('time [ms] -->');
% legend('original RR', 'cleaned RR');











