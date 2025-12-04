function [FQRSout,m1,m2] = guess_fqrs(FQRS,rMQRS,fs)
% this function is used to guess what the FHR can be or completely guess 
% it when the extraction completely failed. It then generate a corresponding 
% sensible RR time series.
% 
% inputs: 
%   FQRS: FQRS time series in ms
%   rMQRS: referene mother QRS
%   fs: sampling frequency
%
% output:
%   m: representative HR mode (in bpm)
%
% Author: Joachim Behar - IPMG Oxford (joachim.behar@eng.ox.ac.uk)
% Last updated: 01-09-2013


% == cut into first 30 and second 30 sec 
FQRSp1 = FQRS(FQRS<30*fs);
FQRSp2 = FQRS(FQRS>=30*fs);

m1 = select_mode(FQRSp1,rMQRS,fs);
m2 = select_mode(FQRSp2,rMQRS,fs);

rMQRSp1 = rMQRS(rMQRS<30*fs);
rMQRSp2 = rMQRS(rMQRS>=30*fs);

rm1 = mode(60./diff(rMQRSp1./fs));
rm2 = mode(60./diff(rMQRSp2./fs));

if abs(rm1-m1)>15 && m1~=0
    FQRSpart1 = 0.1:60/m1:30;
else
    FQRSpart1 = 0.1:0.42:60;
    m1 = 142.8571;
end

if abs(rm2-m2)>15 && m2~=0
    FQRSpart2 = max(FQRSpart1):60/m2:60;
else
    if abs(rm1-m1)>15 && m1~=0
        % if the first mode is ok then keep it for the other half of the
        % seg
        FQRSpart2 = max(FQRSpart1):60/m1:60;
        m2=m1;
    else
        FQRSpart2 = max(FQRSpart1):0.42:60;
        m2 = 142.8571;
    end
end

FQRSout = [FQRSpart1 FQRSpart2(2:end)];

end


function mout = select_mode(FQRS,rMQRS,fs)
% select mode corresponding to FQRS time series. The assumption behind this
% function is that when a certain number or estimated HR following each
% other are the same +/- some physiological variation then it is likely
% that this HR is representative of the FHR on the corresponding segment.
% This is better than taking the mode (most representative HR 
% probabilistically speaking) because it takes into consideration the fact
% that following HR datapoints should not differ substrantially from each
% other. 
% NOTE: this is not to be used for a practical application!
% 
% inputs: 
%   FQRS: FQRS time series in ms
%   rMQRS: MQRS reference- to make sure we are not catching that
%   fs: sampling frequency
%
% output:
%   m: representative HR mode (in bpm)
%
% Author: Joachim Behar - IPMG Oxford (joachim.behar@eng.ox.ac.uk)
% Last updated: 02-09-2013

mHR  = 60./diff(rMQRS./fs); medmHR = median(mHR);
HR  = 60./diff(FQRS./fs);
HRV = diff(HR);
N = length(HRV);
MASK_SIZE = 5;
TOL = 10; % in bpm
m = [];
compt = 0;

for ii=1:N-MASK_SIZE
    if HRV(ii:ii+MASK_SIZE-1)<TOL
        medfHR = median(HR(ii:ii+MASK_SIZE-2));
        if abs(medfHR-medmHR)>10
            compt=compt+1;
            m(compt) = medfHR;
            fprintf('Extracting mode at %f \n',m(compt));
        end
    end
end

if isempty(m)
    mout = 0;
else
    mout = median(m);
end

end





















