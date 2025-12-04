function residual = MECGcancellation(peaks,ECG,fs) %,nbCycles)
% MECG cancellation algorithm inspired from [1].
%
% inputs:
%   fs: sampling frequency
%   nbCycles: number of cycles on which to build the mean MECG template
%   ECG: matrix of abdominal ECG channels.
%   peaks: MQRS markers in seconds. Each marker corresponds to the
%   position of a MQRS.
%
% output:
%   residual: residual containing the FECG.
%
% Author: Joachim Behar - IPMG Oxford (joachim.behar@eng.ox.ac.uk)
% Last updated: 03_02_2013
%
% [1] Martens S et al. A robust fetal ECG detection method for
% abdominal recordings. Physiol. Meas. (2007) 28(4) 373�388
% 
% Last updated 2013.04.27 by Mantas Lukosevicius

% ---- constants ----
%if nargin < 4
%	r = size(peaks,2);
% else
% 	r                   = nbCycles;
% end
Pstart              = 0.25*fs-1;
Tstop               = 0.45*fs;
N                   = length(peaks);       % number of MECG QRS
ecgLength			= length(ECG);
ECG_temp            = zeros(1,ecgLength);

% ---- ECG template ----
templateStart = 1;
while templateStart < N && peaks(templateStart)-Pstart < 1
	templateStart = templateStart+1;
end
templateEnd = N;
while templateEnd >= templateStart && peaks(templateEnd)+Tstop > ecgLength
	templateEnd = templateEnd-1;
end
	
ECG_last_r_cycles   = zeros(0.7*fs,templateEnd-templateStart+1);

for i=templateStart:templateEnd
    peak_nb = peaks(i); 
    ECG_last_r_cycles(:,i) = ECG(peak_nb-Pstart:peak_nb+Tstop)';
end
ECG_mean = mean(ECG_last_r_cycles,2);

% ---- MECG cancellation ----
% ---- MECG cancellation ----
M  = zeros (0.7*fs,3);
M(1:0.2*fs,1)           = ECG_mean(1:Pstart-0.05*fs+1);
M(0.2*fs+1:0.3*fs,2)    = ECG_mean(Pstart-0.05*fs+2:Pstart+0.05*fs+1);
M(0.3*fs+1:end,3)       = ECG_mean(Pstart+2+0.05*fs:Pstart+1+Tstop);
for i=1:N
	%TODO: handle end cases
	if peaks(i)>Pstart
		Pstart_ = Pstart;
		if ecgLength-peaks(i) >= Tstop
			Tstop_ = Tstop;
			M_ = M;
		else
			Tstop_ = ecgLength-peaks(i);
			M_ = M(1:0.7*fs-(Tstop-Tstop_),:);
			if all(M_(:,3)==0)
				M_ = M_(:,1:2);
			end
		end
	else
		Pstart_ = peaks(i)-1;
		Tstop_ = Tstop;
		M_ = M(Pstart-Pstart_+1:end,:);
		if all(M_(:,1)==0)
			M_ = M_(:,2:3);
		end
	end
	
	
	%M_ = M(Pstart-Pstart_+1:0.7*fs-(Tstop-Tstop_),:);
	
	a = (M_'*M_)\M_'*ECG(peaks(i)-Pstart_:peaks(i)+Tstop_)';
% 	if ~any(isnan(a))
% 		ECG_temp(peaks(i)-Pstart_:peaks(i)+Tstop_) = ...
% 		a(1)*M_(:,1)'+a(2)*M_(:,2)'+a(3)*M_(:,3)';		
% 	end

	if length(a) == 3
		ECG_temp(peaks(i)-Pstart_:peaks(i)+Tstop_) = a(1)*M_(:,1)'+a(2)*M_(:,2)'+a(3)*M_(:,3)';		
	else
		ECG_temp(peaks(i)-Pstart_:peaks(i)+Tstop_) = a(1)*M_(:,1)'+a(2)*M_(:,2)';	
	end


end
% compute residual
residual = ECG - ECG_temp;
end
