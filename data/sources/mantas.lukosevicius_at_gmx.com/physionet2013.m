function [fetal_QRSAnn_est,QT_Interval] = physionet2013(tm,ECG)
% inputs:
%   ECG: 4x60000 (4 channels and 1min of signal at 1000Hz) matrix of
%   abdominal ECG channels.
%   tm : Nx1 vector of time in milliseconds
% output:
%   FQRS: FQRS markers in seconds. Each marker indicates the position of one
%   of the FQRS detected by the algorithm.
%   QT_Interval:   1x1 estimated fetal QT duration (enter NaN or 0 if you do wish to calculate)
% 
% Fetal QRS detection: PhysioNet 2013 Challenge solution
% by Mantas Lukosevicius & Vaidotas Marozas 2013

% ---- check size of ECG ----
if size(ECG,2)<size(ECG,1)
	ECG = ECG';
end

%dataPreprocess = @loadPreprocData7;
meanMRR = 1.25;
uStd = 1.; % - pre-mqrs removal
uTanh = 0;
dataPostprocess = @detectRrEvents15_grid;% @detectEvents1;% 
rrPow = 0.02;

nEsns = 5;%10;
nUnits = 1000;% 500;% 

fs = 1000; % sampling frequency
outSize = 1;
inSize = size(ECG,1);
dataLength = size(ECG,2);


%% load ESN
persistent besns startStates filterB1 filterB2

if isempty( besns )

	load( ['trained_' num2str(nEsns) '_esns_' num2str(nUnits) '_5.mat'] ); %besns

	startStates = cell(nEsns,1);
	for iesn = 1:nEsns
		% record starting state (initialized with long zero input)
		besns{iesn} = resetState( besns{iesn} );
		besns{iesn} = runEsn( besns{iesn}, zeros(inSize,500) );
		startStates{iesn} = state( besns{iesn} );
	end
	
	%% Bandpass filtering before MECG removal
	Fstop1 = 1; %1 %2 %5;               % First Stopband Frequency
	Fpass1 = 8; %8 %8 %12;              % First Passband Frequency
	Fpass2 = 45;              % Second Passband Frequency
	Fstop2 = 49;              % Second Stopband Frequency
	Dstop1 = 0.01;            % First Stopband Attenuation
	Dpass  = 0.057501127785;  % Passband Ripple
	Dstop2 = 0.01;            % Second Stopband Attenuation
	dens   = 20;              % Density Factor
	% Calculate the order from the parameters using FIRPMORD.
	[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(fs/2), [0 1 ...
							  0], [Dstop1 Dpass Dstop2]);
	% Calculate the coefficients using the FIRPM function.
	filterB1  = firpm(N, Fo, Ao, W, {dens});
%     freqz(b1,1,10000,1000)
	
	%% Bandpass filtering after MECG removal
    Fstop1 = 9; %9 %8 % 15;%         % First Stopband Frequency
    Fpass1 = 16; %16 %12% 20;%        % First Passband Frequency
    Fpass2 = 46;% 42;%        % Second Passband Frequency
    Fstop2 = 49;% 49;%        % Second Stopband Frequency
    Dstop1 = 0.01;            % First Stopband Attenuation
    Dpass  = 0.057501127785;  % Passband Ripple
    Dstop2 = 0.01;            % Second Stopband Attenuation
    dens   = 20;              % Density Factor
    % Calculate the order from the parameters using FIRPMORD.
    [N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(fs/2), [0 1 ...
                              0], [Dstop1 Dpass Dstop2]);
    % Calculate the coefficients using the FIRPM function.
    filterB2  = firpm(N, Fo, Ao, W, {dens});

end

	%% preprocess data - from loadPreprocData7.m
	data = ECG';

	u = preprocData7( data, filterB1, filterB2, uStd, uTanh, meanMRR );

	ym = zeros( outSize, dataLength );
	%% test ESNs
	for iesn = 1:nEsns
		%% collect train data for ESN
		%no intro
		besns{iesn} = state( besns{iesn}, startStates{iesn} );
		[besns{iesn} y] = runEsn( besns{iesn}, u ); 
		ym = ym + y;
	end
	ym = ym ./ nEsns;

	%% Post processing
	fqrs_det = dataPostprocess( ym, rrPow );
	fetal_QRSAnn_est = round(fqrs_det');
	
	%% guesstimate QT as half of a mean RR :)
	QT_Interval = 0.5 * mean(fqrs_det(2:end)-fqrs_det(1:end-1));

	
