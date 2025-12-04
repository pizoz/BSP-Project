% Fetal ECG detection
% on PhysioNet 2013 Challenge data B
% by Mantas Lukosevicius & Vaidotas Marozas 2013

clear all
close all
clc

% addpath('SMOORT');
% addpath('SMOORT\utils');

filtering = 1;

plotting = 0;	spacing = 5;
saving = 1;		saveExt = '.txt';

% params of data
fs = 1000; % sampling frequency
inSize = 4;
outSize = 1;
runLength = 60000;
introLength = 0; % not used anymore
targetDelay = 0;

dataPreprocess = @loadPreprocData7;
meanMRR = 1.25;
uStd = 1.; % - pre-mqrs removal
uTanh = 0;
dataPostprocess = @detectRrEvents15_grid;% @detectEvents1;% 
rrPow = 0.02;

nEsns = 5;%10;
nUnits = 1000;% 500;% 

disp '======='

%% Create filter
if filtering
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
	b1  = firpm(N, Fo, Ao, W, {dens});
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
    b2  = firpm(N, Fo, Ao, W, {dens});
    %freqz(b2,1,10000,1000)
    
%     b2 = [];
else
	b1 = [];
	b2 = [];
end

%% load the ESN
%load( ['trained_10_esns_' num2str(nUnits) '.mat'] ); %besns
load( ['trained_' num2str(nEsns) '_esns_' num2str(nUnits) '_5.mat'] ); %besns

startStates = cell(nEsns,1);
for iesn = 1:nEsns
	% record starting state (initialized with long zero input)
	besns{iesn} = resetState( besns{iesn} );
	besns{iesn} = runEsn( besns{iesn}, zeros(inSize,500) );
	startStates{iesn} = state( besns{iesn} );
end

%% ========= 'test' ESN on the B datasets
disp --------
for dataFileNr = 0:99;
	tic
    if dataFileNr < 10
        dataFileName = ['b0' num2str(dataFileNr) ];
    else
        dataFileName = ['b' num2str(dataFileNr) ];
    end

	[u, data] = dataPreprocess( ...
		dataFileName, b1, b2, uStd, uTanh, meanMRR );

	ym = zeros( outSize, runLength );
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
	fqrs_det = round(fqrs_det');

	%% save results
	if saving
% 		if dataFileNr < 10
% 			dataFileName = ['b0' num2str(dataFileNr) ];
% 		else
% 			dataFileName = ['b' num2str(dataFileNr) ];
% 		end
		dataFilePath = ['annotations\', dataFileName, saveExt];
		
		save( dataFilePath, 'fqrs_det', '-ASCII' );
	end
	
	%% plot
	if plotting
		set( figure(800+dataFileNr), 'WindowStyle', 'docked' );
		stem( fqrs_det, ones(size(fqrs_det))*4*spacing, '.r');
		hold on
		for lead = 1:inSize
			%plot( data_norm(lead,:) + (inSize-lead)*spacing,'g' );
			plot( data(lead,:) + (inSize-lead)*spacing,'g' );
			plot( u(lead,:) + (inSize-lead)*spacing );
		end
		axis tight
		hold off
	end
	toc
end



%playBeep;
