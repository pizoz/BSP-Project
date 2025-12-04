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

%% ========= 'test' ESN on the B datasets
disp --------
for dataFileNr = 0:99;
	tic
    if dataFileNr < 10
        dataFileName = ['b0' num2str(dataFileNr) ];
    else
        dataFileName = ['b' num2str(dataFileNr) ];
    end

	%% Open data files
	data = rdsamp( ['duomenys\' dataFileName] );
	tm = data(:,1);
	data = data(:,2:1+inSize);
	
	%% do fQRS detection
	%size(data)
	[fqrs_det qt] = physionet2013(tm, data);

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
%			plot( u(lead,:) + (inSize-lead)*spacing );
		end
		axis tight
		hold off
	end
	toc
end



%playBeep;
