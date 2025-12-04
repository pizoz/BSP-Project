function [u, data, data_orig, sumlead, MQRS] = preprocData7( ...
	data, filterB1, filterB2, inStd, inTanh, meanMRR )
%loadPreprocData6 + normalization after MECG removal
% Please save modifications of this file under name ...2, etc.
% by Mantas Lukosevicius & Vaidotas Marozas 2013]
	
% default parameters
	if nargin < 6
		meanMRR = 1.5;
		if nargin < 5
			inTanh = 1;
			if nargin < 4
				inStd = 1;
			end
		end
	end
		
	%data = data'; % comes as 4x60000
	
	fs = 1000;
	inSize = size(data,2);
	dataLength = size(data,1);

	%% === replace underflow outliers:
	outlier_mask = ( data == -32768 );
	if any(any(outlier_mask))
		% find should-be value
		outlier_correction = min(min(data(~outlier_mask))) - 1;
		% replace outliers
		data( outlier_mask ) = outlier_correction;
	end
	%disp([int2str(sum(sum(outlier_mask))),' outliers corrected']);

	% convert data to uV
	data = data / 1000;
	data_orig = data;
	
	for lead = 1:inSize
		data(:,lead) = data(:,lead) - mean(data(:,lead));
		
% 		if abs( data(1,lead) - data(5,lead) ) > 2
			% flatten the begining to eliminate the initial transient when
			data(1:5,lead) = data(5,lead);
% 		end
	end
	
	if nargin > 1 && ~isempty(filterB1)
		%% Bandpass filtering   
		data = filtfilt(filterB1,1,data);  

% 		% mirror the signal to eliminate the initial transient when filtering:
% 		data2 = [flipud(data); data];
% 		data2 = filtfilt(filterB,1,data2);
% 		data = data2(dataLength+1:end,:);
	end

	data = data';
	
	%% remove MECG
	data = normalizeToStd( data, inStd );
	% flip the channels

    maxMedian = 0;
% 	sumlead = zeros(1,dataLength);
	for lead = 1:inSize 
		med = skewness(data(lead,:)); %  median(data(lead,:));
		if med < 0 
			data(lead,:) = -data(lead,:);
		end
%		sumlead = sumlead + data(lead,:) * (-abs(med));
		if abs(med) > maxMedian
			sumlead = data(lead,:);
			maxMedian = abs(med);
		end
	end
	%sumlead = mean(data);

%  	sumlead = data(1,:);
% 	for lead = 2:inSize 
% 		sumlead1 = sumlead + data(lead,:);
% 		sumlead2 = sumlead - data(lead,:);
% 		if sum(abs(sumlead1)) > sum(abs(sumlead2))  
% 			sumlead = sumlead1;
% 		else
% 			sumlead = sumlead2;
% 		end
% 	end
	
	%detect Maternal R peaks
	MQRS = PeakDetection( sumlead, meanMRR/fs, 1 );
	
	%	[temp, MQRS] = fun_RRI_estimation( sumlead, fs ); 

	% ---- MECG cancellation ----
	u = zeros(size(data));
	for lead = 1:inSize             % run algorithm for each channel				
		u(lead,:) = MECGcancellation2( MQRS, data(lead,:), fs );
% 		med = skewness(u(lead,:)); %  median(data(lead,:));
% 		if med < 0 
% 			u(lead,:) = -u(lead,:);
% 		end
	end

	if ~isempty(filterB2)
		%% Bandpass filtering   
		u = filtfilt(filterB2,1,u')';  

% 		% mirror the signal to eliminate the initial transient when filtering:
% 		data2 = [flipud(data); data];
% 		data2 = filtfilt(filterB,1,data2);
% 		data = data2(dataLength+1:end,:);
	end

	u = normalizeToStd( u, inStd );

	%% format data for ESN

	outSize = 1;
	% normalize data to std
	%u = normalizeToStd( u, inStd );
	if inTanh
		u = tanh(u);	
	end
	dataLength = size(u,2);

end