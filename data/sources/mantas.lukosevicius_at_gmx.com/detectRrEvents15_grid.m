function [events, val ] = detectRrEvents15_grid( y, rrPow )
%detectRrEvents14 detects fetal ECG R events form the indicator signal
% detectRrEvents14 + thick tail rrPdf from detectRrEvents0
% for PhysioNet 2013 challenge
% by Mantas Lukosevicius 2013

% if nargin < 2
%     rrPow = 1/100;
% end
% 
	% RR interval duration is assumed to have normal disribution 
	% from cinc_total_rr_hist.m: RR mean = 442.1759, std = 45.99
% 	rrMean = 442;
	rrStd = 46;
	% from cinc_total_75_rr_hist.m: RR mean = 424.8335, std = 83.6907 
	% from cinc_total_rr_hist_ign.m: RR mean = 423.8379, std = 39.867
	%  (has outliers up to ~5000 in A54)
	rrMean = 424;
%	rrStd = 40;
% 	rrStd = 83.7;
	
	% change of RR interval duration is assumed to have normal disribution 
	% from cinc_total_drr_hist.m: dRR mean = -0.06497, std = 18.8144
	%drrMean = 0;
	drrStd = 19;

	runLength = size( y, 2 );
	%yMean = mean(y);
		
	rrMin = 300;
	rrMax = 500; %550; % 2*rrMean; % ~0.9s


	drrMax = drrStd * 3;
	drrPdf = normpdf( -drrMax:drrMax, 0, drrStd )';
	drrPdf = drrPdf / max(drrPdf);

	% to regain the rythm once it's lost
% 	broadFactor = 0.17;
% 	broadPdf = ones(1, rrMax );
	
	broadFactor = 0.35;
	broadPdf = [2:2:(2*rrMean) (2*rrMean):-2:2]/(2*rrMean);

	% to keep the rithm
	rrPdf = normpdf( 1:(2*rrMean), rrMean, rrStd );
	rrPdf = rrPdf / max( rrPdf );
	
	rrPdf = ( rrPdf + broadPdf .* broadFactor) / (1+broadFactor);

	rrPdf = rrPdf(rrMin:rrMax)' .^ rrPow;% ( 1/100 );
		
	
	%% convolve with a peak of target signal
	
	% indicator width 5:
	peakWidth = 5; %should be odd and the same as for target
	peakCenter = (peakWidth-1)/2 + 1;
	peak = zeros(peakWidth,1);
	peak(1:peakCenter) = 1/peakCenter:1/peakCenter:1;
	peak(peakCenter:end) = 1:-1/peakCenter:1/peakCenter;
	peak = peak ./ sum(peak);
	
	yc = conv(peak,y);
	y = yc(peakCenter:end-peakCenter+1);

	%% find a probable sequence of RR events
 	y(y<0) = 0; % for multiplicative scores
	
	b = 1/6;% 1/4;% 1/8;% 1/10;% 1/5;% 1/3;% 1/2;%
	
	valOffset = rrMin-drrMax-1;
	
	val = zeros( rrMax+drrMax-valOffset, runLength ); % (lastRR, currentR)
	from = zeros( rrMax+drrMax-valOffset, runLength );
	valWind = zeros(rrMax-rrMin+1,1+2*drrMax);
	drrPdfM = (drrPdf * ones(1,rrMax-rrMin+1))';
	for i = 1:runLength
		if i < rrMin+1
			%valWind = valWind
			val(:,i) = y(i);
		else

			valWind(2:rrMax-rrMin+1,1:2*drrMax) = valWind(1:rrMax-rrMin,2:1+2*drrMax);
			valWind(1,:) = val(rrMin-valOffset-drrMax:rrMin-valOffset+drrMax,i-rrMin);
			
			if i < rrMax + 1
				valWind(i-rrMin:end,:) = y(i);
				rrRange = rrMin:i-1;
			else
				rrRange = rrMin:rrMax;
			end			
			
			%rrRange = rrMin:min(rrMax,i-1);
			rrIndices = -valOffset+drrMax + rrRange + size(val,1)*(i-1-rrRange);			
			valWind(1:length(rrRange),end) = val( rrIndices );
			
			[vs, fs] = max( valWind .* drrPdfM, [], 2 );
			val(rrMin-valOffset:rrMax-valOffset,i) = rrPdf .* y(i).^b .* vs.^(1-b);
			from(rrMin-valOffset:rrMax-valOffset,i) = (rrMin:rrMax)' + fs - (drrMax+1);
			
		end
	end
			
	% trace back
%  	[ w, j ] = max( val(:,runLength-rrMean-2*rrStd:runLength) );
% 	[ v, i ] = max( w );
 	[ w, i ] = max( val(:,runLength-rrMean-2*rrStd:runLength), [], 2 );
	% bias the last beat according to rrPdf
	%rrBias = b; % 1/6;
	%[ v, j ] = max( w' .* rrPdf( valOffset+1:end ).^rrBias );
	[ v, j ] = max( w );
	i = i(j) + runLength-rrMean-2*rrStd-1;
	j = j+valOffset;
	events = [];
	while i > 0 && j > 0 
		events = [i, events];
		i2 = i-j;
		j = from(j-valOffset,i);
		i = i2;
	end

end