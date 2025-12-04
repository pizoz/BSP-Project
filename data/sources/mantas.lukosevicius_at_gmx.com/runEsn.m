function [ esn, output, states ] = run( esn, input )
%RUN Runs ESN with a given input.
% by Mantas Lukosevicius

scrLoadEsnFields;

if inSize ~= size( input, 1 )
    error('Wrong size of the input. Must be [inputSize,runLength].');
end
runLength = size( input, 2 );

if runLength == 0 
	return
end

noiseLevel2 = noiseLevel * 2;
retainRate = 1 - leakingRate;
output = zeros( outSize, runLength );
x = esn.state;
% lastOut - later

% split weight matrices for better performance
Woutb = Wout(:,1);
Woutu = Wout(:,2:inSize+1);
Woutx = Wout(:,inSize+2:end);
Winb = Win(:,1);
Winu = Win(:,2:end);

if nargout >=3
	% leaky with noise, collect states
	states = zeros( esnSize, runLength );
	for i = 1:runLength 
		u = input( :, i );
		newx = f( act, Winb + Winu * u + W * x ); 
		if noiseLevel2 ~= 0
			newx = newx + ( rand( esnSize, 1 ) - 0.5 ) .* noiseLevel2;
		end
		if retainRate ~= 0
			x = x .* retainRate + newx .* leakingRate;
		else
			x = newx;
		end
		y = f( actOut, Woutb + Woutu * u + Woutx * x ); %Wout * [u; x] );
		output( :, i ) = y;
		states( :, i ) = x;
	end
else
	% leaky with noise, no collect states (separate for speed)
	for i = 1:runLength 
		u = input( :, i );
		newx = f( act, Winb + Winu * u + W * x ); % Win * u + W * x ); 
		if noiseLevel2 ~= 0
			newx = newx + ( rand( esnSize, 1 ) - 0.5 ) .* noiseLevel2;
		end
		if retainRate ~= 0
			x = x .* retainRate + newx .* leakingRate;
		else
			x = newx;
		end
		y = f( actOut, Woutb + Woutu * u + Woutx * x ); %Wout * [u; x] );
		output( :, i ) = y;
	end
end
% save state
esn.state = x;
esn.lastIn = u;
esn.lastOut = y;
