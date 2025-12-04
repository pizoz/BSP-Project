function argout = state( esn, state )
%state gets/sets state of the ESN
% by Mantas Lukosevicius

if nargin >= 2
	esn.state = state;
	argout = esn;
else
	argout = esn.state;
end