% Helper script that loads ESN fields into local variables
% by Mantas Lukosevicius

W = esn.W;
Win = esn.Win;
Wout = esn.Wout;
inSize = size( Win, 2 ) - 1;
esnSize = size( W, 1 );
outSize = size( Wout, 1 );

actOut = esn.outActivation;
act = esn.activation;

leakingRate = esn.leakingRate;
noiseLevel = esn.noiseLevel;