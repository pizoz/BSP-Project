function esn = resetState( esn )
%RESETSTATE Resets the ESN's state to zero values.
% by Mantas Lukosevicius

scrLoadEsnFields;
% --------- Set default state
esn.state = zeros( esnSize, 1 );
esn.lastIn = zeros( inSize, 1 );
esn.lastOut = zeros( outSize, 1 );