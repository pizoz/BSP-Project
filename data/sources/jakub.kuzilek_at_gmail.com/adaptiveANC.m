function [signal50] = adaptiveANC(signal,fs,f);
%[signal50]=adaptiveANC(signal,fs)
%filter the inteferance from the power line
%
%signal-to be filtered
%fs-sample frequency
%signal50-signal with eliminated power line interference
%
%EXAMPLE: signal50=adaptiveANC(signal,250)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signal_upr=[signal(1:fs);signal];   %pomoc pri odstraneni nabehu filtru
N = length(signal_upr);
% f = 33;                         % frekvence, ktera bude odstranena
r1 = sin(2*pi*f/fs*(1:N));
r2 = cos(2*pi*f/fs*(1:N));
% mu = 5.8/sum(signal.*signal);   % konvergencni konstanta - idealni hodnota je kolem 0.04
mu = 0.04;
u = zeros(1,N);
e = zeros(1,N);
W1 = zeros(1,N); W1(1) = 0;
W2 = zeros(1,N); W2(1) = 0;
for n = 1:N-1
    %filtrace
    u(n) = W1(n)*r1(n) + W2(n)*r2(n);
    e(n) = signal_upr(n) - u(n);    % signal bez sumu
    %aktualizace
        W1(n+1) = W1(n) + mu*e(n)*r1(n); 
        W2(n+1) = W2(n) + mu*e(n)*r2(n);
end
e(1:fs)=[]; %odstraneni pridaneho signalu
signal50=e;