function y = FindEMG(x, Fs)
%
% computes frequency of peak PSD in EMG band
%

[PSD, f] = pwelch(x, Fs*4, 2*Fs, Fs*4, Fs);

[pks, locs] = findpeaks(PSD);
vi     = find(f>49 & f<52);
vx     = find(locs>=vi(1) & locs<=vi(end));
[V, I] = max(pks(vx));

ifP = locs(vx(I));
fP = f(ifP);
M1 = max( PSD(f>(fP-10) & f<(fP-5)) );
M2 = max( PSD(f>(fP+5) & f<(fP+10)) );
M0 = PSD(ifP);

Thr = 10;
if(M0>Thr*M1 || M0>Thr*M2)
    y = fP;
else
    y = 0;
end
