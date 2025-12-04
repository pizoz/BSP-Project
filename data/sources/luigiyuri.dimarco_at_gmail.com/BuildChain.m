function Tq = BuildChain(S, Fs, mu, sigma)
%
% function: attempts to build sequence of annotatios given 
%           normal distribution pdf (mu, sigma)
%


[pks, locs] = findpeaks(S);

thr  = prctile(pks, 50);
locs = locs(pks>thr);

pks  = pks(pks>thr);

tx  = locs;
cx  = pks';
Tq  = zeros(1, numel(tx));

vi     = find(tx <= max(600, locs(1))); % search 1st element within 600ms window (arbitrary)

clear pks locs

DT = 2*Fs;
W  = normpdf(1:DT, mu, sigma);  % 2s window (normal pdf)

[~, I] = max(cx(vi));
Tq(1)  = tx(vi(I));

for q=2:numel(tx)
    t0   = Tq(q-1);
    vi   = find(tx>t0 & tx<(t0+DT));
    if(isempty(vi))
        break
    end
    pksW = W(tx(vi)-t0+1).*cx(vi);
    [~, I] = max(pksW);
    Tq(q)  = tx(vi(I));
end
Tq(q:end) = [];