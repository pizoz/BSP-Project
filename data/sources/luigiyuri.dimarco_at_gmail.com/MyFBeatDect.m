function [fqrs, ch_ref, X, Xcl] = MyFBeatDect(ECG, Fs, bs, mqrs, bPlot)
%
% function: computes beat detection
% 
% IN: 
%      ECG:   raw data. N x 4 matrix
%      Fs:    sampling rate
%      bs:    annotated fqrs series (for debug only)
%
% OUT:
%      fqrs:     foetal QRS annotation. Lq x 1 matrix
%      ch_ref:   reference channel (used for fqrs)
%      X:        feature vector
%      Xcl:      classif outcome vector (set A)
%
fqrs    = [];

N = 48+1;
b = ones(N, 1)/N; a = 1;

for k=1:size(ECG,2)
    D      = diff(ECG(:,k)).^2;
    Sk     = filtfilt(b, a, D);
    S(:,k) = Sk/std(Sk);
    [Tq, std_rr(k) rr_iqr(k) rr_med(k) ampl(k)] = MyBuildSeq(S(:,k), Fs, bs);
end
[~, ch_ref] = min(std_rr);

fqrs = MyBuildSeq(S(:,ch_ref), Fs, bs);

for k=1:size(ECG,2)
    if( bPlot )
        figure(1)
        set(gcf, 'Name', sprintf('rECG - QRS CH: %d', ch_ref));
        t = (1:size(ECG,1))/Fs;
        subplot(2,2,k)
        plot(t, ECG(:,k))
        hold on
        plot(t(bs), ECG(bs,k)+1, '.r'),
        if( isempty(fqrs)==0 )
            plot(t(fqrs), ECG(fqrs,k), '.k'),
        end
        set(gca, 'TickLength', [0 0])
        xlim([6 10])
        hold off
        figure(2)
        subplot(2,2,k)
        plot(t(1:end-1), S(:,k))
        hold on
        plot(t(bs), S(bs,k), '.r'),
        if( isempty(fqrs)==0 )
            plot(t(fqrs), S(fqrs,k), '.k'),
        end
        set(gca, 'TickLength', [0 0])
        xlim([6 10])
        hold off
    end
end
