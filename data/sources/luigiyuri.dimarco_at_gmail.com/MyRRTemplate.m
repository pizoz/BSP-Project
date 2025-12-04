function rECG = MyRRTemplate(ECG, mqrs, Fs, bs)
%
% function: generates template rof RR interval
%
% IN:
%     ECG: raw data Nx4 matrix, N=Fs*60s
%     mqrs: detected maternal QRS time series [smpls]
%     Fs:  sample rate
%     bs:  foetal qrs annotations (for plot only)
%
% OUT:
%
DW      = round(0.08*Fs);
Le      = size(ECG,1); 

rECG = ECG;

% align abdominal ecg mQRS
mqrs = MyAlignQRS(ECG, mqrs, Fs);

for k=1:size(ECG,2)
    clear bc len vi1
    l = 1;
    for j=1:numel(mqrs)
        i1 = mqrs(j)-DW; 
        i2 = mqrs(j)+DW;
        if(i1>0 && i2<=Le)
            bc(1:i2-i1+1,l) = ECG(i1:i2,k)-median(ECG(i1:i2,k));
            len(l) = i2-i1+1; % length of RR interval 
            vi1(l) = i1; % 1st sample index
            l=l+1;
        end
    end

    % calculate pair-wise distance on qrs window
    Y  = pdist(bc', 'correlation'); B = squareform(Y); 
    for q=1:numel(vi1)
        i1 = vi1(q); 
        i2 = vi1(q)+2*DW;
        
        % select NB beats with most similar length
        NB = 10;
        [~, Is] = sort(B(q,:), 'ascend');
        TMPL = median(bc(1:i2-i1+1, Is(1:NB)),2);
        rECG(i1:i2, k) = ECG(i1:i2, k) - TMPL; % MyTimeAlign(ECG(i1:i2,k), TMPL);
        
%         plot(ECG(:, k)), 
%         hold on,
%         plot(rECG(:, k), 'k'), 
%         plot(mqrs, ECG(mqrs,k), '.m', 'MarkerSize', 14),
%         plot(bs, ECG(bs,k), '.r', 'MarkerSize', 14),
%         hold off
%         xlim([0 2]*Fs), ylim([-40 40])
%         set(gca, 'TickLength', [0 0])
    end
    
%     subplot(2,2,k)
%     t = (1:size(rECG,1))/Fs;
%     plot(t, rECG(:, k)), 
%     hold on,
%     plot(mqrs/Fs, rECG(mqrs,k), '.r'),
%     hold off
%     xlim([6 10])
%     set(gca, 'TickLength', [0 0])
end
