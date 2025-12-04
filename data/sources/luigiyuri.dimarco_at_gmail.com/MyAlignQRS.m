function qrs = MyAlignQRS(ECG, mqrs, Fs)
%
% function: aligns QRSs on abdominal ECG 
%
% IN:
%     ECG: raw data Nx4 matrix, N=Fs*60s
%     bs: detected QRS time series [smpls]
%     Fs:  sample rate
%
% OUT:
%     qrs:  corrected mQRS time series [smpls]
%
DW      = round(0.05*Fs);
DS_10MS = round(0.01*Fs);
Le      = size(ECG,1); 
qrs     = mqrs; % preset (prior to fine-tuning)

for k=1:size(ECG,2)
    clear bc 
    l = 1;
    for j=1:size(mqrs,1)
        i1 = mqrs(j,k)-DW; i2 = mqrs(j,k)+DW;
        if(i1>0 && i2<=Le)
            xtmpl = ECG(i1:i2,k)-median(ECG(i1:i2,k));
            if(l==1)
                bc(1:i2-i1+1,l) = xtmpl;
            else
                % fine tune position
                TMPL = median(bc,2); dc_ref = -1;
                xtmpl_ref = xtmpl; DX = 0;
                for c=-DS_10MS:DS_10MS
                    if( (i1+c)>0 && (i2+c)<=Le )
                        xtmpl = ECG((i1:i2)+c,k)-median(ECG((i1:i2)+c,k));
                        dc = corr(xtmpl, TMPL);
                        if(dc>dc_ref)
                            dc_ref    = dc;
                            xtmpl_ref = xtmpl;
                            DX        = c; 
                        end
                    end
                end
                bc(1:i2-i1+1,l) = xtmpl_ref;
                qrs(j, k) = qrs(j, k)+DX; 
            end
            l=l+1;
        end
    end
end
