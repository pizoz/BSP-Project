function estimated_data=fillMissed2(ns,ne,fwdata,bwdata,defval)
% estimates the missing values in the signal using a forward and backward
% AR model
% Inputs:
% ns: first index of missing data interval
% ne: last index of missing data interval
% fwdata: samples before the missing part
% bwdata: samples after the missing part 
% Output: estimated_data
% This code is written based on the estimation algorithm of (Esquef et al., 2006)

G=ne-ns+1;%length of missing data interval
fwdata_diff=diff(fwdata);
if length(fwdata_diff)==0
    fwdata_diff=NaN;
end

bwdata_diff=diff(bwdata);
if length(bwdata_diff)==0
    bwdata_diff=NaN;
end 


fwdata_diff_mean=mean(fwdata_diff);
bwdata_diff_mean=mean(bwdata_diff);
interval=bwdata(end)-fwdata(end);


% while floor(interval/(G+1))>min(fwdata_diff_mean,-bwdata_diff_mean)
% while floor(interval/(G+1))>max([fwdata_diff,-bwdata_diff])
while floor(interval/(G+1))-max([fwdata_diff,-bwdata_diff])>0.05*max([fwdata_diff,-bwdata_diff])
    G=G+1;
end

if abs(interval/(G+1)-max([fwdata_diff,-bwdata_diff])) > abs(interval/G-max([fwdata_diff,-bwdata_diff]))
    G=G-1;
end


flag=0;
round=1;
while ~flag && round<=3
    if G==1
        estimated_data=floor(1/2*(fwdata(end)+bwdata(end)));
    else

        % % Forward data estimation
        % if N1>1
        %     fwcoeffs=arburg(fwdata,p);
        %     fwexcit=filter(fwcoeffs,1,fwdata);
        %     fwexcit=[fwexcit [fwexcit(N1:-1:max(N1-G+1,1)) fwexcit(1)*ones(1,max(G-N1,0))]];
        %     fwdata_extend=filter(1,fwcoeffs,fwexcit);
        % else
        %     fwdata_extend=fwdata(1)*ones(1,N1+G);
        % end
% % % %         if ~isnan(fwdata_diff_mean) && ~isnan(bwdata_diff_mean)
% % % %             fwdata_diff_estimate=fwdata_diff_mean*ones(1,G);
% % % %             bwdata_diff_estimate=-bwdata_diff_mean*ones(1,G);
% % % %         elseif isnan(fwdata_diff_mean) && ~isnan(bwdata_diff_mean)
% % % %             bwdata_diff_estimate=-bwdata_diff_mean*ones(1,G);
% % % %             fwdata_diff_estimate=bwdata_diff_estimate;
% % % %         elseif ~isnan(fwdata_diff_mean) && isnan(bwdata_diff_mean)
% % % %              fwdata_diff_estimate=fwdata_diff_mean*ones(1,G);
% % % %              bwdata_diff_estimate=fwdata_diff_estimate;
% % % %         else 
% % % %             fwdata_diff_estimate=defval*ones(1,G);
% % % %             bwdata_diff_estimate=fwdata_diff_estimate;
% % % %         end

        %
        % % Backward data estimation
        % if N2>1
        %     bwcoeffs=arburg(bwdata,p);
        %     bwexcit=filter(bwcoeffs,1,bwdata);
        %     bwexcit=[bwexcit [bwexcit(N2:-1:max(N2-G+1,1)) bwexcit(1)*ones(1,max(G-N2,0))]];
        %     bwdata_extend=filter(1,bwcoeffs,bwexcit);
        % else
        %     bwdata_extend=bwdata(1)*ones(1,N2+G);
        % end

        
        
        % Crossfading forward and backkward estimates
% % % %         W=1:-1/(G-1):0;
% % % % 
% % % %         estimated_data_diff=W.*fwdata_diff_estimate+(1-W).*bwdata_diff_estimate;
% % % %         estimated_data=floor(fwdata(end)+estimated_data_diff*triu(ones(G)));
        estimated_data_diff= interval/(G+1);
        estimated_data=floor(fwdata(end)+estimated_data_diff*[1:G]);
%         estimated_data=floor(W.*(fwdata(end)*ones(1,G)+estimated_data_diff*triu(ones(G)))+...
%         (1-W).*(bwdata(end)*ones(1,G)-estimated_data_diff*rot90(triu(ones(G)))));
    end
    data_test=[fwdata(end) estimated_data bwdata(end)];
%     if any(diff(data_test)<150) & G>1
    if G>1 && any(diff(data_test)<0.5*min([fwdata_diff,-bwdata_diff]))
        G=G-1;
        round=round+1;
    elseif any(diff(data_test)>1.5*max([fwdata_diff,-bwdata_diff]))
        G=G+1;
        round=round+1;
    else
        flag=1;
    end
end

