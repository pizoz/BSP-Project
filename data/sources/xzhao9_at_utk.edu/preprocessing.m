function [X]=preprocessing(tm,ecgs,gflag,nflag,bflag)
% preprocess ecg data
% gflag: flag indicating whether to fill the gaps
% nflat: flag indicating whether to reduce noise level
% bflag: flag indicating whether to remove baseline wandering
%
% Author:  Xiaopeng Zhao
% Updated: Xiaopeng zhao, April 25, 2013
% =====================
% to do
%   evaluate signal quality of each channel


X=ecgs;
N_chn=size(X,2); %num of channels

if gflag
    % Use interpolation to fill the gaps in the data
    ind_missing=isnan(ecgs);
    % interpolation
    for i=1:N_chn % for each channel
        pos=ind_missing(:,i);
        if sum(pos)~=0
            X(pos,i)=interp1(tm(pos==0),ecgs(pos==0,i),tm(pos),'spline');
        end
    end
end


if nflag
    % low-pass filter using 18-point moving average to suppress high frequency
    % noise
    M=9;
    for i=1:N_chn
        %X(:,i)=movAvg(X(:,i),2*M+1);
        X(:,i)=smooth(X(:,i),2*M+1,'moving');
    end
end


if bflag
    % high-pass filter using 401 point moving average to suppress baseline
    % wandering
    % the smoothing function performs better than the movAvg function in
    % tracing the baseline
    M=200; %num of points for moving average
    for i=1:N_chn
        %tmp=movAvg(X(:,i),2*M+1);
        tmp=smooth(X(:,i),2*M+1,'moving');
        X(:,i)=X(:,i)-tmp;
    end
end

% finally, we normlize the data to zero mean and std = 1
 [X,Xs]=mapminmax(X');
 X=X';
%  X=mapstd(X')';
end
