function cmpRes=QRSdet_ann_cmp(qrsD,qrsA, diffMax, debug)
% -------------------------------------------------------------------------------------------------
%   Routine to compare QRS detection annotations
%
% --- Input parameters:
%   qrsD      estimated time position for QRS complex
%   qrsA      reference time position for QRS complex
%   diffMax   max difference to accept annotation matching
% --- Output parameters:
%   cmpRes.nTP        number of True Positives (TP) (matching annotations)
%   cmpRes.nFP        number of False Positives (FP)
%   cmpRes.nFN        number of False Negatives (FN)
%   cmpRes.meanDiff   mean of differences between TP annotations
% -------------------------------------------------------------------------------------------------
% Author: Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
% For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% -------------------------------------------------------------------------------------------------
% This program is free software; you can redistribute it and/or modify it under the terms
% of the GNU General Public License as published by the Free Software Foundation; either
% version 2 of the License, or (at your option) any later version.
%
% This program is distributed "as is" and "as available" in the hope that it will be useful,
% but WITHOUT ANY WARRANTY of any kind; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% --------------------------------------------------------------------------------------------
if(nargin<3), diffMax=0.075; end        % qrsD and qrsA are supposed to be in seconds
if(nargin<4), debug=0; end
nqrsA=length(qrsA);
nqrsD=length(qrsD);
qrsA_TP=zeros(1,nqrsA);
qrsA_diff=zeros(1,nqrsA);
nTP = 0;    % number of matching (TP)
nFN = 0;    % number of False Positive
for i = 1:nqrsA
    TP=0;
    qrsAi=qrsA(i);
    iQs = min(find(qrsD > qrsAi));
    iQp = max(find(qrsD <= qrsAi));
    if(~isempty(iQs))
        qrsDis = qrsD(iQs);
        diffQPs= qrsDis - qrsAi;
        if(diffQPs <= diffMax)
            TP=1;
            qrsA_TP(i) = qrsDis;
            qrsA_diff(i) = diffQPs;
        end
    end
    if(~isempty(iQp))
        qrsDip = qrsD(iQp);
        diffQPp= qrsAi - qrsDip;
        if(diffQPp <= diffMax && diffQPp<=diffQPs)
            TP=1;
            qrsA_TP(i) = qrsDip;
            qrsA_diff(i) = -diffQPp;
        end
    end
    if (TP)
        nTP = nTP + 1;
    else
        nFN = nFN + 1;
    end
end

nFP = max(0, nqrsD - nTP);

indTP=find(qrsA_TP>0);
if(~isempty(indTP))
    qrsA_diffTP=qrsA_diff(indTP);
    meanDiff=mean(qrsA_diffTP);
else
    meanDiff=0;   % patch
end
if(debug)
    fprintf('Number of mother QRSs= %d\n', nqrsA);
    fprintf('Number of fetal  QRSs= %d\n', nqrsD);
    
    fprintf('Number of true positives= %d\n', nTP);
    fprintf('Number of false positives= %d\n', nFP);
    fprintf('Number of false negatives= %d\n', nFN);
end

cmpRes.nTP=nTP;
cmpRes.nFP=nFP;
cmpRes.nFN=nFN;

cmpRes.meanDiff=meanDiff;

end     %== function ================================================================
%
