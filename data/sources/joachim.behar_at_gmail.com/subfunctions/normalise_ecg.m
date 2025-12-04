function [necg,scalings] = normalise_ecg(ecg,start,stop)
% normalisation function. See STEP 1-3 in the code below to understand 
% what it does.
% 
% inputs:
%   ecg: ecg before normalization
%   start/stop:  starting and ending points (in samples) defining
%   representative ecg segments on which to compute scalings and shifts.
%   The whole ecg is not used for that otherwise the scaling would be
%   dependant on the local artefacts etc.
%   
% outputs:
%   necg: the normalized ecg which values will be within [0 1].
%   scalings: scaling factor by which the ecg(s) have been multiplied 
%   (NOTE: that it does not take the detrand or tanh operations into account).
%
%
% FECG extraction toolbox, version 1.0, Sept 2013
% Released under the GNU General Public License
%
% Copyright (C) 2013  Joachim Behar
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2013
% joachim.behar@eng.ox.ac.uk
%
% Last updated : 28-08-2013
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

% == manage inputs
if size(ecg,2)>size(ecg,1)
    ecg = ecg';
end

% == core function
try
    % STEP 1: Normalizes a multivariate dataset data (where each data vector is a row in data) 
    % by scaling and shifting such that in each column the min/max becomes 0/1.
    % If the min and the max in some column coincide, this column is not
    % changed.
    [~,scalings,~] = normalizeData01(ecg(start:stop,:));

    % STEP 2: detrend the ecg
    detrendedInput = detrend(bsxfun(@times,ecg,scalings),'constant');
        
    % STEP 3: take the tanh of the sigal . This step is applied in order to
    % avoid outliers whih could result in the reservoir state or LMS
    % weights to take unexpected values that are not covered by the normal
    % trajectory given the optimized global parameters or output learned.
    necg = tanh(detrendedInput);
catch ME
    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
    necg = ecg;
end

end


function [normData, scalings, shifts] = normalizeData01(data)
% Normalizes a multivariate dataset data (where each data vector is a row in data) 
% by scaling and shifting such that in each column the min/max becomes 0/1.
% If the min and the max in some column coincide, this column is not
% changed. 
% 
% Input arg: 
% data: a dataset, either real-valued array of size N by dim or a cell array of size
%    [nrSamples, 1], where the i-th cell is a real-valued array of size N_i by dim 
%
% Outputs:
% normData: the dataset normalized to columns with min/max = 0/1. Each
%    column in normData is computed from the corresponding column in data by 
%    normalizedColumn = scalefactor * (originalColum + shiftconstant). If
%    the input is a cell structure, the same scalefactors and shiftconstants are
%    applied across all cells, such that the *global* min/max of normData
%    becomes 0/1.
% scalings: a row vector of lenght dim giving the scalefactors
% shifts: a row vector of lenght dim giving the shiftconstants
%
% Created by H. Jaeger, June 21, 2006

if isnumeric(data)
    dim = size(data,2);
    mins = min(data); maxs = max(data);
    scalingsInv = maxs - mins;
    scalings = ones(1,dim);
    shifts = zeros(1,dim);
    normData = data;
    for d = 1:dim
        if scalingsInv(1,d) > 0
            scalings(1,d) = 1/scalingsInv(1,d);
            shifts(1,d) = -mins(1,d);
            normData(:,d) = (data(:,d) + shifts(1,d)) * scalings(1,d);    
        end
    end
elseif iscell(data)
    dim = size(data{1,1},2);
    nrSamples = size(data,1);
    %check if all cells have same dim
    for n = 1:nrSamples
        if size(data{n,1},2) ~= dim
            error('all cells must have same row dim');
        end
    end
    mins = min(data{1,1});
    maxs = max(data{1,1});
    for n = 1:nrSamples
        mins = min(mins, min(data{n,1}));
        maxs = max(maxs, max(data{n,1}));
    end
    scalingsInv = maxs - mins;
    scalings = ones(1,dim);
    shifts = zeros(1,dim);
    normData = data;
    for d = 1:dim
        if scalingsInv(1,d) > 0
            scalings(1,d) = 1/scalingsInv(1,d);
            shifts(1,d) = -mins(1,d);
            for n = 1:nrSamples
                normData{n,1}(:,d) = (data{n,1}(:,d) + shifts(1,d)) * scalings(1,d);  
            end
        end
    end   
else error('input data must be array or cell structure');
end
end
    