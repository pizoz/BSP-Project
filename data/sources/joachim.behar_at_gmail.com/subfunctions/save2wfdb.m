function status = save2wfdb(name,fqrs,path,path2WFDBlib)    
% this function saves the fqrs fiducials returned by the algorithm in WFDB
% anns format ans .txt. Output is saved into path. Path must contain the
% corresponding name.hea for creating the WFDB annotation files.
% 
% inputs
%   name: name of the record that is being processed
%   fqrs: fqrs fiducial position (in seconds)
%   path: path
%   path2WFDBlib: path 2 wfdb library to make sure the WFDB function are
%   accepcible from path
%
% output
%   status: 1 for success and 0 for failure
%
% Safe Foetus Monitoring Toolbox, version 1.0, Sept 2013
% Released under the GNU General Public License
%
% Copyright (C) 2013  Joachim Behar
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2013
% joachim.behar@eng.ox.ac.uk
%
% Last updated : 02-08-2013
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

try
    mainDir = pwd;    
    cd(path)
        addpath(genpath(path2WFDBlib));
        % write QRS annotation to txt file in ms
        dlmwrite(strcat(name,'.txt'),round(1000*fqrs'));        
        % then read it
        ann = dlmread([name '.txt']);
        % then write in WFDB format
        wrann(name,'fqrsmyalgo',ann);
    cd(mainDir)
        status = 1;
catch ME
    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
    status = 0;
end

end



