function out = findFQRS(FECG)
% FQRS detection based on entropy estimation and peak detection
%   (c) Jakub Kuzilek
% 	http://bio.felk.cvut.cz/~kuziljak/     E-mail: jakub.kuzilek@gmail.com
%   Version: 1.0                    Last update:   23/08/2013.
% 					                (Version: 1.0, 23/08/2013)
%
%======================================================
%
% PURPOSE:    This function detects FQRS in ECG recording. Developed for
%             purpose of Physionet Challenge 2013.
%
% MANDATORY INPUT ARGUMENTS
%   FECG ..... fetal ECG Mx1, M - length of data (one lead)
% OPTIONAL INPUT ARGUMENTS
%   none
% OUTPUT ARGUMENTS
%   out .... positions of FQRS
%=========================================================
ent = [];
M = 5;
   for n = 1:length(FECG)-M
       ent(end+1) = wentropy(FECG(n:n+M),'shannon');
   end

[~,out]=findpeaks(-1*ent,'MINPEAKHEIGHT',200,'MINPEAKDISTANCE',300);