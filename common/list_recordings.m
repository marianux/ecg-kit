function rec_names = list_recordings(this_path)
%% Find recordings of known extension in a given path
%   
% Example
% 
%   rec_names = list_recordings(this_path)
% 
% Arguments:
%      +this_path: [cell] REQUIRED
%           
%           The path to look into for recordings.
% 
% 
% Output:
%     + rec_names: The names of the recordings.
% 
% See also ECGtask_QRS_detection, ECGwrapper
% 
% Author: Mariano Llamedo Soria llamedom@frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 28/11/2018
% Birthdate  : 28/11/2018
% Copyright 2008-2018

cKnownFormatExtensions = {'.dat', '.ecg', '.mat', '.xml', '.hes', '.dat', };

aux_recs = cellfun( @(a)(dir([this_path '*' a] ) ), cKnownFormatExtensions, 'UniformOutput', false );

aux_recs = cell2mat(aux_recs(:));

rec_names = unique({aux_recs.name});

