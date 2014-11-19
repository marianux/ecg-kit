%% List al ECGtask availables
% Look for the available ECGtask classes defined in the ECGkit\common
% folder.
%   
% Example
% 
%   list_all_ECGtask()
% 
% Arguments:
% 
% Output:
%   All ECGtask classes found and its object handles.
% 
% See also ECGwrapper
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Birthdate: 22/10/2014
% Last update: 22/10/2014
% Copyright 2008-2014
% 
function [cKnownECGtasks, cKnownECGtasksHdl] = list_all_ECGtask()
common_path = fileparts(mfilename('fullpath'));
common_path = [common_path filesep ];

ECGt_filenames = dir([common_path 'ECGtask*.m' ]);
ECGt_filenames = {ECGt_filenames(:).name};

if( isempty(ECGt_filenames) )
    error('ECGwrapper:NoTasks', 'Could not find any ECGtask in %s', common_path );
else
    % avoid the abstract class
    ECGt_filenames = ECGt_filenames(~strcmpi(ECGt_filenames, 'ECGtask.m'));
    cKnownECGtasksHdl = cellfun(@(a)( eval(a(1:end-2)) ) , ECGt_filenames, 'UniformOutput', false);
    cKnownECGtasks = arrayfun(@(a)( a{1}.name ) , cKnownECGtasksHdl, 'UniformOutput', false);
end            

if( nargout < 1 )
    fprintf(1, disp_option_enumeration( 'Task names found:\n', cKnownECGtasks));
    clear cKnownECGtasks
end
