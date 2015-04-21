%% System commands to initialize the WFDB toolbox
% Create the string prefix to execute WFDB binaries.
% 
% Example
% 
%   wfdb_prefix = WFDB_command_prefix(recording_path)
% 
% See also ECGtask_QRS_detection
% 
% Author: Mariano Llamedo Soria
% <matlab:web('mailto:llamedom@electron.frba.utn.edu.ar','-browser') (email)> 
% Version: 0.1 beta
% Birthdate: 1/1/2014
% Last update: 31/10/2014
% Copyright 2008-2015
% 
function wfdb_prefix = WFDB_command_prefix(recording_path)

    command_sep_str = sys_cmd_separation_string();
    
    % string to change the working directory.
    if( ispc() )
        % in windows, it is important to change the
        % drive unit before the directory change
        % command (CD)
        wfdb_prefix = recording_path;
        wfdb_prefix = [wfdb_prefix(1:2) command_sep_str ];
    else
        wfdb_prefix = [];
    end
    wfdb_prefix = [wfdb_prefix 'cd ' recording_path command_sep_str ];
