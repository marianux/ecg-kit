%% String to issue multiline system commands
% Creates the string to issue multiple system commands.
% 
% Example
% 
%   command_sep_str = sys_cmd_separation_string()
% 
% See also WFDB_command_prefix
% 
% Author: Mariano Llamedo Soria
% <matlab:web('mailto:llamedom@electron.frba.utn.edu.ar','-browser') (email)> 
% Version: 0.1 beta
% Birthdate: 5/11/2014
% Last update: 5/11/2014
% Copyright 2008-2014
% 
function command_sep_str = sys_cmd_separation_string()

    % string to change the working directory.
    if( ispc() )
        % in windows, it is important to change the
        % drive unit before the directory change
        % command (CD)
        command_sep_str = '&';
    else
        command_sep_str = ';';
    end
