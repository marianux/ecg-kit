%% Strings to execute typical I/O commands via system calls
% Creates a struct with the strings to issue system calls.
% 
% Example
% 
%   sys_cmds = sys_command_strings()
% 
% See also HasAdminPrivs, InstallECGkit
% 
% Author: Mariano Llamedo Soria
% <matlab:web('mailto:llamedom@electron.frba.utn.edu.ar','-browser') (email)> 
% Version: 0.1 beta
% Birthdate: 5/11/2014
% Last update: 5/11/2014
% Copyright 2008-2015
% 
function sys_cmds = sys_command_strings()


    if( ispc() )
        % in windows, it is important to change the
        % drive unit before the directory change
        % command (CD)
        sys_cmds.command_sep_str = '&';
        sys_cmds.delete_command = 'del';
        sys_cmds.copy_command = 'copy /Y';
        sys_cmds.move_command = 'move /Y';
    else
        sys_cmds.delete_command = 'rm';
        sys_cmds.copy_command = 'cp -f';
        sys_cmds.move_command = 'mv -f';
        sys_cmds.command_sep_str = ';';
    end
