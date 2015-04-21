%% Checks administrator privileges 
% This script check either if the Matlab session was executed with
% administrator privileges. It is used during the installation to produce
% all-users changes.
% 
% Example
% 
%   bHasPrivs = HasAdminPrivs()
% 
% See also InstallECGkit, UnInstallECGkit, isMatlab, isOctave
% 
% Author: Mariano Llamedo Soria
% <matlab:web('mailto:llamedom@electron.frba.utn.edu.ar','-browser') (email)> 
% Version: 0.1 beta
% Birthdate: 31/10/2014
% Last update: 31/10/2014
% Copyright 2008-2015
% 
function bHasPrivs = HasAdminPrivs(bAskUserForPrivs)

    if( nargin < 1 || isempty(bAskUserForPrivs) )
        bAskUserForPrivs = true;
    end

    bHasPrivs = false;
    
    % get the correct command names in this architecture.
    sys_cmds = sys_command_strings();

    pathdef_filename = fullfile(matlabroot,'toolbox', 'local', 'pathdef.m');
    pathdef_dummy_filename = [pathdef_filename '.delete_me'];

    % copy and delete
    if( ispc() )
        str_command = [ 'is_admin.bat "' pathdef_filename '" "' pathdef_dummy_filename '"' ];
    else
        str_command = [ '. is_admin.sh ' pathdef_filename ' ' pathdef_dummy_filename ];
    end

    if( bAskUserForPrivs )
        [status,~] = system( str_command, '-runAsAdmin' );
    else
        [status,~] = system( str_command );
    end
    
    if( status == 0 )
        bHasPrivs = true;
    end        
    
end
