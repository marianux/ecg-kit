%% Uninstallation script of ECGkit
% Script to install this toolbox.
% 
% ****************************************************************
% ** Check the pretty printed help documentation in help folder **
% ****************************************************************
% 
% Example
%   UnInstallECGkit()
% 
% 
% See also InstallECGkit
% 
% Author: Mariano Llamedo Soria
% <matlab:web('mailto:llamedom@electron.frba.utn.edu.ar','-browser') (email)> 
% Version: 0.1 beta
% Birthdate: 01/09/2012
% Last update: 18/10/2014
% Copyright 2008-2015

function UnInstallECGkit( bIgnoreAdminPrivs )

    if( nargin < 1 || isempty(bIgnoreAdminPrivs) )
        bIgnoreAdminPrivs = false;
    end        

    release_str = 'v0.1.0 beta - 05/05/2015';
    backup_string = 'ECGkit_backup';

    % get the correct command names in this architecture.
    sys_cmds = sys_command_strings();
    
    bMatlab = isMatlab();
    bOctave = isOctave();
    
    bUseDesktop = usejava('desktop');

    
    if( bUseDesktop || bOctave )
    
        disp_string_framed('*Blue', sprintf('ECGkit for Matlab %s', release_str) );
        
        if( bMatlab )
        
          % Tab completion of selected functions
          % Modify the TC.xml and TC.xsd
          bTabUnInstallError = false;
          
          tcxml_filename = fullfile(matlabroot,'toolbox', 'local', 'TC.xml');
          tcxsd_filename = fullfile(matlabroot,'toolbox', 'local', 'TC.xsd');
          
          if( exist([tcxml_filename backup_string ], 'file') )
              str_command = [ sys_cmds.copy_command ' "' tcxml_filename '" "' tcxml_filename backup_string '"' sys_cmds.command_sep_str ];
          else
              str_command = [];
          end

          if( exist([tcxsd_filename backup_string ], 'file') )
              str_command = [ str_command sys_cmds.copy_command ' "' tcxsd_filename '" "' tcxsd_filename backup_string '"' ];
          end
          
          if( isempty(str_command) )
              status = 1;
          else
              [status, ~] = system( str_command, '-runAsAdmin' );
          end
          
          if( status == 0 )
              cprintf('*blue', 'Remember to restart Matlab for function''s "tab completion" changes to take effect.\n' );
  %             cprintf('SystemCommands', '\nIn case something went wrong regarding "tab completion" feature, you can restore the original files:\n%s\n%s\n\nTo the original versions, just removing the ', [tcxml_filename backup_string], [tcxsd_filename backup_string] );
  %             cprintf('*red', '"%s"', backup_string );
  %             cprintf('SystemCommands', ' extension.\n\n' );
          else                    
              disp_string_framed('*red', 'Could not uninstall tab completion');
              cprintf('SystemCommands', 'to uninstall "tab completion" capability it is needed to run this script with administrator/root privileges, or say ' );
              cprintf('*red', '"YES"' );
              cprintf('SystemCommands', ' when.\nRun Matlab as administrator and execute this script again.\n\n' );
              if( ~bIgnoreAdminPrivs )
                  return
              end
          end
          
        end
        
    end
    
    %path related constants.
    root_path = [fileparts(mfilename('fullpath')) filesep ];

    if( bUseDesktop || bOctave )
        fprintf(1, 'Removing paths. ' );    
    end

    %path related constants.
    default_paths = { ...
                        [ root_path ';' ]; ...
                        [ root_path 'common' filesep ';' ]; ...
                        [ root_path 'help' filesep ';' ]; ...
                        [ root_path 'common' filesep 'a2hbc' filesep ';' ]; ...
                        [ root_path 'common' filesep 'a2hbc'  filesep 'scripts' filesep ';' ]; ...
                        [ root_path 'common' filesep 'export_fig' filesep ';' ]; ...
                        [ root_path 'common' filesep 'wavedet' filesep ';' ]; ...
                        [ root_path 'common' filesep 'plot2svg' filesep ';' ]; ...
                        [ root_path 'common' filesep 'ppg' filesep ';' ]; ...
                        [ root_path 'common' filesep 'prtools' filesep ';' ]; ...
                        [ root_path 'common' filesep 'prtools_addins' filesep ';' ]; ...
                        [ root_path 'common' filesep 'bin' filesep ';' ]; ...
                        [ root_path 'common' filesep 'kur' filesep ';' ]; ...
                        [ root_path 'common' filesep 'LIBRA' filesep ';' ]; ...
                    };

    cellfun(@(a)(rmpath(a)), default_paths);    
    
    if( bUseDesktop || bOctave )

        bOk = savepath();
    
        if( bOk == 0 )
            
            fprintf(1, 'done !\n' );

            fprintf(1, 'Thanks for trying the kit.\n\nNow you can safely delete %s\n', root_path );
        
        else
           
            if( ispc() )
                home_path = [getenv('HOMEDRIVE') getenv('HOMEPATH') ];
            elseif( isunix() )
                home_path = getenv('HOME');
            elseif( ismac() )
                home_path = getenv('HOME');
            end

            if( home_path(end) ~= filesep )
                home_path = [home_path filesep];
            end

            bOk = savepath([home_path 'pathdef.m'] );
            
            if( bOk ~= 0 )
                disp_string_framed(2, 'Path Not Saved')
                fprintf(2, 'Save path manually to uninstall the ECG-Kit permanently. Deletes all occurrences of %s in your path.\n', root_path);
            end
            
        end
        
        
    end
    
end
