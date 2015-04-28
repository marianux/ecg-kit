%% Installation script of ECGkit
% Script to install this toolbox.
% 
% ****************************************************************
% ** Check the pretty printed help documentation in help folder **
% ****************************************************************
% 
% Example
% 
%   InstallECGkit()
% 
% or
% 
%   InstallECGkit(false) 
% 
% to avoid checking online for updates.
% 
% See also UnInstallECGkit
% 
% Author: Mariano Llamedo Soria
% <matlab:web('mailto:llamedom@electron.frba.utn.edu.ar','-browser') (email)> 
% Version: 0.1 beta
% Birthdate: 01/09/2012
% Last update: 18/10/2014
% Copyright 2008-2015
% 
function InstallECGkit(bIgnoreUpdate)
    
    backup_string = 'ECGkit_backup';

    if( nargin < 1 || ~islogical(bIgnoreUpdate) )
        bIgnoreUpdate = false;
    end

    fprintf(1, 'Adding paths\n' );
    
    %path related constants.
    root_path = [fileparts(mfilename('fullpath')) filesep ];
    
    %path related constants.
    default_paths = { ...
                        [ root_path ';' ]; ...
                        [ root_path 'help' filesep ';' ]; ...
                        [ root_path 'common' filesep ';' ]; ...
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

    cellfun(@(a)(addpath(a)), default_paths);    
    
    % get the correct command names in this architecture.
    sys_cmds = sys_command_strings();
    
    ECGkitURL_str = 'http://code.google.com/p/ecg-kit/wiki/LatestVersion';
    % Version info.
    release_str = 'v0.1 beta - 22/07/2014';

    % URL check for updates

    % Check admin privileges
    bAdminPrivs = HasAdminPrivs(false);
    
    if( ~bAdminPrivs )
        disp_string_framed('*red', 'No admin/root privileges');
        cprintf('SystemCommands', 'To install the ECGkit it is needed to run this script with administrator/root privileges.\nRun Matlab as administrator or ' );
        cprintf('*red', 'answer "YES"' );
        cprintf('SystemCommands', ' when the system asks for permission.\n\n' );
    end
    
    % Check if the GUI is up and working
    bUseDesktop = usejava('desktop');
    bMatlab = isMatlab();
    bOctave = isOctave();
    
    if( bMatlab )
        matlab_ver = ver('Matlab');

        if( str2double(matlab_ver.Version) < 8 )
            cprintf('Red', 'Check updates manually since the Matlab version is too old to do it automatically.\n' );
            bIgnoreUpdate = true;
        end

    elseif( bOctave )
        disp_string_framed([1,0.5,0], 'Octave is not fully supported');
        cprintf('Red', 'Although you can find problems using the toolbox in Octave, I am trying to make it fully compatible. Please report any problem that you could find.\n\n' );
    end
    
    if ~bIgnoreUpdate
        
        fprintf(1, 'Checking for updates\n' );
        if( bMatlab )
          [strResponse, URLstatus] = urlread( ECGkitURL_str,'TimeOut',5);
        elseif( bOctave )
          [strResponse, URLstatus] = urlread( ECGkitURL_str);
        end
          
        if( URLstatus == 1 )
            aux_idx = strfind(strResponse, '###');
            latest_release_str = strResponse(aux_idx(1)+3:aux_idx(2)-1);
            
            if( strcmpi(latest_release_str, release_str) )
                disp_string_framed('Blue', 'ECGkit is up to the date');
            else
                disp_string_framed('*Red', sprintf('ECGkit %s version available', latest_release_str) );
                fprintf(1, 'Consider updating %s.\n', '<a href = "https://code.google.com/p/ecg-kit/">here</a>' );
                
                cellfun(@(a)(rmpath(a)), default_paths);    
                return
            end
            
        end
    end

    fprintf(1, 'Compiling sources\n' );
    
    %Check compilation of source MEX files
    common_path = [ root_path 'common' filesep];
    source_files = dir([ common_path '*.c'] );
    lsource_files = length(source_files);
    
    if(bOctave)
        prev_folder = cd(common_path);
    end
        
    for ii = 1:lsource_files
        [~, source_file_name] = fileparts( source_files(ii).name);
        mex_file = dir([common_path  source_file_name '.' mexext ]);
        if( isempty(mex_file) || mex_file.datenum <= source_files(ii).datenum  )
        
            if( bMatlab )
              eval(['mex -outdir ''' common_path ''' ''' [common_path source_files(ii).name] '''']);
            elseif(bOctave)
              eval(['mkoctfile --mex ''' [common_path source_files(ii).name] '''']);
            end
        end
    end

    if(bOctave)
        cd(prev_folder);
    end
    
    if( bUseDesktop || bOctave )

        % Tab completion of selected functions
        % Modify the TC.xml and TC.xsd
        fprintf(1, 'Installing tab-completion feature.\n' );
        
        bTabInstallError = false;
        fid = 0;

        tcxml_filename = fullfile(matlabroot,'toolbox', 'local', 'TC.xml');
        tcxsd_filename = fullfile(matlabroot,'toolbox', 'local', 'TC.xsd');
        
        tmp_path = tempdir();

        my_tcxml_filename = [tmp_path 'TC.xml' ];
        my_tcxsd_filename = [tmp_path 'TC.xsd' ];
        
        status = copyfile(tcxml_filename, my_tcxml_filename, 'f');
        status = status & copyfile(tcxsd_filename, my_tcxsd_filename, 'f');
        
        if( status )

            if( exist([tcxml_filename backup_string ], 'file') )
                str_command = [];
            else
                str_command = [ sys_cmds.copy_command ' "' tcxml_filename '" "' tcxml_filename backup_string '"' sys_cmds.command_sep_str ];
            end

            if( ~exist([tcxsd_filename backup_string ], 'file') )
                str_command = [ str_command sys_cmds.copy_command ' "' tcxsd_filename '" "' tcxsd_filename backup_string '"' sys_cmds.command_sep_str ];
            end

            try 

                last_binding = [];
                bStartFound = false;
                bFinishFound = false;
                fid = fopen(my_tcxml_filename, 'r+');

                if( fid > 0 )

                    while( ~feof(fid) && ~bFinishFound ) 

                        str_aux = fgetl(fid);

                        if( ~bStartFound && ~isempty(strfind( str_aux, '<!-- ECGkit start -->')) )
                            bStartFound = true;
                            break
                        elseif( ~isempty(strfind( str_aux, '</TC>')) )
                            bFinishFound = true;
                            break
                        elseif( ~isempty(strfind( str_aux, '</binding>')) )
                            last_binding = ftell(fid);
                        end
                    end

                    if( ~bStartFound && bFinishFound && ~isempty(last_binding) )

                        fseek(fid, last_binding, 'bof');

                        our_fid = fopen( [ root_path 'common' filesep 'TC.xml_additions.xml' ] , 'r');

                        if( our_fid > 0 )
                            fwrite(fid, fread(our_fid, inf, '*uint8') );
                            fclose(our_fid);
                        end
                    end

                    fclose(fid);
                else
                    bTabInstallError = true;
                end

            catch
                if( fid > 0 )
                    fclose(fid);
                end
                bTabInstallError = true;
            end

            try 

                start_idx = [];
                bStartFound = false;
                fid = fopen(my_tcxsd_filename, 'r+');

                if( fid > 0 )

                    while( ~feof(fid) ) 

                        last_binding = ftell(fid);
                        str_aux = fgetl(fid);
                        if( ~bStartFound && ~isempty(strfind( str_aux, '<!-- ECGkit start')) )
                            bStartFound = true;
                            break
                        elseif( isempty(start_idx) && ~isempty(strfind( str_aux, '<xsd:pattern value="-[A-Za-z0-9]+" />')) )
                            fseek(fid, 0, 'bof');
                            str_aux_pre = rowvec(fread(fid, last_binding, '*char'));
                            str_aux = rowvec(fgetl(fid));
                            str_aux_post = rowvec(fread(fid, inf, '*char'));
                            fseek(fid, 0, 'bof');
                            fwrite(fid, [str_aux_pre sprintf('<!-- ECGkit start\n') str_aux sprintf('\nECGkit end -->\n	    <xsd:pattern value="[-_A-Za-z0-9]+" />\n') str_aux_post], 'char' );
                        end
                    end

                    fclose(fid);

                else
                    bTabInstallError = true;
                end

            catch
                if( fid > 0 )
                    fclose(fid);
                end
                bTabInstallError = true;
            end         

            if( bTabInstallError)
                disp_string_framed('*red', 'Could not setup tab completion');
                cprintf('red', 'Some error happened during the tab completion installation. Do it manually or ask for help.\n' );
            else
            
                % apply changes
                str_command = [ str_command ... 
                                sys_cmds.copy_command ' "' my_tcxml_filename '" "' tcxml_filename '"' sys_cmds.command_sep_str ... 
                                sys_cmds.copy_command ' "' my_tcxsd_filename '" "' tcxsd_filename '"'... 
                                ];

                [status, ~] = system( str_command, '-runAsAdmin' );
                if( status == 0 )
                    cprintf('*blue', 'Remember to restart Matlab for function''s "tab completion" changes to take effect.\n' );
                else                    
                    bTabInstallError = true;
                end

            end
            
        end
        
        bOk = savepath();
    
        if( bOk == 0 )

            fprintf(1, 'done !\n' );
            
            disp_string_framed('*Blue', sprintf('ECGkit for Matlab %s', release_str) );
            
            str_aux = [ 'examples' filesep 'examples.m'];
            
            fprintf(1, 'Kit was correctly installed.\n\nYou can start reading the %s, or if you prefer, trying these %s.\nGo to the %s if you need help.\n', ...
                '<a href = "matlab: web(''http://ecg-kit.readthedocs.org/en/master/'', ''-browser'' )">documentation</a>', ...
                [ '<a href = "matlab: opentoline(' str_aux ',1)">examples</a>' ], ...
                '<a href = "matlab: web(''https://groups.google.com/d/forum/ecgkit-users'', ''-browser'' )">forum</a>');
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

            bOk = savepath([h   ome_path 'pathdef.m'] );
            
            if( bOk ~= 0 )
                disp_string_framed(2, 'Path Not Saved')
                fprintf(2, 'Save path manually to install the ECG-Kit permanently.\n');
            end
            
        end

    end
    
end



