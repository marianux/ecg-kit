%% Installation script of ecg-kit
% Script to install this toolbox.
% 
% ****************************************************************
% ** Check the pretty printed help documentation in help folder **
% ****************************************************************
% 
% Example
% 
%   Installecg-kit()
% 
% or
% 
%   Installecg-kit(false) 
% 
% to avoid checking online for updates.
% 
% See also UnInstallecg-kit
% 
% Author: Mariano Llamedo Soria
% <matlab:web('mailto:llamedom@electron.frba.utn.edu.ar','-browser') (email)> 
% Version: 0.1 beta
% Birthdate: 01/09/2012
% Last update: 18/10/2014
% Copyright 2008-2015
% 
function InstallECGkit(bIgnoreUpdate)
    

    % MAC compatibility alert
%     if( ismac() )
%         disp_string_framed('*red', 'MAC is not fully supported at the moment');
%         fprintf(1, 'If you think you can help, please contact us through the %s at the project %s.\n\n', ...
%               '<a href = "matlab: web(''https://groups.google.com/forum/#!forum/ecg-kit-users'', ''-browser'' )">forum</a>', ...
%               '<a href = "matlab: web(''http://marianux.github.io/ecg-kit/'', ''-browser'' )">ecg-kit web''s page</a>' );
%         return
%     end

    backup_string = 'ecg-kit_backup';

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
    
%     ECGkit_URL_str = 'http://ecg-kit.readthedocs.org/en/master/';
    ECGkit_URL_str = 'https://github.com/marianux/ecg-kit/releases';

    % Version info.
    
    fid = fopen([ root_path 'LatestVersion' ]);
    if fid < 0
        % version check not possible - force update message
        this_release_str  = 'v0.1.1 beta - 05/05/2015';
    else
        try
            this_release_str  = fgetl(fid);
            fclose(fid);
        catch MException
            fclose(fid);
            throw(MException);
        end
    end

    % URL check for updates

    % Check admin privileges
    bAdminPrivs = HasAdminPrivs(false);
    
    if( ~bAdminPrivs )
        disp_string_framed('*red', 'No admin/root privileges');
        cprintf('SystemCommands', 'To install tab-completion feature with the ecg-kit, it is suggested to run this script with administrator/root privileges.\nRun Matlab as administrator or ' );
        cprintf('*red', 'answer "YES"' );
        cprintf('SystemCommands', ' when the system asks for permission. Otherwise, tab completion will not be available.\n\n' );
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
        disp_string_framed('[1,0.5,0]', 'Octave is not fully supported');
        cprintf('Red', 'Although you can find problems using the toolbox in Octave, I am trying to make it fully compatible. Please report any problem that you could find.\n\n' );
        
        % Special functions not compatible with Matlab

        octave_files = { ...
                          'octave -ECGt.m' 'ECGtask.m' ...
                        };
                        
        for ii = 1:size(octave_files,1)                  
          
          octave_file_src = fullfile(root_path,'common', octave_files{ii,1});
          octave_file_dst = fullfile(root_path,'common', octave_files{ii,2});
          
          status = copyfile(octave_file_src, octave_file_dst, 'f');
            
          if( status == 0 )

            disp_string_framed('*Red', 'Installation Error' );
            fprintf(1, 'Failed copying "%s"  -->>  "%s".\n', octave_file_src, octave_file_dst );
            cellfun(@(a)(rmpath(a)), default_paths);    
            return
          
          end      
          
        end        
    end
    
    if ~bIgnoreUpdate
        
        fprintf(1, 'Checking for updates\n' );
        URLstatus = 0;
        
        if( bMatlab )
          [strResponse, URLstatus] = urlread( ECGkit_URL_str,'TimeOut',5);
        elseif( bOctave )
          [strResponse, URLstatus] = urlread( ECGkit_URL_str);
        end
          
        if( URLstatus == 1 )
            aux_pattern_pre = '/marianux/ecg-kit/tree/';
            aux_idx = strfind(strResponse, aux_pattern_pre);
            
            if( ~isempty(aux_idx) )
                aux_pattern_pos = '" ';
                aux_idx2 = strfind(strResponse(aux_idx(1):end), aux_pattern_pos);
                
                if( ~isempty(aux_idx2) )
                
                    online_latest_release_str = strResponse(aux_idx(1)+(length(aux_pattern_pre):(aux_idx2(1)-2)) );
                    
                    if( islater(online_latest_release_str, this_release_str) )
                      disp_string_framed('*Red', sprintf('ecg-kit %s version available', online_latest_release_str) );
                      fprintf(1, 'Consider updating %s.\n', '<a href = "matlab: web(''http://marianux.github.io/ecg-kit/'', ''-browser'' )">here</a>' );

                      % dont abort installation, just warn.
%                       cellfun(@(a)(rmpath(a)), default_paths);    
%                       return
                    end
                    
                end
              
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
        
    try
        
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

    catch MException
        
        if( ismac() )
            
            disp_string_framed(2, 'Mex files compilation failed' );

            fprintf(2, '\nAs we don''t include binaries for MAC, you will have to ask for help in order tu debug this error, here is the cause:\n');
            
            rethrow(MException);
            
        else
            
            if( bUseDesktop )
                
                disp_string_framed('[1,0.5,0]', 'Mex files compilation failed' );
                
                cprintf('[1,0.5,0]', '\nBut don''t worry, the kit includes several precompiled binaries for most architectures.\nTry the examples to check if yours is included. Down is the report error, use in case you want to ask for help in the forum:\n');
                % report just in case 

            else
                disp_string_framed(2, 'Mex files compilation failed' );
            end

            report = getReport(MException);
            fprintf(2, '\n%s\n', report);
            
        end
        
    end
         
    if(bOctave)
        cd(prev_folder);
    end
    
    if( bUseDesktop || bOctave )

        if( bMatlab )
        
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
          
          if( exist(tcxml_filename, 'file') && exist(tcxsd_filename, 'file'))
            status = copyfile(tcxml_filename, my_tcxml_filename, 'f');
            status = status & copyfile(tcxsd_filename, my_tcxsd_filename, 'f');
          else
            status = false;
          end
            
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

                          if( ~bStartFound && ~isempty(strfind( str_aux, '<!-- ecg-kit start -->')) )
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
                          if( ~bStartFound && ~isempty(strfind( str_aux, '<!-- ecg-kit start')) )
                              bStartFound = true;
                              break
                          elseif( isempty(start_idx) && ~isempty(strfind( str_aux, '<xsd:pattern value="-[A-Za-z0-9]+" />')) )
                              fseek(fid, 0, 'bof');
                              str_aux_pre = rowvec(fread(fid, last_binding, '*char'));
                              str_aux = rowvec(fgetl(fid));
                              str_aux_post = rowvec(fread(fid, inf, '*char'));
                              fseek(fid, 0, 'bof');
                              fwrite(fid, [str_aux_pre sprintf('<!-- ecg-kit start\n') str_aux sprintf('\necg-kit end -->\n	    <xsd:pattern value="[-_A-Za-z0-9]+" />\n') str_aux_post], 'char' );
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
                  disp_string_framed('[1,0.5,0]', 'Could not setup tab completion');
                  cprintf('[1,0.5,0]', 'Some error happened during the tab completion installation, but don''t worry, this will not diminish the ecg-kit functionality.\nDo it manually or ask for help in the forum.\n' );
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
                      disp_string_framed('[1,0.5,0]', 'Could not setup tab completion');
                      cprintf('[1,0.5,0]', 'Some error happened during the tab completion installation, but don''t worry, this will not diminish the ecg-kit functionality.\nDo it manually or ask for help in the forum.\n' );
                  end

              end
              
          end
        
        end
          
        bOk = savepath();
    
        if( bOk == 0 )

            fprintf(1, 'done !\n' );
            
            disp_string_framed('*Blue', sprintf('ecg-kit for Matlab %s', this_release_str) );
            
            str_aux = [ 'examples' filesep 'examples.m'];
            
            if( bMatlab )
            
              fprintf(1, '%s was correctly installed.\n\nYou can start reading the %s, or if you prefer, trying these %s.\nGo to the %s if you need help.\n', ...
                  '<a href = "matlab: web(''http://marianux.github.io/ecg-kit/'', ''-browser'' )">ecg-kit</a>', ...
                  '<a href = "matlab: web(''http://ecg-kit.readthedocs.org/en/master/'', ''-browser'' )">documentation</a>', ...
                  [ '<a href = "matlab: opentoline(' str_aux ',1)">examples</a>' ], ...
                  '<a href = "matlab: web(''https://groups.google.com/forum/#!forum/ecg-kit-users'', ''-browser'' )">forum</a>');
                  
             elseif( bOctave )
             
                fprintf(1, [ 'ecg-kit was correctly installed.\n\nYou can start reading the documentation, or if you prefer, trying these examples.\nGo to the forum if you need help.\n\n' ...
                    'web-page:      http://marianux.github.io/ecg-kit/\n', ...
                    'documentation: http://ecg-kit.readthedocs.org/en/master/\n', ...
                    'examples:      .\' filesep 'examples\' filesep 'examples.m\n', ...
                    'forum:         https://groups.google.com/forum/#!forum/ecg-kit-users\n' ]);
                  
             end
                  
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
                disp_string_framed(2, 'Path Not Saved');
                fprintf(2, 'Save path manually to install the ecg-kit permanently for future sessions.\n');
            end
            
        end

    end
    
end
