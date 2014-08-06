% InstallECGkit 
% -------------
% 
% Description:
% 
% Function for installing the Kit.
% 
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate  : 01/09/2012
% Last update: 22/07/2014

function InstallECGkit(bIgnoreUpdate)

    if( nargin < 1 || ~islogical(bIgnoreUpdate) )
        bIgnoreUpdate = false;
    end

    fprintf(1, 'Adding paths\n' );
    
    %path related constants.
    root_path = [fileparts(mfilename('fullpath')) filesep ];
    
    %path related constants.
    default_paths = { ...
                        [ root_path ';' ]; ...
                        [ root_path 'common' filesep ';' ]; ...
                        [ root_path 'common' filesep 'a2hbc' filesep ';' ]; ...
                        [ root_path 'common' filesep 'a2hbc'  filesep 'scripts' filesep ';' ]; ...
                        [ root_path 'common' filesep 'export_fig' filesep ';' ]; ...
                        [ root_path 'common' filesep 'wavedet' filesep ';' ]; ...
                        [ root_path 'common' filesep 'ppg' filesep ';' ]; ...
                        [ root_path 'common' filesep 'prtools' filesep ';' ]; ...
                        [ root_path 'common' filesep 'prtools_addins' filesep ';' ]; ...
                        [ root_path 'common' filesep 'bin' filesep ';' ]; ...
                        [ root_path 'common' filesep 'kur' filesep ';' ]; ...
                        [ root_path 'common' filesep 'LIBRA' filesep ';' ]; ...
                    };

    default_paths = char(default_paths)';
    default_paths = (default_paths(:))';
    addpath(default_paths);    
    
    ECGkitURL_str = 'https://code.google.com/p/ecg-kit/wiki/LatestVersion';
    % Version info.
    release_str = 'v0.1 beta - 22/07/2014';

    % URL check for updates
    if ~bIgnoreUpdate && ~verLessThan('matlab','8.0')
        fprintf(1, 'Checking for updates\n' );
        [strResponse, URLstatus] = urlread( ECGkitURL_str,'TimeOut',5);
        if( URLstatus == 1 )
            aux_idx = strfind(strResponse, '###');
            latest_release_str = strResponse(aux_idx(1)+3:aux_idx(2)-1);
            
            if( strcmpi(latest_release_str, release_str) )
                disp_string_framed('Blue', 'ECGkit is up to the date');
            else
                disp_string_framed('*Red', sprintf('ECGkit %s version available', latest_release_str) );
                fprintf(1, 'Consider updating %s.\n', '<a href = "https://code.google.com/p/ecg-kit/">here</a>' );
                
                rmpath(default_paths);
                return
            end
            
        end
    end

    fprintf(1, 'Compiling sources\n' );
    
    %Check compilation of source MEX files
    common_path = [ root_path 'common' filesep];
    source_files = dir([ common_path '*.c'] );
    lsource_files = length(source_files);
    for ii = 1:lsource_files
        [~, source_file_name] = fileparts( source_files(ii).name);
        mex_file = dir([common_path  source_file_name '.' mexext ]);
        if( isempty(mex_file) || mex_file.datenum <= source_files(ii).datenum  )
            eval(['mex -outdir ''' common_path ''' ''' [common_path source_files(ii).name] '''']);
        end
    end
    
    bUseDesktop = usejava('desktop');
    
    if( bUseDesktop )

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
    
        if( bOk == 0 )

            fprintf(1, 'done !\n' );
            
            disp_string_framed('*Blue', sprintf('ECGkit for Matlab %s', release_str) );
            
            fprintf(1, 'Kit was correctly installed.\n\nYou can start reading the %s, or if you prefer, trying these %s.\nGo to the %s if you need help.\n', ...
                '<a href = "matlab: web(''https://code.google.com/p/ecg-kit/w/list'', ''-browser'' )">documentation</a>', ...
                '<a href = "matlab: opentoline(examples.m,1)">examples</a>', ...
                '<a href = "matlab: web(''https://groups.google.com/d/forum/ecgkit-users'', ''-browser'' )">forum</a>');
        else
            
            disp_string_framed(2, 'Path Not Saved')
            fprintf(2, 'Save path manually to install the ECG-Kit permanently.\n');
            
        end

        
    end
    
end
