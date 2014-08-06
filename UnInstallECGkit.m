function UnInstallECGkit()

    release_str = 'v0.1 beta - 22/07/2014';

    bUseDesktop = usejava('desktop');

    %path related constants.
    root_path = [fileparts(mfilename('fullpath')) filesep ];

    if( bUseDesktop )
    
        disp_string_framed('*Blue', sprintf('ECGkit for Matlab %s', release_str) );

        fprintf(1, 'Removing paths. ' );

    end
    
    %path related constants.
    default_paths = { ...
                        [root_path ';' ]; ...
                        [root_path 'common' filesep ';' ]; ...
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
    rmpath(default_paths);

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

            fprintf(1, 'Thanks for trying the kit.\n\nNow you can safely delete %s\n', root_path );
        
        else
           
            disp_string_framed(2, 'Path Not Saved')
            fprintf(2, 'Save path manually to uninstall the ECG-Kit permanently. Deletes all occurrences of %s in your path.\n', root_path);
            
        end
        
    end
    
end
