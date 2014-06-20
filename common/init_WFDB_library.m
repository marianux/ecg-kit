function WFDB_bin_path =  init_WFDB_library(tmp_path_local)

    WFDB_paths = { ...
                    ['wfdb' filesep 'linux-amd64' filesep 'bin' ] ...
                    ['wfdb' filesep 'macosx-x86_64' filesep 'bin' ] ...
                    ['wfdb' filesep 'windows-amd64' filesep 'bin' ] ...
                    };
    WFDB_UNIX_idx = 1;
    WFDB_MAC_idx = 2;
    WFDB_WIN_idx = 3;

    persistent WFDB_initiated

    if( ispc() )
        path_OS_var = 'Path';
        path_sep = ';';
        bin_path  = WFDB_paths{WFDB_WIN_idx};

    elseif( isunix()  )            
        path_OS_var = 'PATH';
        path_sep = ':';
        bin_path  = WFDB_paths{WFDB_UNIX_idx};

    elseif( ismac() )            
        path_OS_var = 'PATH';
        path_sep = ':';
        bin_path  = WFDB_paths{WFDB_MAC_idx};

    else
        warning('Unknown system architecture. WFDB library QRS detectors probably won''t work.')
        path_OS_var = 'PATH';
    end

    % this file is located in common folder.
    common_path = fileparts(mfilename('fullpath'));
    WFDB_bin_path = [common_path filesep bin_path filesep];         
    
    wfdb_val = getenv('WFDB');
    if( isempty(strfind(wfdb_val, tmp_path_local ) ))
        aux_str = tmp_path_local;
        if( aux_str(end) == filesep )
            aux_str = aux_str(1:end-1);
        end
        if( isempty(wfdb_val))
            wfdb_val = aux_str;
        else
            wfdb_val = [wfdb_val path_sep aux_str];
        end
    end
    setenv('WFDB', wfdb_val);     
    
    if( isempty(WFDB_initiated) || ~WFDB_initiated )

        % WFDB paths required to locate recordings and binaries
        this_path = getenv(path_OS_var);
        if( isempty(strfind(this_path, WFDB_bin_path ) ))
            setenv(path_OS_var, [this_path path_sep WFDB_bin_path]);
        end
        
        WFDB_initiated = true;
        
    end
    

    