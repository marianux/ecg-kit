%% (Internal) Init environment variables for using WFDB toolbox
%   
%   WFDB_bin_path =  init_WFDB_library(tmp_path_local)
% 
% Arguments:
% 
%      + tmp_path_local: a path for temporary data
% 
% Output:
% 
%      + WFDB_bin_path: The path to the WFDB command, depending of the architecture used.
% 
% Example:
% 
% See also get_ECG_idx_from_header, ECGwrapper
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function WFDB_bin_path =  init_WFDB_library(tmp_path_local)

    WFDB_paths = { ...
                    ['wfdb' filesep 'linux-amd64' filesep ] ...
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

     elseif( ismac() )            
        path_OS_var = 'PATH';
        libpath_OS_var = 'DYLD_LIBRARY_PATH';
        path_sep = ':';
        bin_path  = WFDB_paths{WFDB_MAC_idx};
        lib_path  = WFDB_paths{WFDB_MAC_idx};

    elseif( isunix()  )            
        path_OS_var = 'PATH';
        libpath_OS_var = 'LD_LIBRARY_PATH';
        path_sep = ':';
        bin_path  = [ WFDB_paths{WFDB_UNIX_idx} 'bin'];
        lib_path  = [ WFDB_paths{WFDB_UNIX_idx} 'lib64'];

    else
        warning('Unknown system architecture. WFDB library QRS detectors probably won''t work.')
        path_OS_var = 'PATH';
    end

    % this file is located in common folder.
    common_path = fileparts(mfilename('fullpath'));
    WFDB_bin_path = [common_path filesep bin_path filesep];         
    
    aux_str = tmp_path_local;
    if( aux_str(end) == filesep )
        aux_str = aux_str(1:end-1);
    end
    wfdb_val = getenv('WFDB');
    if( isempty(strfind(wfdb_val, aux_str ) ))
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
        
        % WFDB paths required to locate shared libraries
        if( isunix() )
            WFDB_lib_path = [common_path filesep lib_path filesep ];
            this_path = getenv(libpath_OS_var);
            setenv(libpath_OS_var, [this_path path_sep WFDB_lib_path ]);
        end
        
        WFDB_initiated = true;
        
    end
    

    