function init_ghostscript()

    persistent GS_initiated

    if( isempty(GS_initiated) || ~GS_initiated )
    
        if( ispc() )
            path_OS_var = 'Path';
            path_sep = ';';

        elseif( isunix()  )            
            path_OS_var = 'PATH';
            path_sep = ':';

        elseif( ismac() )            
            path_OS_var = 'PATH';
            path_sep = ':';

        else
            warning('Unknown system architecture.')
            path_OS_var = 'PATH';
        end

        % this file is located in common folder.
        common_path = fileparts(mfilename('fullpath'));
        GS_bin_path = [common_path filesep 'bin' ];         


        % WFDB paths required to locate recordings and binaries
        this_path = getenv(path_OS_var);
        if( isempty(strfind(this_path, GS_bin_path ) ))
            setenv(path_OS_var, [this_path path_sep GS_bin_path]);
        end
        
        GS_initiated = true;
        
    end
    

    