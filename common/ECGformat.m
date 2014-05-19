function recording_format = ECGformat(recording_name)
    
    recording_format = [];

    cKnownFormats = {'MIT' 'ISHNE', 'HES', 'MAT'};
    lcKnownFormats = length(cKnownFormats);

    [rec_filepath, rec_filename, rec_fileExt] = fileparts(recording_name);

    aux_idx = 1:lcKnownFormats;
    
    % first guess based on typical file extensions
    if( strcmpi(rec_fileExt, '.dat' ) )
        % MIT first
        fg_idx = find(strcmp(cKnownFormats, 'MIT'));
        aux_idx = [fg_idx aux_idx(aux_idx~=fg_idx)];
    elseif( strcmpi(rec_fileExt, '.ecg' ) )
        % ISHNE and AHA first
        fg_idx = find(strcmp(cKnownFormats, 'ISHNE'));
        aux_idx = [fg_idx aux_idx(aux_idx~=fg_idx)];
    elseif( strcmpi(rec_fileExt, '.mat' ) )
        % Matlab first
        fg_idx = find(strcmp(cKnownFormats, 'MAT'));
        aux_idx = [fg_idx aux_idx(aux_idx~=fg_idx)];
    elseif( strcmpi(rec_fileExt, '.hes' ) )
        % Matlab first
        fg_idx = find(strcmp(cKnownFormats, 'HES'));
        aux_idx = [fg_idx aux_idx(aux_idx~=fg_idx)];
    end
    
    for ii = aux_idx

        if( strcmp(cKnownFormats{ii}, 'MIT') )
            
            header_fn = [rec_filepath filesep rec_filename '.hea'];
            if( exist( header_fn, 'file' ) )
                aux_header = readheader(header_fn);
                if( strcmp(aux_header.recname, rec_filename) )
                    recording_format = cKnownFormats{ii};
                    return
                end
            end

        elseif( strcmp(cKnownFormats{ii}, 'ISHNE') )
            
            if(isISHNEformat(recording_name))
                recording_format = cKnownFormats{ii};
                return
            end

        elseif( strcmp(cKnownFormats{ii}, 'HES') )

            if(isHESformat(recording_name))
                recording_format = cKnownFormats{ii};
                return
            end
            
        elseif( strcmp(cKnownFormats{ii}, 'MAT') )

            try
                
                aux_load = load(recording_name, 'header');

                recording_format = cKnownFormats{ii};
                return
                
%                 if( strcmp(aux_load.header.recname, rec_filename) )
%                     recording_format = cKnownFormats{ii};
%                     return
%                 end
            
            end

        end

    end
    
    