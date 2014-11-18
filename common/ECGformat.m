function recording_format = ECGformat(recording_name)
    
    recording_format = [];

    cKnownFormats = {'MIT' 'ISHNE', 'AHA', 'HES', 'MAT', 'Mortara'};
    lcKnownFormats = length(cKnownFormats);

    [rec_filepath, rec_filename, rec_fileExt] = fileparts(recording_name);

    aux_idx = 1:lcKnownFormats;

    % Mortara format is a folder with predefined files inside named
    % HourXXRawData.bin
    if( isdir(recording_name) )
        rec_path = recording_name;
        hours_files = dir([rec_path filesep 'Hour*.bin' ]);
        if( isempty(hours_files) )
            % empty folder, not Mortara
            return
        else
            recording_format = 'Mortara';
            return
        end
    else
        [~, rec_name] = fileparts(recording_name);
        if( length(rec_name) > 7 && strcmpi( rec_name(1:4), 'Hour') && strcmpi( rec_name(end-6:end), 'RawData') )
            recording_format = 'Mortara';
            return
        end
    end
    
    
    % first guess based on typical file extensions
    if( strcmpi(rec_fileExt, '.dat' ) )
        % MIT first
        fg_idx = find(strcmp(cKnownFormats, 'MIT'));
        aux_idx = [fg_idx aux_idx(aux_idx~=fg_idx)];
    elseif( strcmpi(rec_fileExt, '.ecg' ) )
        % ISHNE and AHA first
        fg_idx = find(strcmp(cKnownFormats, 'ISHNE'));
        aux_idx = [fg_idx aux_idx(aux_idx~=fg_idx)];
        fg_idx = find(strcmp(cKnownFormats, 'AHA'));
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

        elseif( strcmp(cKnownFormats{ii}, 'AHA') )
            
            if(isAHAformat(recording_name))
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
            
            catch ME
                
            end

        end

    end
    
    