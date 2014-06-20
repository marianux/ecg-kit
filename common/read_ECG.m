function [ECG heasig ann recording_format] = read_ECG(recording_name, ECG_start_idx, ECG_end_idx, recording_format)

ECG = [];
heasig = [];
ann = [];
cKnownFormats = {'MIT' 'ISHNE', 'HES', 'MAT', 'Mortara'};

if( nargin < 4 || isempty(recording_format) || ~any(strcmpi( recording_format, cKnownFormats)) )
    %Try guessing the ECG format
    recording_format = ECGformat(recording_name);
    if( isempty(recording_format) )
        error('read_ECG:UnknownDataFormat', 'Unknown data format.')
    end
end

if( nargin < 2 || ~isnumeric(ECG_start_idx) )
    ECG_start_idx = [];
end

if( nargin < 3 || ~isnumeric(ECG_end_idx))
    ECG_end_idx = [];
end

if( strcmp(recording_format, 'ISHNE') )
    if( nargout > 1 )
        [ECG heasig ann] = read_ishne(recording_name, ECG_start_idx, ECG_end_idx );
    else
        ECG = read_ishne(recording_name, ECG_start_idx, ECG_end_idx );
    end        
elseif( strcmp(recording_format, 'MAT') )
    
    aux_load = load(recording_name);

    if( nargout > 1 )
        
        heasig = aux_load.header;
        ECG = aux_load.signal;
        
        if( isfield(aux_load, 'ann') )
            ann = aux_load.ann;            
        end
    else
        ECG = aux_load.ECG;
    end
    clear aux_load
    
elseif( strcmp(recording_format, 'Mortara') )
    
    [ECG heasig ] = read_Mortara(recording_name, ECG_start_idx, ECG_end_idx );
    
elseif( strcmp(recording_format, 'AHA') )
    if( nargout > 1 )
        [ECG heasig ann] = read_AHA_format(recording_name, ECG_start_idx, ECG_end_idx );
    else
        ECG = read_AHA_format(recording_name, ECG_start_idx, ECG_end_idx );
    end
elseif( strcmp(recording_format, 'HES') )
    if( nargout > 1 )
        [ECG heasig ann] = read_HES_format(recording_name, ECG_start_idx, ECG_end_idx );
    else
        ECG = read_HES_format(recording_name, ECG_start_idx, ECG_end_idx );
    end
elseif( strcmp(recording_format, 'MIT') )
    %% MIT various data formats

    if( nargout > 2 )
        strAnnExtension = {'ari' 'atr' 'ecg' 'cardmean'};
        bAnnotationFound = false;
        for ii = 1:length(strAnnExtension)
            annFileName = [recording_name(1:end-3) strAnnExtension{ii}];
            if( exist(annFileName, 'file') )
                bAnnotationFound = true;
                break
            end
        end
        if(~bAnnotationFound)
            error( 'read_ECG:AnnNotFound', ['Can�t find annotations for file : ' rowvec(recording_name') ] );
        end

        ann = readannot(annFileName);
    end

    [ recording_path rec_name ] = fileparts(recording_name);
    recording_path = [recording_path filesep];
    
    heasig = readheader([recording_name(1:end-3) 'hea']);    
    
    if( ~isfield(heasig, 'fname') ) 
        heasig.fname = repmat(rec_name, heasig.nsig, 1);
    end
    
    recording_files = unique(cellstr(heasig.fname));

    % Sometime I filtered only ECG signals right here.
%     ECG_idx = get_ECG_idx_from_header(heasig);
    ECG_idx = 1:heasig.nsig;
    
    if( isempty( ECG_start_idx ) )
        ECG_start_idx = 1;
    else
        ECG_start_idx = max(1,ECG_start_idx);
    end

    if( isempty( ECG_end_idx ) )
        %Intento la lectura total por defecto
        ECG_end_idx = heasig.nsamp;
    else
        ECG_end_idx = min(heasig.nsamp, ECG_end_idx);
    end
    
    if( ECG_end_idx == 0 )
        ECG_end_idx = realmax;
    end
    
    for ii = 1:length(recording_files)
        sig_idx = find(strcmpi(recording_files(ii), cellstr(heasig.fname)));
        % this should respect the stride of the data
        nsig_present = length(sig_idx);
        sig2_idx = intersect(sig_idx, ECG_idx);
        nsig2read = length(sig2_idx);
        fmt = heasig.fmt(sig2_idx(1));
        ECG(:,sig2_idx) = read_MIT_ecg( [recording_path heasig.fname(sig2_idx(1),:)], ECG_start_idx, ECG_end_idx, nsig2read, nsig_present, fmt, heasig);
    end    
%     ECG = int16(ECG);
    
end

heasig.ECG_format = recording_format;

function ECG = read_MIT_ecg(recording_name, ECG_start_idx, ECG_end_idx, nsig2read, nsig_present, fmt, heasig)
    
    if( fmt == 16 || fmt == 61 || isnan(fmt) )

        if( fmt == 16 || isnan(fmt) )
            byte_ordering = 'ieee-le';
        else
            byte_ordering = 'ieee-be';
        end
        
        ECG_size = ECG_end_idx - (ECG_start_idx-1);
        
        fidECG = fopen(recording_name, 'r');
        try
            fseek(fidECG, ((ECG_start_idx-1)*nsig_present)*2, 'bof');
            ECG = fread(fidECG, [nsig2read ECG_size ], '*int16', (nsig_present-nsig2read)*2 , byte_ordering)';
            fclose(fidECG);
        catch ME
            fclose(fidECG);
            rethrow(ME);
        end

    elseif(fmt == 212)

        ECG = rdsign212(recording_name, heasig.nsig, ECG_start_idx, ECG_end_idx);

    elseif(fmt == 310)

        ECG = read_310_format(recording_name, ECG_start_idx, ECG_end_idx, heasig  );

    elseif(fmt == 311)

        warning('read_ECG:UntestedRegion', 'Untested!! danger ...')
        ECG = read_311_format(recording_name, ECG_start_idx, ECG_end_idx, heasig  );

    else
        error('read_ECG:UnknownDataFormat', 'Unknown data format.')
    end
    
    
