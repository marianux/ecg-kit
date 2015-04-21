%% Reads an ECG recording
% This function is the main I/O interface with ECG recordings. It can read
% several formats, such as 'MIT' 'AHA' 'ISHNE', 'HES', 'MAT' and 'Mortara'.
% Also has the feature of "auto-detect" file format.
% 
% Example
% 
%   [ECG heasig ann recording_format] = read_ECG(recording_name, ECG_start_idx, ECG_end_idx, recording_format)
% 
% Arguments: (specified as ECGwrapper('arg_name1', arg_val1, ... , 'arg_nameN', arg_valN) )
% 
%           + recording_name : (char) the full filename of the ECG
%                recording. 
% 
%           + ECG_start_idx : (number) the sample to start reading. Empty
%                             or 1 to start at the begining.
% 
%           + ECG_end_idx : (number) the sample to end reading. Empty to
%                           read all.
% 
%           + recording_format: (char) recording format of "recording_name"
%                 file. You can check the following links regarding file
%                 formats: 
% 
%                   - AHA: http://physionet.org/physiotools/old/dbpg/dbu_84.htm
%                   - ISHNE: http://thew-project.org/THEWFileFormat.html
%                   - MIT: http://physionet.org/physiotools/wag/signal-5.htm
% 
%                 HES and mortara are propietary formats and no
%                 information is distributed with the kit. You may need
%                 SuperECG from mortara in order con extract the binary
%                 files required to read this kind of files. MAT format
%                 requieres a typical "*.mat" matlab file with the
%                 following variables inside:
% 
%                   - signal: the signal matrix of size [ header.nsamp header.nsig ]
%                   - header: A struct with all signal properties. See the
%                       example recording included with the kit
%                       \recordings\example_recording.mat for an example of
%                       this struct format, and Physionet documentation for
%                       further details:
%                       http://physionet.org/physiotools/wag/header-5.htm  
%                   - ann: Annotations struct included with the signal in MIT
%                       format as the one returned by "readannot" function.
%                       The "time" field must be present, being this a
%                       numeric vector with the sample indexes of the
%                       annotations. For the "type" and "subtype" fields, see
%                       the MIT format documentation in Physionet:
%                    http://www.physionet.org/physiobank/annotations.shtml
% 
% See also ECGwrapper, ECGformat, read_AHA_format, read_ishne_format, read_HES_format, read_Mortara_format
% 
% Author: Mariano Llamedo Soria
% <matlab:web('mailto:llamedom@electron.frba.utn.edu.ar','-browser') (email)> 
% Version: 0.1 beta
% Birthdate: 5/01/2014
% Last update: 19/11/2014
% Copyright 2008-2015
% 
function [ECG heasig ann recording_format end_sample] = read_ECG(recording_name, ECG_start_idx, ECG_end_idx, recording_format)

ECG = [];
heasig = [];
ann = [];
cKnownFormats = {'MIT' 'AHA' 'ISHNE', 'HES', 'MAT', 'Mortara'};

matformat_definitions();

end_sample = [];

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
        [ECG heasig ann end_sample] = read_ishne_format(recording_name, ECG_start_idx, ECG_end_idx );
    else
        ECG = read_ishne_format(recording_name, ECG_start_idx, ECG_end_idx );
    end        
elseif( strcmp(recording_format, 'MAT') )
    
    aux_load = load(recording_name);

    fnames = fieldnames(aux_load);
    signal_name = intersect( fnames, cMatSignalNames);
    header_name = intersect( fnames, cMatSignalHeaderNames);
    ann_name = intersect( fnames, cMatSignalAnnNames);
    
    if( nargout > 1 )
        
        if( isempty(header_name) )
            error( 'read_ECG:HeaderNotFound', ['Can''t find annotations for file : ' rowvec(recording_name') ] );
        else            
            heasig = aux_load.(header_name{1});
        end
        
        if( ~isempty(ann_name) )
            ann = aux_load.(ann_name{1});            
        end
        
    end

    if( ~isempty(signal_name) )
        ECG = aux_load.(signal_name{1});
    end
    
    ECG = ECG(ECG_start_idx:ECG_end_idx,:);
    
    clear aux_load
    
elseif( strcmp(recording_format, 'Mortara') )
    
    [ECG heasig end_sample] = read_Mortara_format(recording_name, ECG_start_idx, ECG_end_idx );
    
elseif( strcmp(recording_format, 'AHA') )
    if( nargout > 1 )
        [ECG heasig ann end_sample] = read_AHA_format(recording_name, ECG_start_idx, ECG_end_idx );
    else
        ECG = read_AHA_format(recording_name, ECG_start_idx, ECG_end_idx );
    end
elseif( strcmp(recording_format, 'HES') )
    if( nargout > 1 )
        [ECG heasig ann end_sample] = read_HES_format(recording_name, ECG_start_idx, ECG_end_idx );
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
            error( 'read_ECG:AnnNotFound', ['Can''t find annotations for file : ' rowvec(recording_name') ] );
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
    
    samples2read = ECG_end_idx - ECG_start_idx + 1;
    %No leer bloques mas grandes de 200 megabytes
    MaxIOread = 200; %megabytes
    
    if( (samples2read*heasig.nsig*2) > (MaxIOread * 1024^2) )
        samples2read = (MaxIOread * 1024^2) / heasig.nsig / 2;
        ECG_end_idx = samples2read + ECG_start_idx - 1;
        warning(['Limited to read ' num2str(MaxIOread) ' Mb.'])
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
    
    if( exist(recording_name, 'file') )

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
    
    else
        
        error('read_ECG:FileNotFound', disp_string_framed(0, sprintf( 'Could not find %s', recording_name ) ) );
        
    end
    
