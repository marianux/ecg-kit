function detect_QRS_complexes( rec_name_full_path, varargin )

cQRSdetectors = {'all-detectors' 'wavedet' 'pantom' 'aristotle' 'gqrs' 'sqrs' 'wqrs'};

AAMI_anntypes =  ('NSVFQ')';

% cantidad minima de anotaciones para intentar crear nuevas
min_annotations = 1;

if( ~exist(rec_name_full_path, 'file') )
    error([ 'File not found ' rec_name_full_path ])
end

p = inputParser;   % Create instance of inputParser class.
p.addParamValue('detectors', cQRSdetectors{1}, @(x)( (ischar(x) || iscellstr(x)) && ~isempty(intersect(cQRSdetectors, cellstr(x))) ) );
p.addParamValue('gqrs_config', 'gqrs.conf', @(x)( ischar(x) ));

try
    p.parse( varargin{:} );
catch MyError
    rethrow(MyError);
end

detector = p.Results.detectors;
gqrs_config_filename = p.Results.gqrs_config;


if( ispc() ) 

    folder_paths = { ...
                    'D:\Mariano\misc\ECGkit\common' ...
                    };
    
    tmp_path = 'd:\tmp\';
    tmp_path_local = tmp_path;
    
else
    
    folder_paths = { ...
                    '/home/bio/mllamedo/ECGKit/common' ...
                    '/home/bio/mllamedo/ECGKit/common/wavedet' ...
                    '/home/bio/mllamedo/ECGKit/common/prtools_addins' ...
                    '/home/bio/mllamedo/ECGKit/common/prtools' ...
                    };
    
    tmp_path = '/extra/scratch/bio/mllamedo/tmp/detector';
    tmp_path_local = '/scratch/mllamedo/';

end

if(~exist(tmp_path, 'dir'))
    mkdir(tmp_path);
end

if(~exist(tmp_path_local, 'dir'))
    mkdir(tmp_path_local);
end

if( strcmpi('all-detectors', detector) )
    detectors2do = cQRSdetectors(2:end);
else
    if( ischar(detector) )
        detectors2do = {detector};
    else
        detectors2do = detector;
    end
end

added_paths = addpath_if_not_added(folder_paths);

%%     MIT conversion

if(exist(gqrs_config_filename, 'file'))
    gqrs_config_filename = which(gqrs_config_filename);
end
    
[ rec_path, rec_name ]= fileparts(rec_name_full_path);

rec_name(rec_name == ' ') = '_';

rec_path = [rec_path filesep];

file_output = [rec_path rec_name '.mat'];
file_output_local = [tmp_path_local rec_name '.mat'];
    
if( exist([rec_path rec_name '.hea'], 'file') )

    if( exist(file_output, 'file') )
        aux_struct = load( file_output );
    else
        
        [ ECG, header] = read_ECG([rec_path rec_name '.hea']);
        
        aux_struct.signal = ECG;
        aux_struct.header = header;

        save( file_output_local, '-append', '-struct', 'aux_struct');
        
    end

else

    [ ECG, header, ann, rec_format ] = read_ECG(rec_name_full_path);

    if( ~isempty(ann) ) 
        ann = AnnotationFilterConvert(ann, rec_format, 'AAMI');
        ann.anntyp = AAMI_anntypes(ann.anntyp);
        
        file_name = [tmp_path_local filesep rec_name '.atr'];
        pos2mit(file_name, ann.time, ann.anntyp);
        
    end

    % porque en algun momento le puse un entero_ como prefijo
    % para que siempre sea el mismo orden
    header.recname = rec_name;

    file_name = [tmp_path_local rec_name '.dat'];
    fidECG = fopen(file_name, 'w');
    try
        fwrite(fidECG, ECG', 'int16', 0 );
        fclose(fidECG);
    catch MEE
        fclose(fidECG);
        rethrow(MEE);
    end

    writeheader(tmp_path_local, header);   
    
    if( strcmpi('MAT',rec_format) )
        
        aux_struct = load(file_output);
        
    else
        aux_struct.signal = ECG;
        aux_struct.header = header;

        if( ~isempty(ann) )
            aux_struct.ann = ann;
        end

        save( file_output_local, '-append', '-struct', 'aux_struct');
    end
    
    clear ECG header ann
    
end

% Período refractario. Mín RR posible. Como son registros en reposo 0.4
% estaría bien.
wavedet_config = [];
%     wavedet_config.setup.wavedet.refrper = 0.4; % seconds


% es necesario fijar esta variable de entorno para que los programas de la
% WFDB sepan donde buscar los registros.
wfdb_val = getenv('WFDB');

if( isempty(strfind(wfdb_val, rec_path ) ))
    aux_str = rec_path;
    if( isempty(wfdb_val))
        wfdb_val = aux_str(1:end-1);
    else
        wfdb_val = [wfdb_val ':' aux_str(1:end-1)];
    end
end

setenv('WFDB', wfdb_val);    

cant_QRSdetectors = length(detectors2do);

for ii = 1:cant_QRSdetectors

    this_detector = detectors2do{ii};

    fprintf(1, [ '\nProcessing ' this_detector '\n\n' ] );
    
%% WORK LOAD
    
    %% perform QRS detection
    
    switch( this_detector )

        case 'wavedet'
        %% Wavedet delineation

            try
                
                [position_multilead, positions_single_lead] = wavedet_interface(aux_struct.signal, aux_struct.header, [], [], wavedet_config);

                for jj = 1:aux_struct.header.nsig
                    aux_struct.(['wavedet_' num2str(jj)]) = positions_single_lead(jj);
                end
                aux_struct.wavedet_multilead = position_multilead;

                save( file_output_local, '-struct', 'aux_struct');
                
            catch aux_ME

                strAux = sprintf('Wavedet failed in recording %s\n', rec_name);
                strInvocation = GetFunctionInvocation(mfilename, {rec_name_full_path, detector} );
                strAux = sprintf('%s\nMatlab function invocation:\n%s\n', strAux, strInvocation);

                report = getReport(aux_ME);
                fprintf(2, '%s\nError report:\n%s', strAux, report);

            end
            
            
        case 'pantom'
            %% Pan-Tompkins detector
            
            for jj = 1:aux_struct.header.nsig

                try

                    peaks = PeakDetection2(aux_struct.signal(:,jj), aux_struct.header.freq);

                    aux_struct.(['pantom_' num2str(jj)]).time = peaks;

                catch aux_ME

                    strAux = sprintf('Pan & Tompkins failed in recording %s lead %d\n', rec_name, jj);
                    strInvocation = GetFunctionInvocation(mfilename, {rec_name_full_path, detector} );
                    strAux = sprintf('%s\nMatlab function invocation:\n%s\n', strAux, strInvocation);

                    report = getReport(aux_ME);
                    fprintf(2, '%s\nError report:\n%s', strAux, report);
                    
                end
                

            end

            save( file_output_local, '-struct', 'aux_struct');
            
        case { 'aristotle' 'gqrs' 'sqrs' 'wqrs' }
            %% WFDB_comp_interface

            for jj = 0:(aux_struct.header.nsig-1)
                
                file_name_orig = [];
                
                try

                    if( any(strcmpi( {'sqrs' 'wqrs'} , this_detector ) ) )
                        % wqrs tiene una interface diferente
                        system(['cd ' tmp_path_local ';' this_detector ' -r ' rec_name ' -s ' num2str(jj)]);
                        
                    elseif( strcmpi( 'gqrs' , this_detector  ) && exist(gqrs_config_filename, 'file') )
                        
                        system(['cd ' tmp_path_local ';gqrs -c ' gqrs_config_filename ' -r ' rec_name ' -s ' num2str(jj)]);
                        system(['cd ' tmp_path_local ';gqpost -c ' gqrs_config_filename ' -r ' rec_name ' -o ' this_detector num2str(jj) ]);
                        delete([tmp_path_local rec_name '.qrs' ]);
                        
                    else
                        system(['cd ' tmp_path_local ';' this_detector ' -r ' rec_name ' -o ' this_detector num2str(jj) ' -s ' num2str(jj)]);
                    end

                    if( strcmpi( 'sqrs' , this_detector ) )
                        file_name_orig =  [tmp_path_local rec_name '.qrs' ];
                    else
                        file_name_orig =  [tmp_path_local rec_name '.' this_detector num2str(jj) ];
                    end

                catch aux_ME

                    strAux = sprintf( '%s failed in recording %s lead %d\n', this_detector, rec_name, jj);
                    strInvocation = GetFunctionInvocation(mfilename, {rec_name_full_path, detector} );
                    strAux = sprintf('%s\nMatlab function invocation:\n%s\n', strAux, strInvocation);

                    report = getReport(aux_ME);
                    fprintf(2, '%s\nError report:\n%s', strAux, report);

                end

                if( ~isempty(file_name_orig) )
                    % reference comparison
                    anns_test = [];                    
                    try
                        anns_test = readannot(file_name_orig);
                        
                        anns_test = AnnotationFilterConvert(anns_test, 'MIT', 'AAMI');

                        aux_struct.([this_detector '_' num2str(jj+1)]) = anns_test;
                        
                    catch aux_ME
                        if( ~strcmpi(aux_ME.identifier, 'MATLAB:nomem') )
                            strAux = sprintf( '%s failed in recording %s lead %d\n', this_detector, rec_name, jj);
                            strInvocation = GetFunctionInvocation(mfilename, {rec_name_full_path, detector} );
                            strAux = sprintf('%s\nMatlab function invocation:\n%s\n', strAux, strInvocation);

                            report = getReport(aux_ME);
                            fprintf(2, '%s\nError report:\n%s', strAux, report);
                        end
                    end

                end
                
            end
            
            save( file_output_local, '-struct', 'aux_struct');
            
    end


end

% attemp to build a better detection from single-lead detections.

aux_struct = calculate_artificial_QRS_detections(aux_struct);

[AnnNames, all_annotations] = getAnnNames(aux_struct);

save( file_output_local, '-struct', 'aux_struct');

% en caso de trabajar sobre el directorio temporal, no lo muevo. Se
% encargarian de afuera.
if( ~strcmpi(file_output_local, file_output) )
    movefile( file_output_local, file_output, 'f');
end

if( ~isempty(added_paths) )
    cellfun( @rmpath, added_paths);
end

