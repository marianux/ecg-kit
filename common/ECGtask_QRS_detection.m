classdef ECGtask_QRS_detection < ECGtask

% ECGtask for ECGwrapper (for Matlab)
% ---------------------------------
% 
% Description:
% 
% Abstract class for defining ECGtask interface
% 
% Adding user-defined QRS detectors:
% A QRS detector that has the following interface can be added to the task:
% 
%     [positions_single_lead, position_multilead] = your_QRS_detector( ECG_matrix, ECG_header, progress_handle, payload_in);

% where the arguments are:
%    + ECG_matrix, is a matrix size [ECG_header.nsamp ECG_header.nsig]
%    + ECG_header, is a struct with info about the ECG signal, such as:
%         .freq, the sampling frequency
%         .desc, description about the signals.
%    + progress_handle, is a handle to a waitbar object, that can be used
%          to track the progress within your function. See the
%          documentation about this class in this kit.
%    + payload_in, is a user data variable allowed to be sent each call to
%          your function. It is sent, via the payload property of this
%          class, for example:
% 
%         this_ECG_wrappers.ECGtaskHandle.payload = your_variable;
%         this_ECG_wrappers.ECGtaskHandle.payload = {your_var1 your_var2};
%         this_ECG_wrappers.ECGtaskHandle.payload = load(cached_filenames);
% 
% the output of your function must be:
%    + positions_single_lead, a cell array size ECG_header.nsig with the
%          QRS sample locations found in each lead.
%    + position_multilead, a numeric vector with the QRS locations
%          calculated using multilead rules.
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate  : 18/2/2013
% Last update: 18/2/2013
       
    properties(GetAccess = public, Constant)
        name = 'QRS_detection';
        target_units = 'ADCu';
        doPayload = true;
    end

    properties( GetAccess = public, SetAccess = private)
        % if user = memory;
        % memory_constant is the fraction respect to user.MaxPossibleArrayBytes
        % which determines the maximum input data size.
        memory_constant = 0.3;
        
        started = false;
        
    end
    
    properties( Access = private, Constant)
        
        cQRSdetectors = {'all-detectors' 'wavedet' 'pantom' 'gqrs' };
        
    end
    
    properties( Access = private )
        
        detectors2do
        bWFDBdetectors
        tmp_path_local
        WFDB_bin_path
        WFDB_cmd_prefix_str 
        lead_idx = [];
        
    end
    
    properties
        progress_handle
        user_string = '';
        tmp_path
        detectors = 'all-detectors';
        only_ECG_leads = false;
        gqrs_config_filename = [];
        wavedet_config
        payload
        CalculateArtificialDetections = true;
    end
    
    methods
           
        function obj = ECGtask_QRS_detection(obj)

            obj.wavedet_config.setup.wavedet.QRS_detection_only = true;
            
        end
        
        function Start(obj, ECG_header, ECG_annotations)

            if( obj.only_ECG_leads )
                obj.lead_idx = get_ECG_idx_from_header(ECG_header);
            else
%                 'all-leads'
                obj.lead_idx = 1:ECG_header.nsig;
            end
            
            if( isempty(obj.lead_idx) )
                return
            end
            
            if( any(strcmpi('all-detectors', obj.detectors)) )
                obj.detectors2do = obj.cQRSdetectors(2:end);
            else
                if( ischar(obj.detectors) )
                    obj.detectors2do = cellstr(obj.detectors);
                else
                    obj.detectors2do = obj.detectors;
                end
            end

            obj.bWFDBdetectors = any(strcmpi('gqrs', obj.detectors2do));
            
            % local path required to avoid network bottlenecks in distributed filesystems 
            if( isunix() && exist('/scratch/', 'dir') )
                str_username = getenv('USER');
                obj.tmp_path_local = ['/scratch/' str_username filesep];
                if( ~exist(obj.tmp_path_local, 'dir') )
                    if(~mkdir(obj.tmp_path_local))
                        obj.tmp_path_local = '/scratch/';
                    end
                end
                obj.tmp_path = tempdir;
            else
                if( isempty(obj.tmp_path) )
                    obj.tmp_path = tempdir;
                end
                obj.tmp_path_local = obj.tmp_path;
            end
            
            obj.WFDB_bin_path = init_WFDB_library(obj.tmp_path_local);
            
            obj.WFDB_cmd_prefix_str = WFDB_command_prefix(obj.tmp_path_local);
            
            obj.started = true;
        end
        
        function payload_out = Process(obj, ECG, ECG_start_offset, ECG_sample_start_end_idx, ECG_header, ECG_annotations, ECG_annotations_start_end_idx  )
            
            payload_out = [];

            if( ~obj.started )
                obj.Start(ECG_header);
                if( ~obj.started )
                    cprintf('*[1,0.5,0]', 'Task %s unable to be started for %s.\n', obj.name, ECG_header.recname);
                    return
                end
            end
            
            if(obj.bWFDBdetectors)
                % MIT conversion is needed for WFDB detectors.
                
                ECG_header.recname(ECG_header.recname == ' ') = '_';
                
                MIT_filename = [ECG_header.recname '_' num2str(ECG_start_offset) '_' num2str(ECG_header.nsamp+ECG_start_offset-1) ];
                ECG_header.recname = MIT_filename;
                MIT_filename = [obj.tmp_path_local MIT_filename '.dat'];
                fidECG = fopen(MIT_filename, 'w');
                try
                    fwrite(fidECG, ECG', 'int16', 0 );
                    fclose(fidECG);
                catch MEE
                    fclose(fidECG);
                    rethrow(MEE);
                end

                writeheader(obj.tmp_path_local, ECG_header);   

                if( isempty(obj.gqrs_config_filename) )
                    aux_str = obj.WFDB_bin_path;
                    if( aux_str(end) == filesep )
                        obj.gqrs_config_filename = [aux_str 'gqrs.conf' ];
                    else
                        obj.gqrs_config_filename = [aux_str filesep 'gqrs.conf' ];
                    end
                end
                
            end
            
            cant_QRSdetectors = length(obj.detectors2do);
            
            for ii = 1:cant_QRSdetectors

                this_detector = obj.detectors2do{ii};
                
                [this_detector, this_detector_name] = strtok(this_detector, ':');
                if( isempty(this_detector_name) )
                    this_detector_name = this_detector;
                else
                    this_detector_name = this_detector_name(2:end);
                end

                cprintf( 'Blue', [ 'Processing QRS detector ' this_detector_name '\n' ] );
                
                %% perform QRS detection

                switch( this_detector )

                    case 'wavedet'
                    %% Wavedet delineation

                        try

                            [position_multilead, positions_single_lead] = wavedet_interface(ECG, ECG_header, [], obj.lead_idx, obj.wavedet_config, ECG_sample_start_end_idx, ECG_start_offset, obj.progress_handle);

                            for jj = 1:length(obj.lead_idx)
                                % QRS detections in milliseconds
                                str_aux = strrep(strtrim(ECG_header.desc(obj.lead_idx(jj),:)), ' ', '_' );
                                payload_out.(['wavedet_' str_aux ]).time = colvec(positions_single_lead(jj).qrs);
                            end
                            
                            % QRS detections in milliseconds
                            payload_out.wavedet_multilead.time = colvec(position_multilead.qrs);

                        catch aux_ME
                            

                            strAux = sprintf('Wavedet failed in recording %s\n', ECG_header.recname);
                            strAux = sprintf('%s\n', strAux);

                            report = getReport(aux_ME);
                            fprintf(2, '%s\nError report:\n%s', strAux, report);

                        end

                    case 'pantom'
                        %% Pan-Tompkins detector

                        for jj = rowvec(obj.lead_idx)

                            obj.progress_handle.checkpoint(['Lead ' ECG_header.desc(jj,:)])
                            
                            try

                                aux_val = PeakDetection2(double(ECG(:,jj)), ECG_header.freq, [], [], [], [], [], obj.progress_handle);

                                % filter and offset
                                aux_val = aux_val( aux_val >= ECG_sample_start_end_idx(1) &  aux_val <= ECG_sample_start_end_idx(2) ) + ECG_start_offset - 1;
                                
                                str_aux = strrep(strtrim(ECG_header.desc(jj,:)), ' ', '_' );
                                % QRS detections in milliseconds
                                payload_out.(['pantom_' str_aux ]).time = colvec(aux_val);
                                
                            catch aux_ME
                                

                                strAux = sprintf('Pan & Tompkins failed in recording %s lead %s\n', ECG_header.recname, ECG_header.desc(jj,:) );
                                strAux = sprintf('%s\n', strAux);

                                report = getReport(aux_ME);
                                fprintf(2, 'Error report:\n%s', report);

                            end

                        end
                        
                    case 'gqrs'
                    %% WFDB_comp_interface (sqrs, wqrs, aristotle, ecgpuwave)
                        
                        % QRS detection.                        
                        for jj = rowvec(obj.lead_idx)

                            file_name_orig =  [obj.tmp_path_local ECG_header.recname '.' this_detector num2str(jj) ];

                            try

                                if( exist(obj.gqrs_config_filename, 'file') )
                                    % using the configuration file to
                                    % post-process
                                    [status, ~] = system([ obj.WFDB_cmd_prefix_str 'gqrs -c ' obj.gqrs_config_filename ' -r ' ECG_header.recname ' -s ' num2str(jj-1)]);
                                    if( status ~= 0 ); error('ECGtask_QRS_detection:WFDB_error', 'Error calling WFDB function.'); end
                                    
                                    [status, ~] = system([ obj.WFDB_cmd_prefix_str 'gqpost -c ' obj.gqrs_config_filename ' -r ' ECG_header.recname ' -o ' this_detector num2str(jj) ]);
                                    if( status ~= 0 ); error('ECGtask_QRS_detection:WFDB_error', 'Error calling WFDB function.'); end
                                    
                                    delete([obj.tmp_path_local ECG_header.recname '.qrs' ]);
                                else
                                    % run only gqrs
                                    [status, ~] = system([ obj.WFDB_cmd_prefix_str  'gqrs -r ' ECG_header.recname ' -o ' this_detector num2str(jj) ' -s ' num2str(jj-1)]);
                                    if( status ~= 0 ); error('ECGtask_QRS_detection:WFDB_error', 'Error calling WFDB function.'); end
                                end

                            catch aux_ME
                                

                                strAux = sprintf( '%s failed in recording %s lead %d\n', this_detector, ECG_header.recname, jj);

                                report = getReport(aux_ME);
                                fprintf(2, '%s\nError report:\n%s', strAux, report);

                            end

                            if( exist(file_name_orig, 'file') )
                                % reference comparison
                                anns_test = [];                    
                                try
                                    anns_test = readannot(file_name_orig);
                                    
                                    if( isempty(anns_test) )

                                        str_aux = strrep(strtrim(ECG_header.desc(jj,:)), ' ', '_' );
                                        payload_out.([this_detector '_' str_aux ]).time = [];

                                    else
                                        
                                        anns_test = AnnotationFilterConvert(anns_test, 'MIT', 'AAMI');

                                        % filter and offset
                                        bAux = anns_test.time >= ECG_sample_start_end_idx(1) & anns_test.time <= ECG_sample_start_end_idx(2);
                                        for fname = rowvec(fieldnames(anns_test) )
                                            aux_val = anns_test.(fname{1});
                                            anns_test.(fname{1}) = aux_val( bAux );
                                        end
                                        anns_test.time = anns_test.time + ECG_start_offset - 1;
                                        
                                        str_aux = strrep(strtrim(ECG_header.desc(jj,:)), ' ', '_' );
                                        payload_out.([this_detector '_' str_aux ]).time = colvec(anns_test.time);

                                    end
                                    
                                catch aux_ME
                                    
                                    
                                    if( strcmpi(aux_ME.identifier, 'MATLAB:nomem') )
                                        payload_out.([this_detector '_' num2str(jj)]).time = [];
                                    else
                                        strAux = sprintf( '%s failed in recording %s lead %s\n', this_detector, ECG_header.recname, ECG_header.desc(jj,:) );

                                        report = getReport(aux_ME);
                                        fprintf(2, '%s\nError report:\n%s', strAux, report);
                                    end
                                end

                            end

                        end

                    case 'user'
                        %% user-defined QRS detector

                        try

                            if( exist(this_detector_name) == 2 )
                                
                                ud_func_pointer = eval(['@' this_detector_name]);

                                obj.progress_handle.checkpoint([ 'User defined function: ' this_detector_name])

                                [positions_single_lead, position_multilead] = ud_func_pointer( double(ECG(:,obj.lead_idx)), trim_ECG_header(ECG_header, obj.lead_idx), obj.progress_handle, obj.payload);

                                for jj = 1:length(obj.lead_idx)

                                    aux_val = positions_single_lead{jj};
                                    % filter and offset
                                    aux_val = aux_val( aux_val >= ECG_sample_start_end_idx(1) &  aux_val <= ECG_sample_start_end_idx(2) ) + ECG_start_offset - 1;

                                    str_aux = strrep(strtrim(ECG_header.desc(obj.lead_idx(jj),:)), ' ', '_' );
                                    % QRS detections in milliseconds
                                    payload_out.([ this_detector_name '_' str_aux ]).time = colvec(aux_val);

                                end

                                if( ~isempty(position_multilead) )
                                    % filter and offset
                                    position_multilead = position_multilead( position_multilead >= ECG_sample_start_end_idx(1) &  position_multilead <= ECG_sample_start_end_idx(2) ) + ECG_start_offset - 1;
                                    payload_out.([ this_detector_name '_multilead' ]).time = colvec(position_multilead);
                                end
                                
                            else
                                disp_string_framed(2, sprintf('Function "%s" is not reachable in path.', this_detector_name));  
                                fprintf(1, 'Make sure that exist(%s) == 2\n',this_detector_name);
                            end
                            
                        catch aux_ME

                            disp_string_framed(2, sprintf('Detector "%s" failed in recording %s lead %s', this_detector_name, ECG_header.recname, ECG_header.desc(jj,:) ) );                                

                            report = getReport(aux_ME);
                            fprintf(2, 'Error report:\n%s', report);

                        end

                end
    
            end
            
            % calculate artificial leads
            if( obj.CalculateArtificialDetections )
                fn = fieldnames(payload_out);
                % more than one ECG signal
                if( length(fn) > 1 )
                    payload_out = calculate_artificial_QRS_detections(payload_out, ECG_header, [1 ECG_header.nsamp] + ECG_start_offset - 1 );
                end
            end
            
            % delete intermediate tmp files
            if(obj.bWFDBdetectors)
                delete([obj.tmp_path_local ECG_header.recname '.*' ]);
            end
            
            
        end
        
        function payload = Finish(obj, payload, ECG_header)
            
            if( obj.CalculateArtificialDetections )
                payload.series_quality.ratios = mean(payload.series_quality.ratios, 2);
            end
                
        end
        
        function payload = Concatenate(obj, plA, plB)

            if( isempty(plA) )
                
                payload = plB;
                
            else
                
                fields = rowvec(unique( [rowvec(fieldnames(plA)) rowvec(fieldnames(plB))] ));
                
                diff_fields = setdiff(fields, {'series_quality'});
                
                for fn = diff_fields
                    
                    if( isfield(plA, fn{1}) && isfield(plB, fn{1}) )
                        payload.(fn{1}).time = [ colvec(plA.(fn{1}).time); colvec(plB.(fn{1}).time) ];
                    
                    elseif( ~isfield(plA, fn{1}) && isfield(plB, fn{1}) )
                        payload.(fn{1}).time = colvec(plB.(fn{1}).time);
                    end
                    
                end
                
                aux_idx = find(strcmpi(fields, 'series_quality'));
                
                if( ~isempty(aux_idx) )
                    
                    if( isfield(plA, 'series_quality') && isfield(plB, 'series_quality') )
                        if( size(plA.series_quality.ratios,1) == size(plB.series_quality.ratios,1) )
                            payload.series_quality.ratios = [plA.series_quality.ratios plB.series_quality.ratios];
                            payload.series_quality.estimated_labs = cellfun(@(a,b)( [colvec(a);colvec(b)] ) , plA.series_quality.estimated_labs, plB.series_quality.estimated_labs, 'UniformOutput', false);
                            payload.series_quality.AnnNames = plA.series_quality.AnnNames;
                        else
                            payload.series_quality.ratios = plA.series_quality.ratios;
                            payload.series_quality.estimated_labs = plA.series_quality.estimated_labs;
                            payload.series_quality.AnnNames = plA.series_quality.AnnNames;
                        end
                    end

                end
                
            end
            

        end

        %% property restriction functions
        
        function set.detectors(obj,x)
            
            if( (ischar(x) || iscellstr(x)) )
                x = cellstr(x);
                aux_val = colvec(intersect(obj.cQRSdetectors, x));
                aux_idx = find(cellfun(@(a)(~isempty(strfind(a, 'user:'))), x));
                obj.detectors = [aux_val; colvec(x(aux_idx)) ];
            else
                warning('ECGtask_QRS_detection:BadArg', 'Invalid detectors.');
            end
        end
        
        function set.only_ECG_leads(obj,x)
            if( islogical(x) )
                obj.only_ECG_leads = x;
            else
                warning('ECGtask_QRS_detection:BadArg', 'Invalid lead configuration.');
            end
        end
        
        function set.gqrs_config_filename(obj,x)
            if( ischar(x) )
                
                if(exist(x, 'file'))
                    obj.gqrs_config_filename = x;
                else
                    warning('ECGtask_QRS_detection:BadArg', 'obj.gqrs_config_filename does not exist.');
                end
                
            else
                warning('ECGtask_QRS_detection:BadArg', 'obj.gqrs_config_filename must be a string.');
            end
        end
        
        function set.tmp_path_local(obj,x)
            if( ischar(x) )
                if(exist(x, 'dir'))
                    obj.tmp_path_local = x;
                else
                    if(mkdir(x))
                        obj.tmp_path_local = x;
                    else
                        warning('ECGtask_QRS_detection:BadArg', ['Could not create ' x ]);
                    end
                end
                
            else
                warning('ECGtask_QRS_detection:BadArg', 'tmp_path_local must be a string.');
            end
        end
        
        function set.tmp_path(obj,x)
            if( ischar(x) )
                if(exist(x, 'dir'))
                    obj.tmp_path = x;
                else
                    if(mkdir(x))
                        obj.tmp_path = x;
                    else
                        warning('ECGtask_QRS_detection:BadArg', ['Could not create ' x ]);
                    end
                end
                
            else
                warning('ECGtask_QRS_detection:BadArg', 'tmp_path_local must be a string.');
            end
        end
        
    end
    
    methods ( Access = private )
        
        
    end
    
end
