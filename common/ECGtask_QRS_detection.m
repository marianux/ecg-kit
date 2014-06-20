classdef ECGtask_QRS_detection < ECGtask

% ECGtask for ECGwrapper (for Matlab)
% ---------------------------------
% 
% Description:
% 
% Abstract class for defining ECGtask interface
% 
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
    end
    
    properties( Access = private, Constant)
        
        cQRSdetectors = {'all-detectors' 'wavedet' 'pantom' 'gqrs' };
        cLeadConfiguration = {'all-leads' 'ECG-leads-only' 'User-defined-leads' };
    end
    
    properties( Access = private )
        
        detectors2do
        bWFDBdetectors
        tmp_path_local
        tmp_path
        WFDB_bin_path
        WFDB_cmd_prefix_str 
        
    end
    
    properties
        progress_handle
        user_string = '';
        detectors = 'all-detectors';
        lead_config = 'all-leads';
        lead_idx = [];
        gqrs_config_filename = [];
        wavedet_config
        CalculateArtificialDetections = true;
    end
    
    methods
           
        function obj = ECGtask_QRS_detection(obj)

            obj.wavedet_config.setup.wavedet.QRS_detection_only = true;
            
        end
        
        function Start(obj, ECG_header, ECG_annotations)

            if( strcmpi('all-leads', obj.lead_config) )
                obj.lead_idx = 1:ECG_header.nsig;
            elseif( strcmpi('ECG-leads-only', obj.lead_config) )
                obj.lead_idx = get_ECG_idx_from_header(ECG_header);
            else
                if( isempty(obj.lead_idx) || ~isnumeric(obj.lead_idx) || ~(all(obj.lead_idx > 0) && all(obj.lead_idx <= ECG_header.nsig)) )
                    warning('ECGtask_QRS_detection:BadArg', 'Invalid lead indexes. Indexes between 1 and %d\n', ECG_header.nsig);
                    obj.lead_idx = 1:ECG_header.nsig;
                end
            end
            
            if( strcmpi('all-detectors', obj.detectors) )
                obj.detectors2do = obj.cQRSdetectors(2:end);
            else
                if( ischar(obj.detectors) )
                    obj.detectors2do = {obj.detectors};
                else
                    obj.detectors2do = obj.detectors;
                end
            end

            obj.bWFDBdetectors = any(strcmpi('gqrs', obj.cQRSdetectors));
            
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
                obj.tmp_path_local = tempdir;
                obj.tmp_path = tempdir;
            end
            
            obj.WFDB_bin_path = init_WFDB_library(obj.tmp_path_local);
            
            obj.WFDB_cmd_prefix_str = WFDB_command_prefix(obj.tmp_path_local);
            
        end
        
        function payload = Process(obj, ECG, ECG_start_offset, ECG_sample_start_end_idx, ECG_header, ECG_annotations, ECG_annotations_start_end_idx, payload_in  )
            
            payload = [];

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

                fprintf(1, [ '\nProcessing ' this_detector '\n\n' ] );

                %% perform QRS detection

                switch( this_detector )

                    case 'wavedet'
                    %% Wavedet delineation

                        try

                            [position_multilead, positions_single_lead] = wavedet_interface(ECG, ECG_header, [], obj.lead_idx, obj.wavedet_config, ECG_sample_start_end_idx, ECG_start_offset, obj.progress_handle);

                            for jj = 1:length(obj.lead_idx)
                                % QRS detections in milliseconds
                                payload.(['wavedet_' num2str(obj.lead_idx(jj))]).time = colvec(positions_single_lead(jj).qrs);
                            end
                            
                            % QRS detections in milliseconds
                            payload.wavedet_multilead.time = colvec(position_multilead.qrs);

                        catch aux_ME

                            strAux = sprintf('Wavedet failed in recording %s\n', ECG_header.recname);
                            strAux = sprintf('%s\n', strAux);

                            report = getReport(aux_ME);
                            fprintf(2, '%s\nError report:\n%s', strAux, report);

                        end


                    case 'pantom'
                        %% Pan-Tompkins detector

                        for jj = rowvec(obj.lead_idx)

                            obj.progress_handle.checkpoint(['Lead ' ECG_header.desc(ii,:)])
                            
                            try

                                aux_val = PeakDetection2(double(ECG(:,jj)), ECG_header.freq, [], [], [], [], [], obj.progress_handle);

                                aux_val = aux_val( aux_val >= ECG_sample_start_end_idx(1) &  aux_val <= ECG_sample_start_end_idx(2) ) - ECG_start_offset + 1;
                                
                                % QRS detections in milliseconds
                                payload.(['pantom_' num2str(jj)]).time = colvec(aux_val);
                                
                            catch aux_ME

                                strAux = sprintf('Pan & Tompkins failed in recording %s lead %d\n', ECG_header.recname, jj);
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
                                    status = system([ obj.WFDB_cmd_prefix_str 'gqrs -c ' obj.gqrs_config_filename ' -r ' ECG_header.recname ' -s ' num2str(jj-1)]);
                                    if( status ~= 0 ); error('ECGtask_QRS_detection:WFDB_error', 'Error calling WFDB function.'); end
                                    
                                    status = system([ obj.WFDB_cmd_prefix_str 'gqpost -c ' obj.gqrs_config_filename ' -r ' ECG_header.recname ' -o ' this_detector num2str(jj) ]);
                                    if( status ~= 0 ); error('ECGtask_QRS_detection:WFDB_error', 'Error calling WFDB function.'); end
                                    
                                    delete([obj.tmp_path_local ECG_header.recname '.qrs' ]);
                                else
                                    % run only gqrs
                                    system([ obj.WFDB_cmd_prefix_str  'gqrs -r ' ECG_header.recname ' -o ' this_detector num2str(jj) ' -s ' num2str(jj-1)]);
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

                                        payload.([this_detector '_' num2str(jj)]).time = [];

                                    else
                                        
                                        anns_test = AnnotationFilterConvert(anns_test, 'MIT', 'AAMI');

                                        bAux = anns_test.time >= ECG_sample_start_end_idx(1) & anns_test.time <= ECG_sample_start_end_idx(2);
                                        for fname = rowvec(fieldnames(anns_test) )
                                            aux_val = anns_test.(fname{1});
                                            anns_test.(fname{1}) = aux_val( bAux );
                                        end
                                        anns_test.time = anns_test.time - ECG_start_offset + 1;

                                        payload.([this_detector '_' num2str(jj)]).time = colvec(anns_test.time);

                                    end
                                    
                                catch aux_ME
                                    if( strcmpi(aux_ME.identifier, 'MATLAB:nomem') )
                                        payload.([this_detector '_' num2str(jj)]).time = [];
                                    else
                                        strAux = sprintf( '%s failed in recording %s lead %d\n', this_detector, ECG_header.recname, jj);

                                        report = getReport(aux_ME);
                                        fprintf(2, '%s\nError report:\n%s', strAux, report);
                                    end
                                end

                            end

                        end

                end
    
            end
            
            % delete intermediate tmp files
            if(obj.bWFDBdetectors)
                delete([obj.tmp_path_local ECG_header.recname '.*' ]);
            end
            
            
        end
        
        function payload = Finish(obj, payload, ECG_header)

            if( obj.CalculateArtificialDetections )
                payload = calculate_artificial_QRS_detections(payload, ECG_header );
            end
            
        end
        
        function payload = Concatenate(obj, plA, plB)

            if( isempty(plA) )
                
                payload = plB;
                
            else
                for fn = rowvec(fieldnames(plA))
                    
                    if( isfield(plA, fn{1}) && isfield(plB, fn{1}) )
                        payload.(fn{1}).time = [ plA.(fn{1}).time; plB.(fn{1}).time ];
                    
                    elseif( ~isfield(plA, fn{1}) && isfield(plB, fn{1}) )
                        payload.(fn{1}).time = plB.(fn{1}).time;
                    end
                    
                end            
            end
            

        end

        %% property restriction functions
        
        function set.detectors(obj,x)
            if( (ischar(x) || iscellstr(x)) && ~isempty(intersect(obj.cQRSdetectors, cellstr(x))) )
                obj.detectors = x;
            else
                warning('ECGtask_QRS_detection:BadArg', 'Invalid detectors.');
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
