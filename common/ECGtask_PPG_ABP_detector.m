classdef ECGtask_PPG_ABP_detector < ECGtask

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
% Birthdate  : 18/7/2014
% Last update: 18/7/2014
       
    properties(GetAccess = public, Constant)
        name = 'PPG_ABP_detector';
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
        
        cPPGdetectors =      {'all-detectors' 'wavePPG' 'wabp' };
        cLeadConfiguration = {'all-leads' 'PPG-ABP-only' 'User-defined-leads' };
        
    end
    
    properties( Access = private )
        
        detectors2do
        bWFDBdetectors
        tmp_path_local
        WFDB_bin_path
        WFDB_cmd_prefix_str 
    end
    
    properties
        progress_handle
        tmp_path
        detectors = 'all-detectors';
        lead_config = 'PPG-ABP-only';
        PPG_ABP_idx = [];
        gqrs_config_filename = [];
        wavedet_config
        
        
    end
    
    methods
           
        function obj = ECGtask_PPG_ABP_detector(obj)
            
        end
        
        function Start(obj, ECG_header, ECG_annotations)
            
            if( strcmpi('all-leads', obj.lead_config) )
                obj.lead_config = 'all-leads';
                obj.PPG_ABP_idx = 1:ECG_header.nsig;
            elseif( strcmpi('PPG-ABP-only', obj.lead_config) )
                obj.PPG_ABP_idx = get_PPG_ABP_idx_from_header(ECG_header);
            else
                if( isempty(obj.PPG_ABP_idx) || ~isnumeric(obj.PPG_ABP_idx) || ~(all(obj.PPG_ABP_idx > 0) && all(obj.PPG_ABP_idx <= ECG_header.nsig)) )
                    warning('ECGtask_QRS_detection:BadArg', 'Invalid lead indexes. Indexes between 1 and %d\n', ECG_header.nsig);
                    obj.lead_config = 'all-leads';
                    obj.PPG_ABP_idx = 1:ECG_header.nsig;
                else
                    obj.lead_config = 'User-defined-leads';
                end
            end
            
            if( isempty(obj.PPG_ABP_idx) )
                return
            end
            
            if( strcmpi('all-detectors', obj.detectors) )
                obj.detectors2do = obj.cPPGdetectors(2:end);
            else
                if( ischar(obj.detectors) )
                    obj.detectors2do = {obj.detectors};
                else
                    obj.detectors2do = obj.detectors;
                end
            end

            obj.bWFDBdetectors = any(strcmpi('wabp', obj.cPPGdetectors));
            
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
        
        function payload = Process(obj, PPG_ABP, ECG_start_offset, ECG_sample_start_end_idx, ECG_header, ECG_annotations, ECG_annotations_start_end_idx, payload_in  )
            
            payload = [];
            
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
                    fwrite(fidECG, PPG_ABP', 'int16', 0 );
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
            
            cant_PPG_ABP_detectors = length(obj.detectors2do);
            
            for ii = 1:cant_PPG_ABP_detectors

                this_detector = obj.detectors2do{ii};

                cprintf( 'Blue', [ 'Processing PPG/ABP detector ' this_detector '\n' ] );
                
                %% perform pulse detection

                switch( this_detector )

                    case 'wavePPG'
                    %% wavePPG delineation

                        for jj = rowvec(obj.PPG_ABP_idx)
                            
                            try

                                PPG_ABP_peak_idx = PPG_pulses_detector( double(PPG_ABP(:,jj)), ECG_header.freq );

                                PPG_ABP_peak_idx = PPG_ABP_peak_idx( PPG_ABP_peak_idx >= ECG_sample_start_end_idx(1) & PPG_ABP_peak_idx <= ECG_sample_start_end_idx(2)) + ECG_start_offset - 1 ;

                                % QRS detections in milliseconds
                                str_aux = strrep(strtrim(ECG_header.desc(jj,:)), ' ', '_' );
                                payload.(['wavePPG_' str_aux]).time = colvec(PPG_ABP_peak_idx);

                            catch aux_ME
                                
                                strAux = sprintf('WavePPG failed in recording %s\n', ECG_header.recname);
                                strAux = sprintf('%s\n', strAux);

                                report = getReport(aux_ME);
                                fprintf(2, '%s\nError report:\n%s', strAux, report);

                            end
                            
                        end
                        
                    case 'wabp'
                    %% wabp WFDB interface

                        file_name_orig =  [obj.tmp_path_local ECG_header.recname '.wabp' ];
                    
                        % QRS detection.                        
                        for jj = rowvec(obj.PPG_ABP_idx)

                            if( exist(file_name_orig, 'file') )
                                delete(file_name_orig)
                            end
                            
                            try

                                % run wabp
                                [status, ~] = system([ obj.WFDB_cmd_prefix_str  'wabp -r ' ECG_header.recname ' -s ' num2str(jj-1)]);
                                if( status ~= 0 ); error('ECGtask_PPG_detector:WFDB_error', 'Error calling WFDB function.'); end

                            catch aux_ME
                                
                                strAux = sprintf( '%s failed in recording %s in %s\n', this_detector, ECG_header.recname, ECG_header.desc(jj,:) );

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
                                        payload.([this_detector '_' str_aux ]).time = [];

                                    else
                                        
                                        anns_test = AnnotationFilterConvert(anns_test, 'MIT', 'AAMI');

                                        PPG_ABP_peak_idx = anns_test.time;

                                        PPG_ABP_peak_idx = PPG_ABP_peak_idx( PPG_ABP_peak_idx >= ECG_sample_start_end_idx(1) & PPG_ABP_peak_idx <= ECG_sample_start_end_idx(2)) + ECG_start_offset - 1 ;
                                        
                                        str_aux = strrep(strtrim(ECG_header.desc(jj,:)), ' ', '_' );
                                        payload.([this_detector '_' str_aux ]).time = colvec(PPG_ABP_peak_idx);

                                    end
                                    
                                catch aux_ME
                                    
                                    if( strcmpi(aux_ME.identifier, 'MATLAB:nomem') )
                                        payload.([this_detector '_' ECG_header.desc(jj,:)]).time = [];
                                    else
                                        strAux = sprintf( '%s failed in recording %s in %s\n', this_detector, ECG_header.recname, ECG_header.desc(jj,:) );

                                        report = getReport(aux_ME);
                                        fprintf(2, '%s\nError report:\n%s', strAux, report);
                                    end
                                end

                            end

                        end

                end
    
            end
            
            % Add QRS detections quality metrics, Names, etc.
            payload = calculateSeriesQuality(payload, ECG_header, [1 ECG_header.nsamp] + ECG_start_offset - 1 );
            
            % delete intermediate tmp files
            if(obj.bWFDBdetectors)
                delete([obj.tmp_path_local ECG_header.recname '.*' ]);
            end
            
            
        end
        
        function payload = Finish(obj, payload, ECG_header)
  
            if( isfield(payload, 'series_quality') && isfield(payload.series_quality, 'ratios') )
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
            if( (ischar(x) || iscellstr(x)) && ~isempty(intersect(obj.cPPGdetectors, cellstr(x))) )
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
