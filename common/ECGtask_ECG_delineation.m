classdef ECGtask_ECG_delineation < ECGtask

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
        name = 'ECG_delineation';
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
        
        cQRSdelineators = {'all-delineators' 'wavedet' };
        
    end
    
    properties( Access = private )
        
        delineators2do
        bWFDBdelineators
        tmp_path_local
        tmp_path
        command_sep
        WFDB_bin_path
        
    end
    
    properties
        progress_handle
        user_string = '';
        delineators = 'all-delineators';
        wavedet_config
        CalculateArtificialDetections = true;
    end
    
    methods
           
        function obj = ECGtask_ECG_delineation(obj)

            % wavedet configuration structure
            obj.wavedet_config.setup.wavedet.QRS_detection_only = false;
            
        end
        
        function Start(obj, ECG_header, ECG_annotations)

            if( strcmpi('all-delineators', obj.delineators) )
                obj.delineators2do = obj.cQRSdelineators(2:end);
            else
                if( ischar(obj.delineators) )
                    obj.delineators2do = {obj.delineators};
                else
                    obj.delineators2do = obj.delineators;
                end
            end

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
            
        end
        
        function payload = Process(obj, ECG, ECG_start_offset, ECG_sample_start_end_idx, ECG_header, ECG_annotations, ECG_annotations_start_end_idx, payload_in  )
            
            payload = [];
            
            if( nargin < 8 )
                payload_in = [];
            else
                % payload_in is used in this task to input an external QRS
                % detector, or manually corrected detections.
                if( isstruct(payload_in) )
                    
                    fnames = fieldnames(payload_in);
                    aux_idx = find(cell2mat( cellfun(@(a)(~isempty(strfind(a, 'corrected_'))), fnames, 'UniformOutput', false)));
                    if( isempty(aux_idx) )
                        payload_in = [];
                    else
                        payload_in = payload_in.(fnames{aux_idx(1)}).time;
                        payload_in = payload_in( payload_in >= ECG_sample_start_end_idx(1) & payload_in < ECG_sample_start_end_idx(2) );
                    end
                    
                else
                    payload_in = [];
                end
            end
            
            % get signals indexes that actually are ECG
            ECG_signals_idx = get_ECG_idx_from_header(ECG_header);
            
            if( isempty(ECG_signals_idx) )
                warning('ECGtask_delineation:NoECGsignals', disp_option_enumeration( 'No ECG signals found:', cellstr(ECG_header.desc) ) );
                ECG_signals_idx = 1:ECG_header.nsig;
            end
            
            cant_QRSdelineators = length(obj.delineators2do);
            
            for ii = 1:cant_QRSdelineators

                this_delineator = obj.delineators2do{ii};

                fprintf(1, [ '\nProcessing ' this_delineator '\n\n' ] );

                %% perform ECG delineation

                switch( this_delineator )

                    case 'wavedet'
                    %% Wavedet delineation

                        try

                            [position_multilead, positions_single_lead] = wavedet_interface(ECG, ECG_header, payload_in, ECG_signals_idx, obj.wavedet_config, ECG_sample_start_end_idx, ECG_start_offset, obj.progress_handle);

                            payload.wavedet_single_lead = positions_single_lead;
                            
                            % QRS detections in milliseconds
                            payload.wavedet_multilead = position_multilead;

                        catch aux_ME

                            strAux = sprintf('Wavedet failed in recording %s\n', ECG_header.recname);
                            strAux = sprintf('%s\n', strAux);

                            report = getReport(aux_ME);
                            fprintf(2, '%s\nError report:\n%s', strAux, report);

                        end


%                     case XXXX
                        %% Room for future delineators.


                end
    
            end
            
        end
        
        function payload = Finish(obj, payload, ECG_header)

            if( obj.CalculateArtificialDetections )
                % future implementation of an strategy for the construction of
                % artificial delineations
%                 payload = calculate_artificial_QRS_detections(payload, ECG_header );
            end
            
        end
        
        function payload = Concatenate(obj, plA, plB)

            if( isempty(plA) )
                
                payload = plB;
                
            else
                for fn = rowvec(fieldnames(plA))
                    
                    if( isfield(plA, fn{1}) && isfield(plB, fn{1}) )
                        payload.(fn{1}) = [ colvec(plA.(fn{1})); colvec(plB.(fn{1})) ];
                    
                    elseif( ~isfield(plA, fn{1}) && isfield(plB, fn{1}) )
                        payload.(fn{1}) = colvec(plB.(fn{1}));
                    end
                    
                end            
            end
            

        end

        %% property restriction functions
        
        function set.delineators(obj,x)
            if( (ischar(x) || iscellstr(x)) && ~isempty(intersect(obj.cQRSdelineators, cellstr(x))) )
                obj.delineators = x;
            else
                warning('ECGtask_QRS_detection:BadArg', 'Invalid delineators.');
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
