classdef ECGtask_heartbeat_classifier < ECGtask

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
        name = 'ECG_heartbeat_classifier';
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
        
        fig_hdl = 1;
        cKnownModesOfOperation = {'auto' 'slightly-assisted', 'assisted'};
        
    end
    
    properties( Access = private )
        
        
    end
    
    properties
        
        progress_handle
        payload
        mode 
        tmp_path
        
        % minimum amount of heartbeats required
        min_heartbeats_required = 500;
        % maximum amount of heartbeats to process each iteration
        max_heartbeats_per_iter = 4000;
        % qrs locations provided to the task, result of previous QRS
        % detections / correction
        annotations = [];
                
    end
    
    methods
           
        function obj = ECGtask_heartbeat_classifier(obj)

            if( isempty(obj.mode) )
                obj.mode = obj.cKnownModesOfOperation{1};
            end
            
        end
        
        function Start(obj, ECG_header, ECG_annotations)

            if( ECG_header.nsig == 1 )
                return
            end

            lead_idx = get_ECG_idx_from_header(ECG_header);

            if( length(lead_idx) <= 1 )
                return
            end
            
            obj.annotations = ParsePayloadAnnotations(obj);
            
            if( (~isempty(obj.annotations) && length(obj.annotations) > obj.min_heartbeats_required ) || ( isfield(ECG_annotations, 'time') && ~isempty(ECG_annotations.time))  )
                obj.started = true;
            end
            
        end
        
        function payload_out = Process(obj, ECG, ECG_start_offset, ECG_sample_start_end_idx, ECG_header, ECG_annotations, ECG_annotations_start_end_idx )
            
            payload_out = [];
            
            if( ~obj.started )
                obj.Start(ECG_header);
                if( ~obj.started )
                    cprintf('*[1,0.5,0]', 'Task %s unable to be started for %s.\n', obj.name, ECG_header.recname);
                    return
                end
            end
            
%             aux_val = ECG_sample_start_end_idx + ECG_start_offset - 1;
%             disp_string_title(1, sprintf( 'Classifying from %d to %d', aux_val(1), aux_val(2) ) );
            
            [lead_idx, ECG_header ] = get_ECG_idx_from_header(ECG_header);
            
            [aux_val, ~, lablist ] = a2hbc('ECG', ECG(:,lead_idx), ...
                            'ECG_header', ECG_header, ...
                            'ECG_annotations', ECG_annotations, ...
                            'op_mode', obj.mode );

            bAux = ECG_annotations.time > ECG_sample_start_end_idx(1) & ECG_annotations.time < ECG_sample_start_end_idx(2);
            ECG_annotations.anntyp = lablist(colvec(aux_val(bAux)),1);
            ECG_annotations.time = ECG_annotations.time(bAux) + ECG_start_offset - 1;
            
            payload_out = ECG_annotations;
            
        end
        
        function payload = Finish(obj, payload, ECG_header)
           
        end
        
        function payload = Concatenate(obj, plA, plB)

            if( isempty(plA) )
                
                payload = plB;
                
            else
                for fn = rowvec(fieldnames(plA))
                    
                    if( isfield(plA, fn{1}) && isfield(plB, fn{1}) )
                        payload.(fn{1}) = [ plA.(fn{1}); plB.(fn{1}) ];
                    elseif( ~isfield(plA, fn{1}) && isfield(plB, fn{1}) )
                        payload.(fn{1}) = plB.(fn{1});
                    end
                    
                end            
            end
            
        end

        %% property restriction functions
        
        function set.mode(obj,x)
            if( ischar(x) && ~isempty(strcmpi(x, obj.cKnownModesOfOperation)) )
                obj.mode = x;
            else
                warning('ECGtask_QRS_detection:BadArg', 'Invalid detectors.');
            end
        end

        
    end
    
    methods ( Access = private )
        
        function aux_val = ParsePayloadAnnotations(obj)

            aux_val = [];
            
            if( isstruct(obj.payload) )
                
                aux_fn = fieldnames(obj.payload);
                
                jj = 1;
                
                while( isempty( aux_val ) )
                    
                    % prefer manually reviewed annotations corrected 
                    aux_idx = find(cell2mat( cellfun(@(a)(~isempty(strfind(a, 'corrected_'))), aux_fn, 'UniformOutput', false)));

                    if( isempty(aux_idx) )

                        if( isfield(obj.payload, 'series_quality') )
                            % choose the best ranked automatic detection
                            [~, aux_idx] = sort(obj.payload.series_quality.ratios, 'descend' );
                            if( jj <= length(aux_idx))
                                aux_val = obj.payload.(obj.payload.series_quality.AnnNames{aux_idx(jj),1} ).(obj.payload.series_quality.AnnNames{aux_idx(jj),2});
                                delineation_chosen = obj.payload.series_quality.AnnNames{aux_idx(jj),1};
                            else
                                break
                            end
                        else
                            % choose the first of wavedet
                            aux_idx = find(cell2mat( cellfun(@(a)(~isempty(strfind(a, 'wavedet_'))), aux_fn, 'UniformOutput', false)));
                            if( isempty(aux_idx) )
                                disp_string_framed(2, 'Could not identify QRS detections in the input payload.');
                                return
                            else
                                if( jj <= length(aux_idx))
                                    aux_val = obj.payload.(aux_fn{aux_idx(jj)}).time;
                                    delineation_chosen = aux_fn{aux_idx(jj)};
                                else
                                    break
                                end
                            end
                        end
                    else
                        % choose the first manually audited
                        if( jj <= length(aux_idx))
                            aux_val = obj.payload.(aux_fn{aux_idx(jj)}).time;
                            delineation_chosen = aux_fn{aux_idx(jj)};
                        else
                            break
                        end
                    end

                    cprintf('blue', [' + Using ' delineation_chosen ' detections.\n' ] );
                    aux_val = sort(unique(aux_val));
                    jj = jj + 1;
                end
                
            elseif( isnumeric(obj.payload) )
                aux_val = sort(unique(obj.payload));
            end
            
        end        
        
    end
    
end
