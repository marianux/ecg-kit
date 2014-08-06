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
        user_string = '';
        caller_variable = 'payload'
        payload
        mode 
        tmp_path        
        
    end
    
    methods
           
        function obj = ECGtask_heartbeat_classifier(obj)

            if( isempty(obj.mode) )
                obj.mode = obj.cKnownModesOfOperation{1};
            end
            
        end
        
        function Start(obj, ECG_header, ECG_annotations)

        end
        
        function payload_out = Process(obj, ECG, ECG_start_offset, ECG_sample_start_end_idx, ECG_header, ECG_annotations, ECG_annotations_start_end_idx )
            
            payload_out = [];
            
            aux_val = ECG_sample_start_end_idx + ECG_start_offset - 1;
            disp_string_title(1, sprintf( 'Classifying from %d to %d', aux_val(1), aux_val(2) ) );
            
            if( isempty(obj.payload) )
                Ann_struct = ECG_annotations;
            elseif( isstruct(obj.payload) )
                
                aux_fn = fieldnames(obj.payload);
                % prefer manually reviewed annotations corrected 
                aux_idx = find(cell2mat( cellfun(@(a)(~isempty(strfind(a, 'corrected_'))), aux_fn, 'UniformOutput', false)));
                
                if( isempty(aux_idx) )
                    aux_idx = find(cell2mat( cellfun(@(a)(~isempty(strfind(a, 'wavedet_'))), aux_fn, 'UniformOutput', false)));
                    if( isempty(aux_idx) )
                        disp_string_framed(2, 'Could not identify QRS detections in the input payload.');
                        return
                    else
                        aux_val = obj.payload.(aux_fn{1}).time;
                    end
                else
                    aux_val = obj.payload.(aux_fn{aux_idx(1)}).time;
                end
                
                aux_val = aux_val - ECG_start_offset + 1 ;
                bAux = aux_val >= ECG_sample_start_end_idx(1) & aux_val <= ECG_sample_start_end_idx(2);
                Ann_struct.time = aux_val(bAux);
            else
                obj.payload = obj.payload - ECG_start_offset + 1 ;
                bAux = obj.payload >= ECG_sample_start_end_idx(1) & obj.payload <= ECG_sample_start_end_idx(2);
                Ann_struct.time = obj.payload(bAux);
            end
            
            [aux_val, ~, lablist ] = a2hbc('ECG', ECG, ...
                            'ECG_header', ECG_header, ...
                            'ECG_annotations', Ann_struct, ...
                            'op_mode', obj.mode );

            Ann_struct.anntyp = lablist(colvec(aux_val),1);
            Ann_struct.time = Ann_struct.time + ECG_start_offset - 1;
            
            payload_out = Ann_struct;
            
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
        
        
    end
    
end
