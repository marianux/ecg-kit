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
    end
    
    properties( Access = private, Constant)
        
        fig_hdl = 1;
        cKnownModesOfOperation = {'auto' 'slightly-assisted', 'assisted'};
        
    end
    
    properties( Access = private )
        
        tmp_path
        
    end
    
    properties
        
        progress_handle
        user_string = '';
        caller_variable = 'payload'
        mode 
        
    end
    
    methods
           
        function obj = ECGtask_heartbeat_classifier(obj)

            if( isempty(obj.mode) )
                obj.mode = obj.cKnownModesOfOperation{1};
            end
            
        end
        
        function Start(obj, ECG_header, ECG_annotations)

        end
        
        function payload = Process(obj, ECG, ECG_start_offset, ECG_sample_start_end_idx, ECG_header, ECG_annotations, ECG_annotations_start_end_idx, payload_in )
            
            payload = [];
            
            aux_val = ECG_sample_start_end_idx + ECG_start_offset - 1;
            disp_string_title(1, sprintf( 'Correcting from %d to %d', aux_val(1), aux_val(2) ) );
            
            if( isempty(payload_in) )
                Ann_struct = ECG_annotations;
            elseif( isstruct(payload_in) )
                
                aux_fn = fieldnames(payload_in);
                % prefer manually reviewed annotations corrected 
                aux_idx = find(cell2mat( cellfun(@(a)(~isempty(strfind(a, 'corrected_'))), aux_fn, 'UniformOutput', false)));
                
                if( isempty(aux_idx) )
                    aux_idx = find(cell2mat( cellfun(@(a)(~isempty(strfind(a, 'wavedet_'))), aux_fn, 'UniformOutput', false)));
                    if( isempty(aux_idx) )
                        disp_string_framed(2, 'Could not identify QRS detections in the input payload.');
                        return
                    else
                        Ann_struct.time = payload_in.(aux_fn{1}).time;
                    end
                else
                    Ann_struct.time = payload_in.(aux_fn{aux_idx(1)}).time;
                end
                
            else
                Ann_struct = payload_in;
            end
            
            [aux_val, ~, lablist ] = a2hbc('ECG', ECG, ...
                            'ECG_header', ECG_header, ...
                            'ECG_annotations', Ann_struct, ...
                            'op_mode', obj.mode );

            Ann_struct.anntyp = lablist(colvec(aux_val),1);
            
            payload = Ann_struct;
            
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
