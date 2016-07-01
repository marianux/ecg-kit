classdef ECGtask_heartbeat_classifier < ECGtask

% ECGtask for ECGwrapper (for Matlab)
% ---------------------------------
% 
% Description:
% 
% Abstract class for defining ECGtask interface
% 
% 
% Author: Mariano Llamedo Soria (llamedom at frba.utn.edu.ar)
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
            
            
            if(( isfield(ECG_annotations, 'time') && ~isempty(ECG_annotations.time)) ) 
                obj.started = true;
            else
                
                % payload property is used in this task to input an external QRS
                % detector, or manually corrected detections.
                aux_val = get_QRS_from_payload(obj.payload);
                
                if(~isempty(aux_val) && length(aux_val.time) > obj.min_heartbeats_required ) 
                    obj.started = true;
                else
                    warning('ECGtask_HB_classifier:BadPayload', 'Not enough heartbeats to perform classification. Try decreasing "min_heartbeats_required" or increasing the amount of heartbeats.\n');
                end
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
            
            % payload property is used in this task to input an external QRS
            % detector, or manually corrected detections.
            aux_QRS_loc = get_QRS_from_payload(obj.payload, ECG_start_offset, ECG_sample_start_end_idx);
            
            if(isempty(aux_QRS_loc))
                annotations_used = ECG_annotations;
            else
                annotations_used = aux_QRS_loc;
            end
            
%             aux_val = ECG_sample_start_end_idx + ECG_start_offset - 1;
%             disp_string_title(1, sprintf( 'Classifying from %d to %d', aux_val(1), aux_val(2) ) );
            
            [lead_idx, ECG_header ] = get_ECG_idx_from_header(ECG_header);
            
            [aux_val, ~, lablist ] = a2hbc('ECG', ECG(:,lead_idx), ...
                            'ECG_header', ECG_header, ...
                            'ECG_annotations', annotations_used, ...
                            'op_mode', obj.mode );

            bAux = annotations_used.time > ECG_sample_start_end_idx(1) & annotations_used.time < ECG_sample_start_end_idx(2);
            payload_out.anntyp = lablist(colvec(aux_val(bAux)),1);
            payload_out.time = annotations_used.time(bAux) + ECG_start_offset - 1;
            payload_out.lablist = lablist;
            
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
                warning('ECGtask_HB_classifier:BadArg', 'Invalid mode.');
            end
        end

        
    end
    
    methods ( Access = private )
          
    end
    
end
