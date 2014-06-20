classdef ECGtask_QRS_corrector < ECGtask

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
        name = 'QRS_corrector';
        target_units = 'ADCu';
        doPayload = true;
    end

    properties( GetAccess = public, SetAccess = private)
        % if user = memory;
        % memory_constant is the fraction respect to user.MaxPossibleArrayBytes
        % which determines the maximum input data size.
        memory_constant = 0.8;
    end
    
    properties( Access = private, Constant)
        
        fig_hdl = 1;
        
    end
    
    properties( Access = private )
        
        tmp_path
        
    end
    
    properties
        
        progress_handle
        user_string = '';
        caller_variable = 'payload'
        
    end
    
    methods
           
        function obj = ECGtask_QRS_detection(obj)

            obj.wavedet_config.setup.wavedet.QRS_detection_only = true;
            
        end
        
        function Start(obj, ECG_header, ECG_annotations)

        end
        
        function payload = Process(obj, ECG, ECG_start_offset, ECG_sample_start_end_idx, ECG_header, ECG_annotations, ECG_annotations_start_end_idx, payload_in )

            obj.progress_handle.hide()
            
            payload = [];
            
            aux_val = ECG_sample_start_end_idx + ECG_start_offset - 1;
            disp_string_title(1, sprintf( 'Correcting from %d to %d', aux_val(1), aux_val(2) ) );
            
            if( isempty(payload_in) )
                Ann_struct = ECG_annotations;
            else
                Ann_struct = payload_in;
            end
            
            QRScorrector('ECG', ECG, 'ECG_header', ECG_header, 'QRS_annotations', Ann_struct, 'OutputVarName', 'payload', 'Figure', figure(obj.fig_hdl) );

            disp_string_framed(1, 'User interaction required' );
            
            fprintf(1, 'This ECGtask allow user interaction. Press CTRL + righ arrow (->) to save results and press F5 (Run) to continue.\n')
            keyboard
            
            if( ishandle(obj.fig_hdl) && (isempty(payload) || ~isstruct(payload)) )
                disp_string_framed(2, 'Payload variable not saved' );
                fprintf(1, 'Press CTRL + righ arrow (->) to save results and press F5 (Run) to continue.\n')
                keyboard
            end

            close(obj.fig_hdl)
            
            obj.progress_handle.show()
            
        end
        
        function payload = Finish(obj, payload, ECG_header)
           
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

        
    end
    
    methods ( Access = private )
        
        
    end
    
end
