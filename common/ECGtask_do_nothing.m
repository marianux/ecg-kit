classdef ECGtask_do_nothing < ECGtask

% Null ECGtask (for Matlab)
% ---------------------------------
% 
% Description:
% 
% Abstract class for defining ECGtask interface
% 
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate  : 20/2/2013
% Last update: 20/2/2013
       
    properties(GetAccess = public, Constant)
        name = 'Null task';
        target_units = 'uV';
        doPayload = false;
    end
    
    properties(GetAccess = public, SetAccess = private)
        % if user = memory;
        % memory_constant is the fraction respect to user.MaxPossibleArrayBytes
        % which determines the maximum input data size.
        memory_constant = 0.4;
    end
    
    properties
        progress_handle
        user_string
    end

    methods
        
        function obj = ECGtask_do_nothing(obj)
            
        end
        
        function Start(obj, ECG_header)
            % not implemented

        end
        
        function payload = Process(obj, ECG, ECG_sample_start_end_idx, ECG_header, ECG_annotations, ECG_annotations_start_end_idx )
            
            % not implemented
            payload = [];
            
        end
        
        
        function Finish(obj)
            % not implemented

        end
        
        function payload = Concatenate(obj, plA, plB)
            
            % not implemented
            payload = [];
            
        end

    end
    
end
