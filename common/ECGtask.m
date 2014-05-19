classdef ECGtask < handle

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
       
    properties(Abstract = true, GetAccess = public, Constant)
        name
        target_units
        doPayload
        % if user = memory;
        % memory_constant is the fraction respect to user.MaxPossibleArrayBytes
        % which determines the maximum input data size.
        
    end

    properties(Abstract = true, GetAccess = public, SetAccess = private)
        memory_constant
    end
    
    properties(Abstract = true)
        progress_handle
        user_string
    end

    methods (Abstract)
                
        Start(obj, ECG_header, ECG_annotations)
        
        payload = Process(ECG, ECG_sample_start_end_idx, ECG_header, ECG_annotations, ECG_annotations_start_end_idx )

        Finish(obj)
        
        payload = Concatenate(plA, plB)

    end
    
end
