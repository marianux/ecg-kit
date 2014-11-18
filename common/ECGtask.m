%% Defines the class interface for the ECGtask derived classes
% 
% Description:
% Abstract class for defining ECGtask interface
% 
% Arguments: % 
% 
% Output:
% 
% Examples:
% 
% 
% See also ECGtask_QRS_detection, ECGtask_ECG_delineation
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 20/3/2012
% Copyright 2008-2014
% 
classdef ECGtask < handle
       
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
        started
    end
    
    properties(Abstract = true)
        progress_handle
        user_string
        tmp_path
    end

    methods (Abstract)
                
        Start(obj, ECG_header, ECG_annotations)
        
        payload = Process(ECG, ECG_sample_start_end_idx, ECG_header, ECG_annotations, ECG_annotations_start_end_idx )

        Finish(obj)
        
        payload = Concatenate(plA, plB)

    end
    
end
