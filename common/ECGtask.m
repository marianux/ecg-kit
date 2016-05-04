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
% Copyright 2008-2015
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
        tmp_path
    end

    methods (Abstract)
                
% ECG_header is the header for this signal
% ECG_annotations are the QRS locations available for this signal
        Start(obj, ECG_header, ECG_annotations)
        
% ECG is the ECG signal
% ECG_start_offset is the location of ECG(1,:) within the whole signal
% ECG_sample_start_end_idx are the start and end samples within ECG to
%    generate valid results.
% ECG_header is the header for this signal
% ECG_annotations are the QRS locations available for this signal
% ECG_annotations_start_end_idx are the start and end indexes corresponding
%    to the first and last element of ECG_annotations in the current
%    iteration.
        payload = Process(ECG, ECG_start_offset, ECG_sample_start_end_idx, ECG_header, ECG_annotations, ECG_annotations_start_end_idx )

% payload is the result of the task, after the final processing done by this method        
% ECG_header is the header for this signal
        payload = Finish(obj, payload, ECG_header)
        
        
% plA, plB are two payloads produced by the Process method.
        payload = Concatenate(plA, plB)

    end
    
end
