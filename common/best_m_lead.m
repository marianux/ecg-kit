%% Filter the best QRS detections according to the *m* quality criterion
% This is an auxiliar function for the ECGtask_QRS_detections_post_process.
% Lead selection strategy based on the best (quality metric) *m*.
% This is a lead selection algorithm, so this algorithm filters the
% first ranked lead of all the algorithms outputs included in the
% payload struct. In contrast with other postprocessing algorithm, this
% doesn't create a new detection serie, it is just a filter.
% 
%   all_detections = best_m_lead(struct_in, ECG_header, start_end_this_segment)
% 
% Arguments:
% 
%      + struct_in: the structure which results from loading the result of
%                   the QRS detection task:
% 
%                   cached_filenames = ECG_w.GetCahchedFileName('QRS_detection');
%                   struct_in = load(cached_filenames{1});
% 
%      +ECG_header: [struct] OPTIONAL. 
% 
%             Description of the ECG typically available in the
%             ECG_header. Structure with fields:
% 
%               -freq: Sampling rate in Hz. (1)
% 
%               -nsig: Number of ECG leads. (size(ECG,2))
% 
%               -nsamp: Number of ECG samples. (size(ECG,1))
% 
%      + start_end_this_segment: an array with the first and last sample
%                                indexes.
% 
% Output:
% 
%      + all_detections: struct with the detections  of the best lead
%                        according to the *m* quality metric.
%                        
% 
% Example
% 
%     ECG_w.ECGtaskHandle = 'QRS_detections_post_process';
%     % Mix QRS detections strategy function
%     ECG_w.ECGtaskHandle.post_proc_func = 'calculate_artificial_QRS_detections';
%     ECG_w.ECGtaskHandle.payload = load(cached_filenames{1});
%     ECG_w.ECGtaskHandle.CalculatePerformance = true;
%     ECG_w.Run;
% 
% See also ECGtask_QRS_detections_post_process, mixartif, wavedetmix, qrs_detection_and_correction
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Birthdate: 12/05/2018
% Last update: 12/05/2018
% Copyright 2008-2018
% 
function all_detections = best_m_lead(struct_in, ECG_header, start_end_this_segment)

    all_detections = [];

    if( nargin < 3 || isempty(start_end_this_segment) ) 
        start_end_this_segment = [1 ECG_header.nsamp];
    end

    if( isfield(struct_in, 'series_quality') )
        % choose the best ranked automatic detection
        best_detections_idx = max_index(struct_in.series_quality.ratios);
        best_det_name = struct_in.series_quality.AnnNames{best_detections_idx,1};
        all_detections.(best_det_name) = struct_in.(best_det_name);
        
    else
        disp_string_framed(2, 'Could not identify QRS detections in the input payload.');
    end
    
  