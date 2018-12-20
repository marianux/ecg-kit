%% Creates set of QRS detections, based on multi-lead-and-algorithm QRS detections.
% This is an auxiliar function for the ECGtask_QRS_detections_post_process.
% *mixartif* use the output of *all* leads and detectors to perform a
% multi-lead-and-algorithm composition. The result of this algorithm is a
% new set of detections based on the concatenation of the ""best""
% detections found for each 20-seconds window in a recording. So, this
% algorithm generates *new* QRS detection series, as described in [add
% reference]. 
% 
%   all_detections = mixartif(struct_in, ECG_header, start_end_this_segment )
% 
% Arguments:
% 
%      + struct_in: the structure which results from loading the result of
%                   the QRS detection task :
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
%      + all_detections: struct with the artificial detections named with
%                        prefix "wavedetMix_ECGmix". 
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
% See also ECGtask_QRS_detections_post_process, wavedetmix, best_m_lead, qrs_detection_and_correction
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Birthdate: 01/01/2014
% Last update: 14/07/2015
% Copyright 2008-2015
% 
function all_detections = mixartif(struct_in, ECG_header, start_end_this_segment)

    all_detections = [];

    if( nargin < 3 || isempty(start_end_this_segment) ) 
        start_end_this_segment = [1 ECG_header.nsamp];
    end
    % attemp to build a better detection from single-lead detections.

    [~, all_annotations] = getAnnNames(struct_in);

    [ ratios, estimated_labs ] = CalcRRserieRatio(all_annotations, ECG_header, start_end_this_segment);

    [~, best_detections_idx] = sort(ratios, 'descend');

    % generate artificial annotations combining K best annotations
    aux_idx = best_detections_idx(1:min(10, length(best_detections_idx)) );
    artificial_annotations = combine_anns(all_annotations(aux_idx), estimated_labs(aux_idx), ECG_header );

    for ii = 1:length(artificial_annotations)
        aux_str = ['mixartif_ECGmix' num2str(ii)];
        all_detections.(aux_str) = artificial_annotations(ii);
    end



