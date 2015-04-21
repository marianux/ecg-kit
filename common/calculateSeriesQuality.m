%% Estimate the quality of QRS complex detections 
% Calculate the quality of QRS detections based on the the likelihood of
% measurements from RR interval series with respect to a trained model. See
% for further details. See 'A Pattern-Recognition Approach for
% Lead-Selection in Heartbeat Detection' CINC 2014
%   
% Example
% 
%   payload_out = calculateSeriesQuality(payload_out, ECG_header, start_end_this_segment)
% 
% Arguments:
%   + payload_out: structure with the QRS detections series to estimate the
%       quality. One fieldname contains the QRS detections series in a
%       'time' field. For example
%       payload_out.('QRS_detector1_name_lead_I').time contains the
%       detections of QRS_detector1 in lead1. See ECGtask_QRS_detection for
%       examples.%       
%   + ECG_header is a structure (ECG_header prop. in ECGwrapper object)
%       describing the signal.
%   + start_end_this_segment: array samples where starts and finish this
%       signal segment [start_sample end_sample].
% 
% Output:
%   + payload_out: The same input structure with extra field
%       'series_quality' describing the quaility of the QRS detection
%       series in three fields: ratios, a measure of quality of the QRS
%       detections, estimated_labs, for internal use, and AnnNames which
%       matches the previous fields to the fieldname found in payload_out.
% 
% See also ECGtask_QRS_detection
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Birthdate: 01/01/2012
% Last update: 18/10/2014
% Copyright 2008-2015
function payload_out = calculateSeriesQuality(payload_out, ECG_header, start_end_this_segment)

    % estimate quality of QRS detections performed
    [AnnNames, all_annotations] = getAnnNames(payload_out);

    [ ratios, estimated_labs ] = CalcRRserieRatio(all_annotations, ECG_header, start_end_this_segment);
    [ratios, best_detections_idx] = sort(ratios, 'descend');
    AnnNames = AnnNames(best_detections_idx,:);

    payload_out.series_quality.ratios = ratios;
    payload_out.series_quality.estimated_labs = estimated_labs(best_detections_idx);
    payload_out.series_quality.AnnNames = AnnNames;
