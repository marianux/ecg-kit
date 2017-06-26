%% (Internal) Merge QRS detections into a single structure
% 
% This function merge together two QRS detection structures.
% 
%       detC = merge_QRS_detections(detA, detB)
% 
% Arguments:
% 
%	   detA, detB: the detections to merge
% 
%	   detC: merged detections
% 
% Example
% 
% See also ECGtask_QRS_detection
% 
% Author: Mariano Llamedo Soria (llamedom at frba.utn.edu.ar)
% Version: 0.1 beta
% Birthdate  : 21/4/2017
% Last update: 21/4/2017
% Copyright 2008-2017
% 
function detC = merge_QRS_detections(detA, detB, ECG_header)

    if( isstruct(detA) )
        
        bStruct = true;
        
        if( isfield(detA, 'series_quality') && isfield(detA.series_quality, 'sampfreq' ) )
            detA_sampfreq = detA.series_quality.sampfreq;
        else
            % assumed the same sampling rate
            detA_sampfreq = [];
        end
        
        detection_namesA = setdiff(fieldnames(detA), {'series_quality' 'series_performance'} );

        aux_val = {};
        for fname = rowvec(detection_namesA)
            aux_val = [aux_val; {colvec(detA.(fname{1}).time)} ];
        end
    else
        bStruct = false;
        
        aux_val = {detA};
    end
    
    if( isstruct(detB) )

        if( isfield(detB, 'series_quality') && isfield(detB.series_quality, 'sampfreq' ) )
            detB_sampfreq = detB.series_quality.sampfreq;
        else
            % assumed the same sampling rate
            detB_sampfreq = [];
        end

        if( ~isempty(detA_sampfreq) &&  ~isempty(detB_sampfreq) )
            % ratio to convert detB @ detB_sampfreq to detA @ detA_sampfreq
            aux_sampfreq_ratio = detA_sampfreq / detB_sampfreq ;
        else
            % assume the same samp freq
            aux_sampfreq_ratio = 1;
        end
            
        detection_namesB = setdiff(fieldnames(detB), {'series_quality' 'series_performance'} );

        for fname = rowvec(detection_namesB)
            aux_val = [aux_val; { unique( round( colvec(detB.(fname{1}).time) * aux_sampfreq_ratio )) } ];
        end

    else
        aux_val = [aux_val; {detB} ];
    end
    
    merged_detections = aux_val{1};
    
    for ii = 2:length(aux_val)
        
        merged_detections = soft_set_union(merged_detections, aux_val{ii}, round(0.1*ECG_header.freq) );
        
    end
    
    if( bStruct )
        % build the detection struct
        detC.merged.time = merged_detections;

        detC = calculateSeriesQuality(detC, ECG_header, [1 max(merged_detections)+1 ] );
    else
        detC = merged_detections;
    end
    