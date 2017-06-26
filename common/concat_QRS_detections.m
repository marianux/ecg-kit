%% (Internal) Concatenate QRS detections into a single structure
% 
% This function concatenate two QRS detection structures.
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
% See also ECGtask_QRS_detection, merge_QRS_detections
% 
% Author: Mariano Llamedo Soria (llamedom at frba.utn.edu.ar)
% Version: 0.1 beta
% Birthdate  : 21/4/2017
% Last update: 21/4/2017
% Copyright 2008-2017
% 
function detC = concat_QRS_detections(detA, varargin)

    detC = detA;

    if( isfield(detC, 'series_quality') && isfield(detC.series_quality, 'sampfreq' ) )
        detC_sampfreq = detC.series_quality.sampfreq;
    else
        % assumed the same sampling rate
        detC_sampfreq = [];
    end

    for ii = 1:length(varargin)
        
        detB = varargin{ii};

        if( ~isempty(detC_sampfreq) && isfield(detB, 'series_quality') && isfield(detB.series_quality, 'sampfreq' ) )
            % ratio to convert detB @ ann_sampfreq to detC @ detC_sampfreq
            aux_sampfreq_ratio = detC_sampfreq / detB.series_quality.sampfreq;
        else
            % assumed the same sampling rate
            aux_sampfreq_ratio = 1;
        end
        
        detection_namesA = setdiff(fieldnames(detC), {'series_quality' 'series_performance'} );
        detection_namesB = setdiff(fieldnames(detB), {'series_quality' 'series_performance'} );

        [same_names, aux_idx ] = intersect(detection_namesB, detection_namesA);

        % name disambiguation
        for jj = 1:length(same_names)
            new_name = [same_names{jj} '_' ];
            detB.series_quality.AnnNames(aux_idx(jj),1) = {new_name};
            aux_val = unique( round( detB.(same_names{jj}).time * aux_sampfreq_ratio));
            detC.(new_name).time = aux_val;
        end

        diff_names = setdiff(detection_namesB, detection_namesA);

        % discrepancies
        for fname = rowvec(diff_names)
            aux_val = unique( round( detB.(fname{1}).time * aux_sampfreq_ratio ));
            detC.(fname{1}).time = aux_val;
        end

        if( isfield(detC, 'series_quality') ) 
            detC.series_quality.ratios = [ detC.series_quality.ratios ; detB.series_quality.ratios ];
            detC.series_quality.estimated_labs = [ detC.series_quality.estimated_labs ; detB.series_quality.estimated_labs ];
            detC.series_quality.AnnNames = [ detC.series_quality.AnnNames ; detB.series_quality.AnnNames ];
        end

        if( isfield(detC, 'series_performance') ) 
            detC.series_performance.conf_mat = cat(3, detC.series_performance.conf_mat, detB.series_performance.conf_mat);
            detC.series_performance.conf_mat = [ detC.series_performance.error ; detB.series_performance.error ];
        end    

    end
    