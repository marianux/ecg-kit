%% (Internal) Filtered calculation of RR time series from QRS detections
% 
% This function calculates a gap-removed version of the RR interval
% sequence.
% 
%       RR = RR_calculation(QRS_detections, sampling_rate)
% 
% Arguments:
% 
%	   QRS_detections: the QRS detections 
% 
%	   sampling_rate: the sampling frequency of the detections.
% 
% Example
% 
% See also resample_sequences, MedianFiltSequence
% 
% Author: Mariano Llamedo Soria (llamedom at frba.utn.edu.ar)
% Version: 0.1 beta
% Birthdate  : 25/4/2017
% Last update: 25/4/2017
% Copyright 2008-2017
% 
function [ RR, RR_filt] = RR_calculation(QRS_detections, sampling_rate, RR_filt)

    RR = diff(QRS_detections);
    RR = [ RR(1,:); RR ];

    if( nargin < 3 || length(RR_filt) ~= length(QRS_detections) )
        % resample the RR sequence. Heavy in computation
        RR_filt = MedianFiltSequence(QRS_detections, RR, 2*sampling_rate);
    end
    
    % remove gaps using next heartbeats
    gap_relative_time = 3; % times the median RR interval
    aux_val = RR ./ RR_filt;
    aux_idx = find(aux_val >= gap_relative_time);
    
    for ii = rowvec(aux_idx)
        
        aux_idx = find( QRS_detections > QRS_detections(ii) & aux_val < gap_relative_time, 1, 'first');
        
        if ( isempty(aux_idx) )
            aux_idx = find( QRS_detections < QRS_detections(ii) & aux_val < gap_relative_time, 1, 'last');
        end
        
        RR(ii) = RR( aux_idx );
        
    end
    
    
