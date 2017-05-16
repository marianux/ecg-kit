%% (Internal) Filter and unevenly spaced sequence
% 
% This function calculates a filtered version of an unevenly sampled
% sequence, like the interbeat period sequence (RR sequence).
% 
%       RR_filt = MedianFiltSequence(x, sampling_rate)
% 
% Arguments:
% 
%	   x: the QRS detections 
% 
%	   sampling_rate: the sampling frequency of the detections.
% 
% Example
% 
% See also resample_sequences, RR_calculation
% 
% Author: Mariano Llamedo Soria (llamedom at frba.utn.edu.ar)
% Version: 0.1 beta
% Birthdate  : 5/5/2017
% Last update: 5/5/2017
% Copyright 2008-2017
% 
function y_filt = MedianFiltSequence(x, y, filter_win)

% oversample the original sequence
oversampling_ratio = round(median(diff(x))/4);

x_interp = colvec(1:oversampling_ratio:max(x));
y_interp = spline(x, y, x_interp);

y_filt = MedianFilt(y_interp, round(filter_win/oversampling_ratio));

aux_idx = ceil( x / oversampling_ratio );
y_filt = y_filt(aux_idx);
