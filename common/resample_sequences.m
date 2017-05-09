%% (Internal) Resample a cell array of sequences {(xx,yy)} to the same xx
%
% This function resamples the sequences included in all_seq to the same
% sync_x time in order to produce a synchronized sync_seq matrix, build
% from a resampled version of each yy sequence of values.
% 
%       [sync_seq, sync_x] = resample_sequences(all_seq)
% 
% Arguments:
% 
%       all_seq: Cell array of xx and yy values. 
% 
% Output:
% 
%       sync_seq: The resampled version of each yy included in all_seq
% 
%       sync_x: the new xx time reference common to all sync_seq values.
% 
% Example
% 
% See also RR_calculation, MedianFiltSequence
% 
% Author: Mariano Llamedo Soria (llamedom at frba.utn.edu.ar)
% Version: 0.1 beta
% Birthdate  : 9/5/2017
% Last update: 9/5/2017
% Copyright 2008-2017
% 
function [sync_seq, sync_x] = resample_sequences(all_seq)

    sync_seq = [];
    sync_x = [];    

    x_step = median(cellfun(@(a)(median(diff(a(:,1)))), all_seq));
    x_max = max(cellfun(@(a)(max(a(:,1))), all_seq));
    
    sync_x = colvec(1:x_step:x_max);
    
    sync_seq = cell2mat(rowvec((cellfun(@(a)(spline(a(:,1), a(:,2), sync_x)), all_seq, 'UniformOutput', false ))));

