%% (Internal) Find unique elements in a vector within a win_size tolerance window
%
%     [unique_elements idx ] = alineate_positions( val1, win_size )
% 
% 
% Arguments:
% 
%   + val1: data elements, 
% 
%   + win_size: tolerance to consider val1(i) == val1(i+1).
% 
% Output:
% 
%   + sft_intersect1: the elements in the soft intersection
% 
%   + idx1, idx2: the indexes of val1(idx1) which are in the soft
%   intersection.
% 
% Example:
% 
% 
% See also soft_set_intersect
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate: 14/05/2016
% Last update: 14/05/2016
% Copyright 2008-2016
% 
function [unique_elements, idx ] = alineate_positions( val1, win_size )

aux_mat = [ [aux_on;nan(cant_anns - length(aux_on),1)] [aux_peak;nan(cant_anns - length(aux_peak),1)] [aux_off;nan(cant_anns - length(aux_off),1)] ];

aux_dist = cat(3, abs( repmat(aux_mat(:,1),1,cant_anns) - repmat(aux_mat(:,2),1,cant_anns)' ), abs( repmat(aux_mat(:,1),1,cant_anns) - repmat(aux_mat(:,3),1,cant_anns)' ), abs( repmat(aux_mat(:,2),1,cant_anns) - repmat(aux_mat(:,3),1,cant_anns)' ) );

sort(aux_)
