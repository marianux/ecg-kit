%% (Internal) Set union with tolerance
%
%     [sft_set_union, idx1, idx2]= soft_set_union( val1, val2, win_size )
% 
% 
% Arguments:
% 
%   + val1, val2: data sets 
% 
%   + win_size: tolerance to consider val1(i) == val2(i). In fact we
%   consider equal two elements if val1(val1 >= (val2(i) - win_size) & val1 <= (val2(i) + win_size) )
% 
% Output:
% 
%   + sft_set_union: the elements which are in val1, and val2, but not
%   those that are within a win_size tolerance.
% 
%   + idx1-2: the indexes of val1(idx1) and val2(idx2) which are defined as
%   soft set union.
% 
% Example:
% 
% 
% See also soft_intersect, soft_set_difference
% 
% Author: Mariano Llamedo Soria (llamedom at frba.utn.edu.ar)
% Version: 0.1 beta
% Birthdate: 17/12/2010
% Last update: 17/12/2010
% Copyright 2008-2015
% 
function [sft_set_union, idx1, idx2]= soft_set_union( val1, val2, win_size )

% Find the set intersect between index sequence val1-2 within a windows
% win_size.

[sft_intersect1, idx1] = soft_intersect( val1, val2, win_size );

[aux_val, aux_idx] = setdiff(val1, sft_intersect1);

sft_set_union = [sft_intersect1; aux_val];
idx1 = [idx1; aux_idx];

[sft_intersect2, idx2] = soft_intersect( val2, val1, win_size );

[aux_val, aux_idx] = setdiff(val2, sft_intersect2);

sft_set_union = [sft_set_union; aux_val];
idx2 = [idx2; aux_idx];

