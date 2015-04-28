%% (Internal) Set difference with tolerance
%
%     [sft_set_difference idx]= soft_set_difference( val1, val2, win_size )
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
%   + sft_set_difference: the elements which are in val1, but not in soft_intersect( val1, val2, win_size )
% 
%   + idx: the indexes of val1(idx) which are in the sft_set_difference
% 
% Example:
% 
% 
% See also soft_intersect
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate: 17/12/2010
% Last update: 17/12/2010
% Copyright 2008-2015
% 
function [sft_set_difference, idx]= soft_set_difference( val1, val2, win_size )

% Find the set difference between index sequence val1-2 within a windows
% win_size.

sft_intersect = soft_intersect( val1, val2, win_size );

[sft_set_difference, idx] = setdiff(val1, sft_intersect);
