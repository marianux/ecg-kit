%% (Internal) Intersection of two sets with tolerance
%
%     [sft_intersect1 idx1 idx2 ] = soft_intersect( val1, val2, win_size )
% 
% 
% Arguments:
% 
%   + val1, val2: data elements to be intersected, 
% 
%   + win_size: tolerance to consider val1(i) == val2(i). In fact we
%   consider equal two elements if val1(val1 >= (val2(i) - win_size) & val1 <= (val2(i) + win_size) )
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
% See also soft_set_difference
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate: 17/12/2010
% Last update: 17/12/2010
% Copyright 2008-2015
% 
function [sft_intersect1 idx1 idx2 ] = soft_intersect( val1, val2, win_size )

% Find the intersection between index sequence val1-2 within a windows
% win_size.
idx1 = [];
idx2 = [];

sft_intersect1 = unique(cell2mat(arrayfun(@(a)( colvec( val1(val1 >= (a - win_size) & val1 <= (a + win_size) ) ) ), colvec(val2), 'UniformOutput', false)));

if( nargout > 1)
    [ ~, idx1] = intersect(val1, sft_intersect1);
end

if( nargout > 2)
    sft_intersect2 = unique(cell2mat(arrayfun(@(a)( colvec( val2(val2 >= (a - win_size) & val2 <= (a + win_size) ) ) ), colvec(sft_intersect1), 'UniformOutput', false)));
    [ ~, idx2] = intersect(val2, sft_intersect2);
end
