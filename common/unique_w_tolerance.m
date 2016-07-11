%% (Internal) Find unique elements in a vector within a win_size tolerance window
%
%     [unique_elements idx ] = unique_w_tolerance( val1, win_size )
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
% Author: Mariano Llamedo Soria (llamedom at frba.utn.edu.ar)
% Version: 0.1 beta
% Birthdate: 14/05/2016
% Last update: 14/05/2016
% Copyright 2008-2016
% 
function [unique_elements, idx ] = unique_w_tolerance( val1, win_size )

if( nargin < 1 || isempty(val1) )
    unique_elements = [];
    idx = [];
    return;
end

if( nargin < 2 || isempty(win_size) )
    win_size = 1;
end

[val1, idx] = unique(val1);

if( length( val1) == 1 )
    unique_elements = val1;
    return;
end

bAux = true;
while( any(bAux) )
    bAux = abs(val1(1:end-1) - val1(2:end)) <= win_size;
    bAux = [bAux; bAux(end)];
    % keep first occurrence
    % avoid removing consecutive occurrences, keep only the first.
    % Use a logic function to do this.
    bAux1 = toeplitz(double([ false(3,1); bAux ]));
    bAux1 = logical(bAux1(3:end-1,1:4));
    bAux1 = bAux1(:,1) & (~bAux1(:,2)) | bAux1(:,1) & bAux1(:,3) & (~bAux1(:,4));
    
    val1(bAux1) = [];
    idx(bAux1) = [];

end

unique_elements = val1;

