%% (Internal) Reshape a matrix into a row vector
% Return x as a row vector. This function is useful when a function returns a
% column vector or matrix and you want to immediately reshape it in a functional
% way. Suppose f(a,b,c) returns a column vector, matlab will not let you write
% f(a,b,c)(:)' - you would have to first store the result. With this function you
% can write rowvec(f(a,b,c)) and be assured that the result is a row vector.   
% 
%    x = colvec(x)
% 
% Arguments:
% 
%      + x: data
% 
% Output:
% 
%      + x: columnized data
% 
% Example:
% 
% See also colvec
% 
% This file is from pmtk3.googlecode.com
% Copyright 2008-2015
% 
function x = rowvec(x)
    x = x(:)';
end
