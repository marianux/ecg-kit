%% (Internal) Serialize in columns data
% Reshape the input into a column vector
% This function is useful when a function returns a row vector or matrix
% and you want to immediately reshape it in a functional way. Suppose
% f(a,b,c) returns a row vector, matlab will not let you write f(a,b,c)(:)
% - you would have to first store the result. With this function you can
% write colvec(f(a,b,c)) and be assured that the result is a column vector.

% This file is from pmtk3.googlecode.com
% 
%   x = colvec(x)
% 
% Arguments:
% 
%      + x: data
% 
% Output:
% 
%      + x: vertically serialized data.
% 
% Example:
% 
% See also rowvec
% 
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function x = colvec(x)


x = x(:);
end
