%% (Internal) Index of the maximum element in a vector
%   
%   max_idx = max_index(vec)
% 
% Arguments:
% 
%             
% Output:
% 
% 
% Example:
% 
% See also min_index
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function max_idx = max_index(vec)

[~, max_idx] = max(vec);
