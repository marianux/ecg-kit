%% (Internal) Index of the minimum element in a vector
%   
%   min_idx = min_index(vec)
% 
% Arguments:
% 
%             
% Output:
% 
% 
% Example:
% 
% See also max_index
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 17/05/2016
% Birthdate  : 17/05/2016
% Copyright 2008-2016
% 
function min_idx = min_index(vec)

[~, min_idx] = min(vec);
