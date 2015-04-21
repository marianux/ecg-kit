%% (Internal) Add the path only if not was already added.
%   
%   paths_added = addpath_if_not_added(str_path, position)
% 
% Arguments:
% 
%      + str_path: paths to add
% 
%      + position: where to add. See addpath doc.
% 
% Output:
% 
%      + paths_added: Only the paths added
% 
% Example:
% 
% See also addpath
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function paths_added = addpath_if_not_added(str_path, position)

paths_added = [];

if( nargin < 2 || isempty(position) )
    position = '-end';
end

if( ischar(str_path) )
    str_path = cellstr(str_path);
end

all_path = path;
for each_str = rowvec(str_path)
    aux_str = each_str{1};
    aux_val = strfind(all_path, aux_str);
    if(isempty(aux_val))
        paths_added = [paths_added; each_str];
        addpath(aux_str, position);
    end
end
