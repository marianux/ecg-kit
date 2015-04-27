%% (Internal) Detect zero-crosses in a signal
%   
% Returns the index of the input vector in which the first zero crossing is located.
% 
%   index = myzerocros(x)
% 
% Arguments:
% 
%      + x: signal
% 
% Output:
% 
%      + index: indexes where zero-crosses occurs.
% 
% Example:
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function index = myzerocros(x)

sx = sign(x);
diffsx = diff(sx);
diffsx = [diffsx(1); diffsx];
index = find( diffsx ~= 0);

