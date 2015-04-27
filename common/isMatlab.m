%% (Internal) Check if the kit is running on Matlab
%   
%   bAux = isMatlab()
% 
% Output:
% 
%      + bAux: Boolean if it is Matlab.
% 
% Example:
% 
% See also isOctave
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Birthdate: 18/2/2013
% Last update: 18/2/2013
% Copyright 2008-2015
%
function bAux = isMatlab()

bAux = false;

matlab_ver = ver('Matlab');

if( ~isempty(matlab_ver) && strcmpi(matlab_ver.Name, 'MATLAB') )
    bAux = true;
end
