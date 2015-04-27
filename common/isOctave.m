%% (Internal) Check if the kit is running on Octave
%   
%   bAux = isOctave()
% 
% Output:
% 
%      + bAux: Boolean if it is Octave.
% 
% Example:
% 
% See also isMatlab
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Birthdate: 18/2/2013
% Last update: 18/2/2013
% Copyright 2008-2015
%
function bAux = isOctave()

bAux = false;

octave_ver = ver('Octave');

if( ~isempty(octave_ver) && strcmpi(octave_ver.Name, 'Octave') )
    bAux = true;
end
