%% (Internal) Check if the version a is later than version b
%   
%   bAux = islater( ver_a, ver_b )
% 
% Output:
% 
%      + bAux: Boolean if it is later.
% 
% Example:
% 
% See also isMatlab, isOctave
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Birthdate: 18/2/2013
% Last update: 18/2/2013
% Copyright 2008-2015
%
function bAux = islater( ver_a, ver_b )

    bAux = false;
    
    aux_tokens_a = regexp(ver_a, 'v(\d+)\.(\d+)\.(\d+).*', 'tokens');
    aux_tokens_b = regexp(ver_b, 'v(\d+)\.(\d+)\.(\d+).*', 'tokens');
    
    if( ~isempty(aux_tokens_a) && ~isempty(aux_tokens_b) )
    
        va = str2double( aux_tokens_a{1});
        vb = str2double( aux_tokens_b{1});

        bAux = diff([va;vb] * [1e6;1e3;1]) < 0;
        
    end
        
end

