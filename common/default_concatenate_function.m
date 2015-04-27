% Description:
% 
% The default behavior of the concatenate function is to concatenate
% payloads vertically or row-wise.
% 
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate  : 17/4/2015
% Last update: 17/4/2015
% Copyright 2008-2015
% 
function payload = default_concatenate_function(plA, plB)

    if( isempty(plA) )
        payload = plB;
    else
        payload = [plA; plB];
    end
