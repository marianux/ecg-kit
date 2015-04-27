%% (Internal) Logit function
%   
%   y = logit_function( x )
% 
% Arguments:
% 
%             
% Output:
% 
% 
% Example:
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function y = logit_function( x )

y = 1./(1+exp(-x));
