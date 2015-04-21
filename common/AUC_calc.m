%% Compute area under the ROC curve (AUC).
% Compute area under the ROC curve.
%   
% Example
% 
%   AUC_calc(prroc_new(dsResult1))
% 
% Arguments:
% 
% Output:
%   The AUC  
% 
% See also prroc_new
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Birthdate: 22/10/2014
% Last update: 22/10/2014
% Copyright 2008-2015
% 
function AUC = AUC_calc(aux_roc)

    delta_x = diff(aux_roc.xvalues);
    delta_y = diff(1-aux_roc.error);
    AUC = delta_x * (1-aux_roc.error(1:end-1))' + (delta_x * delta_y')/2;
