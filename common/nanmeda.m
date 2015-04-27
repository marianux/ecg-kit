%% (Internal) Calculate the median of absolute deviations from the median (MEDA)
%   
%   meda = nanmeda(data)
% 
% Arguments:
% 
%      + data: data
% 
% Output:
% 
%      + meda:  median of absolute deviations from the median.
% 
% Example:
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function meda = nanmeda(data)

median_data = nanmedian(data);
meda = nanmedian(abs(bsxfun(@minus, data, median_data)));
