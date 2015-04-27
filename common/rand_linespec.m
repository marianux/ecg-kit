%% (Internal) Example of user-created QRS detector
%
%           line_spec = rand_linespec()
% 
% Output:
% 
%      + line_spec : the string with a random linespec for using with
%      plotting functions
% 
% Example:
% 
% 
% See also plot
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate  : 30/7/2014
% Last update: 30/7/2014
% Copyright 2008-2015
% 
function line_spec = rand_linespec()

mrk={'+','o','*','.','x','s','d','^','v','<','>','p','h'}.';
linestyle = {'-','--',':','-.'};
colspec={'m','c','r','g','b','k'}.';

line_spec = [linestyle{randsample(1:length(linestyle),1)} colspec{randsample(1:length(colspec),1)} mrk{randsample(1:length(mrk),1)} ];
