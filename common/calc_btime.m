%% (Internal) Creates a string with the base time
%   
%   btime = calc_btime( start_sample, sampling_rate )
% 
% Arguments:
% 
%      + start_sample: 
% 
%      + sampling_rate: 
% 
% Output:
% 
%      + btime: a string HH:MM:SS:mm
% 
% Example:
% 
% See also read_mortara
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function btime = calc_btime( start_sample, sampling_rate )

hours = floor((start_sample-1)/sampling_rate/60/60);
mins = floor((start_sample-1)/sampling_rate/60 - hours *60);
secs = floor((start_sample-1)/sampling_rate - mins *60 - hours * 60 * 60);
milli = round(((start_sample-1)/sampling_rate - mins *60 - hours * 60 * 60 - secs) * 1000);
btime = sprintf('%0d:%0d:%0d:%03d',hours, mins, secs, milli);

