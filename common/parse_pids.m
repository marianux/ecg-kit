%% (Internal) Identify how many PIDs are in total and which is this PID, based on a string formatted this_pid/cant_pids
%
%  [this_pid, cant_pids] = parse_pids( aux_str )
% 
% Arguments:
% 
%      + aux_str: string to parse
%             
% Output:
% 
%      + this_pid: current PID
% 
%      + cant_pids: total amount of PIDs
% 
% Example:
% 
% 
% See also ECGwrapper
% 
% Author: Mariano Llamedo Soria (llamedom at frba.utn.edu.ar)
% Version: 0.1 beta
% Birthdate  : 30/7/2014
% Last update: 30/7/2014
% Copyright 2008-2015
% 
function [this_pid, cant_pids] = parse_pids( aux_str )

if( ischar(aux_str) )
     delim = ' /-';
    [this_pid,aux_str] = strtok(aux_str,delim);
    cant_pids = strtok(aux_str, delim);
    this_pid = str2double(this_pid);
    cant_pids = str2double(cant_pids);
else
    if( length(aux_str) > 1 )
        this_pid = aux_str(1);
        cant_pids = aux_str(2);
    else
        this_pid = aux_str(1);
        cant_pids = aux_str(1);
    end
end
