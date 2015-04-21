%% (Internal) Generate a PIDs work list
%   
% 
%   [pid_starts, pid_ends] = TaskPartition( task_size, cant_pid)
% 
% Arguments:
% 
%      + task_size: A positive integer with the size of the job
% 
%      + cant_pid: The amount of PIDs to work
% 
% Output:
% 
%      + pid_starts: An index array of size cant_pid x 1 with the starting
%      indexes for each PID
% 
%      + pid_ends: An index array of size cant_pid x 1 with the ending
%      indexes for each PID
% 
% Example:
% 
%       [pid_starts, pid_ends] = TaskPartition( 10, 2)
% 
%       pid_starts = [ 1 6 ] 
%       pid_ends =   [ 5 10 ] 
% 
% See also ECGwrapper
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function [pid_starts, pid_ends] = TaskPartition( task_size, cant_pid)

%no es recomendable hacerlo mas grande de cant_recs la particion.
cant_pid = min(cant_pid, task_size);

things2do = fix(task_size / cant_pid);

remainder = rem(task_size, cant_pid);

cantThingsXpid = repmat(things2do, cant_pid,1);

cantThingsXpid(1:remainder) = cantThingsXpid(1:remainder) + 1;

pid_ends = cumsum(cantThingsXpid);
pid_starts = [1;pid_ends(1:end-1)+1];
