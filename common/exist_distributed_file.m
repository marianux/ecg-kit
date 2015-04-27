%% (Internal) Check the existence of a file in a distributed (slow) filesystem
%   
%   bFileFound = exist_distributed_file(file_name, retries)
% 
% Arguments:
% 
%      + file_name: 
% 
%      + retries: times to check the existence
%             
% Output:
% 
%      + bFileFound:
% 
% Example:
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function bFileFound = exist_distributed_file(file_name, retries)

bFileFound = false;
while( ~bFileFound && retries >= 0)
    if( exist(file_name, 'file') )
        bFileFound = true;
    else            
       retries = retries - 1;
       pause(10)
    end
end
