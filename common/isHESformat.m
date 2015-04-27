%% (Internal) Check if a recording is in HES format.
%   
%   bRetval = isHESformat(filename)
% 
% Arguments:
% 
%      + filename: the recording
% 
% Output:
% 
%      + bRetval: Boolean if it is of this format.
% 
% Example:
% 
% See also ECGformat, read_ECG, isISHNEformat
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Birthdate: 18/2/2013
% Last update: 18/2/2013
% Copyright 2008-2015
% 
function bRetval = isHESformat(filename)

bRetval = false;

fid = fopen(filename);

if( fid > 0 )

    try 
        
        fseek(fid, 894, 'bof');
        
        magic_num = rowvec(fread(fid, 6,'*char'));

        if( feof(fid) || ~strcmpi(magic_num, 'zHYeWs') )
            fclose(fid);
            return
        end
        
        fclose(fid);

        bRetval = true;
        
    catch ME
        
        fclose(fid);
        
        rethrow(ME)
    end
    
end
