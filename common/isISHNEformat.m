%% (Internal) Check if a recording is in ISHNE format.
%   
%   bRetval = isISHNEformat(filename)
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
% See also ECGformat, read_ECG, isHL7aformat
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Birthdate: 18/2/2013
% Last update: 18/2/2013
% Copyright 2008-2015
% 
function bRetval = isISHNEformat(filename)

bRetval = false;

fid = fopen(filename);

if( fid > 0 )

    try 
        
        magic_num = fread(fid, 8,'*char');

        if( feof(fid) || ~strcmpi(rowvec(magic_num(1:5)), 'ISHNE') )
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
