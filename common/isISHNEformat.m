function bRetval = isISHNEformat(filename)

% Check if a recording is in ISHNE format.
% 
% Author: Mariano Llamedo Soria (llamedom at {unizar.es;electron.frba.utn.edu.ar}
% Birthdate: 18/2/2013
% Last update: 18/2/2013
% 
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
