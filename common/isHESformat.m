function bRetval = isHESformat(filename)

% Check if a recording is in HES format.
% 
% Author: Mariano Llamedo Soria (llamedom at {unizar.es;electron.frba.utn.edu.ar}
% Birthdate: 18/2/2013
% Last update: 18/2/2013
% 
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
