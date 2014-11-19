function bRetval = isAHAformat(filename)

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
        fseek(fid, 0, 'eof');
        bytes_totales = ftell(fid);
        %esto es as� en la DB, si fuera necesario se podria cambiar.
        nsamp = (bytes_totales / 5) - fix(bytes_totales / 5);
        
        if( nsamp == 0 )
            
            nsamp = (bytes_totales / 5);            
            
            fseek(fid, 4, 'bof');
            
            iAnnotations = fread(fid, nsamp, '*char', 4);
            
            % two int samples (ECG1-2) and a fifth byte to code each
            % heartbeat, if not heartbeat, there is a '.' for most of the
            % samples.
            % More than half fifth samples are not heartbeats.
            if( sum(iAnnotations == '.') > 0.5*nsamp )
                bRetval = true;
            end
            
        end
        
        fclose(fid);
    catch ME
        fclose(fid);
        rethrow(ME)
    end
    
end
