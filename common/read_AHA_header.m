function heasig = read_AHA_header(filename)

% Reads ECG headers in AHA format. Implements the documentation
% available in the help of the application provided with the database.
% 
% Arguments:
%   + filename: recording to extract header block.
% 
% Output:
%   + heasig: header with the ECG properties. 
% 
% Limits:
% This routine identifies ECG leads as described in the documentation
% available in the THEW web page. Whereas, the author is not sure after
% reading this documentation that the leads aVL and aVR are identified
% correctly (See lead_transformation). Any feedback is very welcome.
% 
% Author: Mariano Llamedo Soria (llamedom at {unizar.es;electron.frba.utn.edu.ar}
% Birthdate: 17/12/2010
% Last update: 17/12/2010
% 

heasig = [];

%default settings for the whole DB
fidECG = fopen( filename, 'r');

if( fidECG > 0 )

    try
        fseek(fidECG, 0, 'eof');
        bytes_totales = ftell(fidECG);
        %esto es as� en la DB, si fuera necesario se podria cambiar.
        heasig.nsamp = (bytes_totales / 5);
        fclose(fidECG);
    catch ME
        fclose(fidECG);
        rethrow(ME)
    end
    
    heasig.nsig = 2;
    heasig.freq = 250;
    heasig.gain = repmat(400,1,heasig.nsig);
    heasig.adczero = zeros(1,heasig.nsig);
end
