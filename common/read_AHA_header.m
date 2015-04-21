%% Reads ECG header in AHA format
% Reads ECG annotations in AHA format. Implements the documentation
% available in the help of the application provided with the database.
% 
% Reads ECG headers in AHA format. Implements the documentation
% available in the help of the application provided with the database.
% 
% Arguments:
%   + filename: recording to extract header block.
% 
% Output:
%   + heasig: header with the ECG properties. 
% 
% 
% See also read_AHA_format, read_AHA_ann, read_ECG, ECGwrapper
% 
% Author: Mariano Llamedo Soria
% <matlab:web('mailto:llamedom@electron.frba.utn.edu.ar','-browser') (email)> 
% Version: 0.1 beta
% Birthdate: 17/12/2010
% Last update: 19/11/2014
% Copyright 2008-2015
% 
function heasig = read_AHA_header(filename)

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
    % hardcoded fields
    [~, heasig.recname ] = fileparts(filename);
    heasig.btime = '00:00:00';
    heasig.bdate = '01/01/2000';
    heasig.desc = char(strcat(repmat({'ECG'}, heasig.nsig,1), cellstr(num2str(colvec(1:heasig.nsig))) ));
    
end
