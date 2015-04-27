%% (Internal) Read the MIT 310 format
%
%     [ECG samples_read] = read_310_format( filename, start_sample, end_sample, heasig )
% 
% Reads ECG recordings in 310 format. Implements the documentation
% available in WFDB.
% 
% Arguments:
%   + filename: recording to be read.
%   + start_sample: (opt) start sample to read. Default 1.
%   + end_sample: (opt) end sample to read. Default min(All recording, ECG block of 200 Mbytes)
% 
% Output:
%   + ECG: the ECG block
%   + heasig: header with the ECG properties. 
%   + ann: annotations for the ECG recordings.
% 
% Limits:
% This routine is limited to read blocks smaller than 200 Mbytes for
% performance reasons. You can disable this limit by doing:
% MaxIOread = Inf; %megabytes
% 
% 
% Example:
% 
% 
% See also read_ECG
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate: 17/12/2010
% Last update: 17/12/2010
% Copyright 2008-2015
% 
function [ECG samples_read] = read_310_format( filename, start_sample, end_sample, heasig )

samples_read = [];
ECG = [];

%No leer bloques mas grandes de 200 megabytes
MaxIOread = 200; %megabytes

if( nargin < 4 || isempty( heasig ) )
    %intento leer heasig.
    heasig = [];
end

if( nargin < 2 || isempty( start_sample ) )
    start_sample = 1;
else
    start_sample = max(1,start_sample);
end

if( heasig.nsig ~= 3 )
    error('Considerado solo para 3 señales.');
end

fidECG = fopen( filename, 'r');

if( fidECG > 0 )

    if( nargin < 3 || isempty( end_sample ) )
        %Intento la lectura total por defecto
        samples2read = heasig.nsamp - (start_sample-1);
    else
        samples2read = min(heasig.nsamp, end_sample) - (start_sample-1);
    end
    
    if( (samples2read*heasig.nsig*2) > (MaxIOread * 1024^2) )
        samples2read = (MaxIOread * 1024^2) / heasig.nsig / 2;
        warning(['No es recomendable leer mas de ' num2str(MaxIOread) ' Mb. Realice varias lecturas.'])
    end
    
    try 
        fseek(fidECG, start_sample-1, 'bof');
        aux = fread(fidECG, [2 samples2read], '*uint16')';
        ECG(:,1) = bitand( aux(:,1), uint16(2046) );
        ECG(:,2) = bitand( aux(:,2), uint16(2046) );
        ECG(:,3) = bitshift(bitand( aux(:,1), uint16(63488) ),-10) + bitshift(bitand( aux(:,2), uint16(63488) ),-5);
        ECG( ECG > 1024 ) = ECG( ECG > 1024 ) - 2048;
        fclose(fidECG);
    catch ME
        fclose(fidECG);
        rethrow(ME)
    end

    samples_read = size(ECG,1);
    
end
