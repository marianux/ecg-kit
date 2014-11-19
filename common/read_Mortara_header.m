%% Reads ECG header in Mortara format. 
% Reads ECG header in Mortara format. 
% 
% Arguments:
%   + filename: recording to be read or folder containing HourXXRawData.bin
%   files.
% 
% Output:
%   + heasig: header with the ECG properties. 
% 
% Limits:
% This routine is limited to read blocks smaller than 200 Mbytes for
% performance reasons. You can disable this limit by doing:
% MaxIOread = Inf; %megabytes
% 
% See also read_Mortara_format, read_ECG, ECGwrapper
% 
% Author: Mariano Llamedo Soria
% <matlab:web('mailto:llamedom@electron.frba.utn.edu.ar','-browser') (email)> 
% Version: 0.1 beta
% Birthdate: 29/05/2014
% Last update: 19/11/2014
% Copyright 2008-2014
% 
function heasig = read_Mortara_header( filename )

heasig.freq = 1000; %Hz
heasig.desc = char({'I','II','V1','V2','V3','V4','V5','V6'});
heasig.nsig = 8;

if( isdir(filename) )
    rec_path = filename;
else
    [rec_path, rec_name] = fileparts(filename);
end

hours_files = dir([rec_path filesep 'Hour*.bin' ]);

hours_available = length(hours_files);

samples_per_hour = zeros(hours_available,1);

expected_samples = 60*60*heasig.freq;

for ii = 1:hours_available
    
    hour_filename = [rec_path filesep 'Hour' num2str(ii) 'RawData.bin' ];
    fid = fopen(hour_filename);

    if( fid > 0 )
        fseek(fid, 0, 'eof');
        bytes_totales = ftell(fid);
        samples_per_hour(ii) = bytes_totales / 2 / heasig.nsig;
        fclose(fid);
        if( ii ~= expected_samples && samples_per_hour(ii) ~= expected_samples)
            error('read_Mortara:Impossible2read', 'Not enough samples in %s. Expected %d - found %d. Problem copying ?\n', hour_filename, expected_samples, samples_per_hour(ii));
        end
    else
        error('read_Mortara:Impossible2read', 'Problem opening %s.\n', hour_filename);
    end
    
end

heasig.nsamp = sum(samples_per_hour);

heasig.adczero = zeros(heasig.nsig,1);
heasig.gain = ones(heasig.nsig,1);
heasig.units = repmat('uV',heasig.nsig,1);

if( rec_path(end) == filesep )
    rec_path = rec_path(1:end-1);
end
[~, rec_name] = fileparts(rec_path);

heasig.recname = rec_name;
heasig.btime = '00:00:00';
heasig.bdate = '01/01/2000';

        
    

        
    
    