function [ECG heasig end_sample] = read_Mortara( filename, start_sample, end_sample )
% Reads ECG recordings in Mortara format. 
% 
% Arguments:
%   + filename: recording to be read or folder containing HourXXRawData.bin
%   files.
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
% Author: Mariano Llamedo Soria (llamedom at {unizar.es;electron.frba.utn.edu.ar}
% Birthdate: 29/05/2014
% Last update: 29/05/2014
% 

heasig.freq = 1000; %Hz
heasig.desc = char({'I','II','V1','V2','V3','V4','V5','V6'});
leads_in_mortara_files = 8; 

ECG = [];
last_sample = [];

%No leer bloques mas grandes de 200 megabytes
MaxIOread = 200; %megabytes

if( nargin < 2 || isempty( start_sample ) )
    start_sample = 1;
else
    start_sample = max(1,start_sample);
end

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
        samples_per_hour(ii) = bytes_totales / 2 / leads_in_mortara_files;
        fclose(fid);
        if( ii ~= expected_samples && samples_per_hour(ii) ~= expected_samples)
            error('read_Mortara:Impossible2read', 'Not enough samples in %s. Expected %d - found %d. Problem copying ?\n', hour_filename, expected_samples, samples_per_hour(ii));
        end
    else
        error('read_Mortara:Impossible2read', 'Problem opening %s.\n', hour_filename);
    end
    
end

total_samples = sum(samples_per_hour);

if( nargin < 3 || isempty( end_sample ) )
    %Intento la lectura total por defecto
    end_sample = total_samples;
end

samples2read = end_sample - start_sample + 1;

if( (samples2read*leads_in_mortara_files*2) > (MaxIOread * 1024^2) )
    samples2read = (MaxIOread * 1024^2) / leads_in_mortara_files / 2;
    end_sample = start_sample + samples2read - 1;
    warning(['No es recomendable leer mas de ' num2str(MaxIOread) ' Mb. Realice varias lecturas.'])
end

end_samples = cumsum(samples_per_hour);
start_samples = [1; end_samples(1:end-1)+1 ];
files2read = find(start_sample <= end_samples & end_sample >= start_samples );

for ii = rowvec(files2read)

    hour_filename = [rec_path filesep 'Hour' num2str(ii) 'RawData.bin' ];
    fidECG = fopen( hour_filename, 'r');

    if( fidECG > 0 )

        if(ii == files2read(1)) 
            fseek(fid, (start_sample - start_samples(ii)) * 2 * leads_in_mortara_files, 'cof');
        end
        
        if(ii == files2read(end)) 
            if( length(files2read) == 1 )
                % read samples from a single file
                this_samples2read = samples2read;
            else
                % read samples from several files
                this_samples2read = end_sample - start_samples(ii) + 1;
            end
        else
            this_samples2read = samples_per_hour(ii);
        end
        
        ECG = [ECG; fread(fid, [leads_in_mortara_files this_samples2read], '*int16')'];

        fclose(fid);
    else
        error('read_Mortara:Impossible2read', 'Problem opening %s.\n', hour_filename);
    end
    
end

[heasig.nsamp, heasig.nsig] = size(ECG);

if( heasig.nsamp ~= samples2read )
    error('read_Mortara:Impossible2read', 'Problem reading %s.\n', filename);
end

heasig.adczero = zeros(heasig.nsig,1);
heasig.gain = ones(heasig.nsig,1);
heasig.units = repmat('uV',heasig.nsig,1);

if( rec_path(end) == filesep )
    rec_path = rec_path(1:end-1);
end
[~, rec_name] = fileparts(rec_path);

heasig.recname = rec_name;
heasig.btime = calc_btime( start_sample, heasig.freq );
heasig.bdate = '01/01/2000';

        
    

        
    
    