function [ECG heasig ann last_sample] = read_HES_format( filename, start_sample, end_sample )
% Reads ECG recordings in HES (Biosigna) format. Implements the documentation
% available in the help of the application provided with the database.
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
% Author: Mariano Llamedo Soria (llamedom at {unizar.es;electron.frba.utn.edu.ar}
% Birthdate: 17/12/2010
% Last update: 18/2/2013
% 
ann = [];
heasig = [];
ECG = [];
last_sample = [];

tablas_y_constantes;

%No leer bloques mas grandes de 200 megabytes
MaxIOread = 200; %megabytes

if( nargin < 2 || isempty( start_sample ) )
    start_sample = 1;
else
    start_sample = max(1,start_sample);
end

fidECG = fopen( filename, 'r');

if( fidECG > 0 )

    fseek(fidECG, 900, 'bof');
    
    num_of_leads = fread(fidECG, 1, 'uint8');
    num_of_leads_simult = fread(fidECG, 1, 'uint8');
    
    if( num_of_leads ~= num_of_leads_simult)
        error('Sampleo no simultáneo, revisar.')
    end
    
    heasig.nsig = num_of_leads_simult;
    
    lead_description_idx = fread(fidECG, heasig.nsig, 'uint8');
    [dummy lead_description_table_idx] = intersect(Lead_description_idx, lead_description_idx);
    heasig.desc = char(cLead_description_table(lead_description_table_idx,1));

    fseek(fidECG, 992, 'bof');
    
    aux = fread(fidECG, 5, 'uint16');
    
    sample_interval_usec = aux(1);
    heasig.freq = 1/sample_interval_usec*1e6;
    heasig.gain = repmat(1/aux(3),1,heasig.nsig);
    heasig.units = repmat('nV',heasig.nsig,1);
    heasig.adcres = repmat(aux(4),1,heasig.nsig);
    heasig.adczero = repmat(2^(aux(4)-1),1,heasig.nsig);
    heasig.nsamp = fread(fidECG, 1, 'uint32');
    
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
    
    ECG = nan(samples2read,heasig.nsig);
    
    try 
        
        fseek(fidECG, 1024+((start_sample-1)*heasig.nsig)*2, 'bof');

        ECG = fread(fidECG, [heasig.nsig samples2read], '*int16')';

        fclose(fidECG);

    catch ME
        fclose(fidECG);
        rethrow(ME)
    end

    last_sample = size(ECG,1) + start_sample - 1;
    
    if( nargout > 2 )
        ann = read_HES_ann([ filename(1:end-4) '.lst' ]);
        ann.time = round(ann.time * heasig.freq);
    end

end
