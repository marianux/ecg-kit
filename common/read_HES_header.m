function heasig = read_HES_header(filename)

% Reads ECG headers in HES (Biosigna's) format. Implements the documentation
% provided with Biosignas database.
% 
% Arguments:
%   + filename: recording to extract header block.
% 
% Output:
%   + heasig: header with the ECG properties. 
% 
% Limits:
% Unknown. Any feedback is very welcome.
% 
% Author: Mariano Llamedo Soria (llamedom at {unizar.es;electron.frba.utn.edu.ar}
% Birthdate: 7/1/2011
% Last update: 7/1/2011
% 


heasig = [];

tablas_y_constantes;

%default settings for the whole DB
fidECG = fopen( filename, 'r');

if( fidECG > 0 )

    try

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
        heasig.gain = repmat(aux(3),1,heasig.nsig);
        heasig.units = repmat('nV/adc_unit',heasig.nsig,1);
        heasig.adcres = repmat(aux(4),1,heasig.nsig);
        heasig.adczero = repmat(2^(aux(4)-1),1,heasig.nsig);
        heasig.nsamp = fread(fidECG, 1, 'uint32');
        
        fclose(fidECG);        
        
    catch ME
        fclose(fidECG);
        rethrow(ME)
    end
    
end
