%% Reads ECG header in HES format
% Reads ECG headers in HES (Biosigna's) format. Implements the documentation
% provided with Biosignas database (not available with the ECHkit).
% 
% Arguments:
%   + filename: recording to be read.
% 
% Output:
%   + ann: annotations for the ECG recordings.
% 
% See also read_HES_format, read_HES_ann, read_ECG, ECGwrapper
% 
% Author: Mariano Llamedo Soria
% <matlab:web('mailto:llamedom@electron.frba.utn.edu.ar','-browser') (email)> 
% Version: 0.1 beta
% Birthdate: 17/12/2010
% Last update: 19/11/2014
% Copyright 2008-2015
% 
function heasig = read_HES_header(filename)

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
        heasig.units = repmat('nV',heasig.nsig,1);
        heasig.adcres = repmat(aux(4),1,heasig.nsig);
        heasig.adczero = repmat(2^(aux(4)-1),1,heasig.nsig);
        heasig.nsamp = fread(fidECG, 1, 'uint32');
        
        fclose(fidECG);        

        [~, heasig.recname] = fileparts(filename);
        
    catch ME
        fclose(fidECG);
        rethrow(ME)
    end
    
end
