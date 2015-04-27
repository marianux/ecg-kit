%% (Internal) Example of user-created QRS detector
%
% Description:
% An ECG delineator that predict wave locations just by chance, i.e.
% without using the ECG. 
% 
% Interface to follow in order to use your own detector with Wrapper and
% task objects.% 
% 
% See also ECGtask_ECG_delineation
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate  : 30/7/2014
% Last update: 30/7/2014
% Copyright 2008-2015
% 
function [positions_single_lead, position_multilead] = example_worst_ever_ECG_delineator( ECG_matrix, ECG_header, progress_handle, payload_in)


cAnnotationFields = { 'Pon' 'P' 'Poff' 'QRSon' 'qrs' 'Q' 'R' 'S' 'QRSoff' 'Ton' 'T' 'Toff' };

positions_single_lead = [];

for ii = 1:ECG_header.nsig
    
    progress_handle.checkpoint(['Guessing waves on lead ' ECG_header.desc(ii,:)])

    aux_struct = [];
    for fn = cAnnotationFields
        aux_struct.(fn{1}) = colvec(sort(randsample(1:ECG_header.nsamp, round(ECG_header.nsamp/ECG_header.freq) )));
    end
    
    if( ii == 1)
        positions_single_lead = aux_struct;
    else
        positions_single_lead(ii) = aux_struct;
    end
    
end

progress_handle.checkpoint('Calculating multilead delineation')

if( ~isempty(positions_single_lead) )
    position_multilead = positions_single_lead(1);
else
    position_multilead = [];
end

