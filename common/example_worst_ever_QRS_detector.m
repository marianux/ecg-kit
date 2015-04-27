function [positions_single_lead, position_multilead] = example_worst_ever_QRS_detector( ECG_matrix, ECG_header, progress_handle, payload_in)

% Example of user-created QRS detector
% ------------------------------------
% Description:
% A QRS fortune-teller that predicts QRS locations just by chance, i.e.
% without using the ECG. 
% 
% Interface to follow in order to use your own detector with Wrapper and
% task objects.% 
% 
% 
% See also ECGtask_QRS_detection
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate  : 30/7/2014
% Last update: 30/7/2014

positions_single_lead = cell(1, ECG_header.nsig);

for ii = 1:ECG_header.nsig

    progress_handle.checkpoint(['Guessing detections on lead ' ECG_header.desc(ii,:)])
    
    positions_single_lead{ii} = colvec(sort(randsample(1:ECG_header.nsamp, round(ECG_header.nsamp/ECG_header.freq) )));
    
end

progress_handle.checkpoint('Calculating multilead detections')

position_multilead = colvec(round(mean(cell2mat(positions_single_lead),2)));
