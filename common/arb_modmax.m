% Description:
% 
% An arbitrary task version of the modmax function for using with ECG
% wrapper objects% 
% 
% See also ECGtask_arbitrary_function
% 
% 
% Author: Mariano Llamedo Soria (llamedom at frba.utn.edu.ar)
% Version: 0.1 beta
% Birthdate  : 5/4/2017
% Last update: 5/4/2017
% Copyright 2008-2017
% 
function payload = arb_modmax( ECG_matrix, ECG_header, ECG_start_offset, progress_handle, payload_in)

aux_modmax = modmax(ECG_matrix, payload_in.xlims, payload_in.thr, 1, round(payload_in.detection_threshold*ECG_header.freq) );

payload = colvec(aux_modmax + ECG_start_offset - 1);

