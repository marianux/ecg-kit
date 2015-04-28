%% (Internal) Trim a header info struct to a subset of signals
%
%     [ ECG_header ] = trim_ECG_header( ECG_header, ECG_idx )
% 
% Arguments:
% 
%   + ECG_header: original header struct.
% 
%   + ECG_idx: signal indexes to select.
% 
% Output:
% 
%   + ECG_header: trimed header struct.
% 
% See also ECGtask_arbitrary_function, ECGtask_ECG_delineation, ECGtask_QRS_detection
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate: 21/7/2010
% Last update: 20/02/2013
% Copyright 2008-2015
% 
function [ ECG_header ] = trim_ECG_header( ECG_header, ECG_idx )

       for fname = rowvec( fieldnames(ECG_header) )
           aux_val = ECG_header.(fname{1});
           [nsig, n2] = size(aux_val);
           if( nsig == ECG_header.nsig )
               % filter only the ECG_idx signals
               if( n2 > 1 )
                   ECG_header.(fname{1}) = aux_val(ECG_idx,:);
               else
                   ECG_header.(fname{1}) = aux_val(ECG_idx);
               end
           end
       end
       
       ECG_header.nsig = length(ECG_idx);
        

end

