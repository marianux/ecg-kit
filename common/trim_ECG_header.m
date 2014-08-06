function [ ECG_header ] = trim_ECG_header( ECG_header, ECG_idx )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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

