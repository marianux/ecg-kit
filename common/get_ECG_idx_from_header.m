%% (Internal) Guess ECG signals indexes in a multimodal recording
%   
%   [ECG_idx, ECG_header ]= get_ECG_idx_from_header(ECG_header)
% 
% Arguments:
% 
%      + ECG_header: header info
% 
% Output:
% 
%      + ECG_idx: indexes of the ECG signals in the multimodal matrix
% 
%      + ECG_header: header of a ECG only recording
% 
% Example:
% 
% 
% See also get_PPG_ABP_idx_from_header, ECGwrapper
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function [ECG_idx, ECG_header ]= get_ECG_idx_from_header(ECG_header)
    
    ECG_idx = [];

    str_12_leads_desc = {'I' 'II' 'III' 'AVR' 'AVL' 'AVF' 'V1' 'V2' 'V3' 'V4' 'V5' 'V6' ... %Standard leads.
                         'MLII' 'MLI' 'MLIII' 'MV1' 'MV2' 'MV3' 'MV4' 'MV5' 'MV6' ... % Modifyied leads.
                         'X' 'Y' 'Z' 'VX' 'VY' 'VZ' ... % Orthogonal leads
                         'E-S' 'A-S' 'A-I' ... % Easy leads
                         };
                     
    % add ECG patterns to be matched in ECG_header.desc
    ECG_pattern_desc = {'ECG' 'BIPOLAR' 'LEAD'};

    [~, ~, aux_idx ] = intersect(str_12_leads_desc, upper(strtrim(cellstr(ECG_header.desc))) );
    
    if( ECG_header.nsig ~= length(aux_idx) )
        
        ECG_idx = colvec(aux_idx);
        
        aux_idx = find(any(cell2mat(cellfun(@(b)(cell2mat(cellfun(@(a)(~isempty(strfind(a, b))), rowvec(upper(cellstr(ECG_header.desc))), 'UniformOutput', false))), colvec(ECG_pattern_desc), 'UniformOutput', false)),1));
        
        if( isempty(aux_idx) && isempty(ECG_idx) )
            if( isempty(ECG_idx) )
                cprintf('[1,0.5,0]', disp_option_enumeration('Could not find any ECG signal, check the lead description of the recording:', cellstr(ECG_header.desc) ) );
                fprintf(1, '\n')
            end
        else
            ECG_idx = [ECG_idx; colvec(aux_idx)];
            Not_ECG_idx = setdiff(1:ECG_header.nsig, ECG_idx);
            if( ~isempty(Not_ECG_idx) )
%                 warning('get_ECG_idx_from_header:SignalsDiscarded', disp_option_enumeration('Some signal/s present are not ECG:', cellstr(ECG_header.desc(Not_ECG_idx,:)) ))
            end
        end
        
        ECG_idx = sort(ECG_idx);
        
    else
        % all leads are standard ECG leads
        ECG_idx = 1:ECG_header.nsig;
    end
    
    if(length(ECG_idx) ~= ECG_header.nsig && nargout > 1)
       %trim the heasig if needed 
        ECG_header = trim_ECG_header(ECG_header, ECG_idx);
    end
    