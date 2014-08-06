function [PPG_ABP_idx, ECG_header ] = get_PPG_ABP_idx_from_header(ECG_header)
    
    PPG_ABP_idx = [];
    
    str_12_leads_desc = { 'PPG' 'PLETH' 'ABP' 'ART' 'BP' 'PULSEWAVE' 'SPO2' };
    % add ECG patterns to be matched in ECG_header.desc
    PPG_ABP_pattern_desc = {'PLETH' 'ABP' 'ART' 'BP' 'PPG' };

    [~, ~, aux_idx ] = intersect(str_12_leads_desc, upper(strtrim(cellstr(ECG_header.desc))) );
    
    if( ECG_header.nsig ~= length(aux_idx) )
        
        PPG_ABP_idx = colvec(aux_idx);
        
        aux_idx = find(any(cell2mat(cellfun(@(b)(cell2mat(cellfun(@(a)(~isempty(strfind(a, b))), rowvec(upper(cellstr(ECG_header.desc))), 'UniformOutput', false))), colvec(PPG_ABP_pattern_desc), 'UniformOutput', false)),1));
        
        if( isempty(aux_idx) && isempty(PPG_ABP_idx) )
            if( isempty(PPG_ABP_idx) )
                cprintf('[1,0.5,0]', disp_option_enumeration('Could not find any PPG/ABP signal, check the lead description of the recording:', cellstr(ECG_header.desc) ) );
                fprintf(1, '\n')
            end
        else
            PPG_ABP_idx = unique([PPG_ABP_idx; colvec(aux_idx)]);
            Not_ECG_idx = setdiff(1:ECG_header.nsig, PPG_ABP_idx);
            if( ~isempty(Not_ECG_idx) )
%                 warning('get_ECG_idx_from_header:SignalsDiscarded', disp_option_enumeration('Some signal/s present are not ECG:', cellstr(ECG_header.desc(Not_ECG_idx,:)) ))
            end
        end
        
    else
        % all leads are standard ECG leads
        PPG_ABP_idx = 1:ECG_header.nsig;
    end
    
    if(length(PPG_ABP_idx) ~= ECG_header.nsig && nargout > 1)
       %trim the heasig if needed 
       
        ECG_header = trim_ECG_header(ECG_header, PPG_ABP_idx);
        
    end
    