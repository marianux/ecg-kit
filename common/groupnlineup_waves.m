%% (Internal) Alineate wave positions (onset, peak and offsets commonly) 
%
%     slpos_out = groupnlineup_waves(slpos)
% 
% 
% Arguments:
% 
%   + slpos: delineation structure, 
% 
% Output:
% 
%   + slpos_out: delineation structure with all heartbeats lined up with its
%   respective waves
% 
% Example:
% 
% 
% See also alineate_positions
% 
% Author: Mariano Llamedo Soria (llamedom at frba.utn.edu.ar)
% Version: 0.1 beta
% Birthdate: 17/05/2016
% Last update: 17/05/2016
% Copyright 2008-2016
% 
function slpos_out = groupnlineup_waves(slpos)

cant_leads = length(slpos);
[slpos_mat, slpos_names] = positions2matrix(slpos);

for this_lead = 1:cant_leads

    cant_QRS = length(slpos(this_lead).qrs);
    this_pos_mat = slpos_mat{this_lead};
    cant_anns = size(this_pos_mat,1);
    
    if( isempty(this_pos_mat) )
        continue
    end
    
    prev_QRS = 1;
    
    for ii = 1:cant_QRS
        
        if( ii == cant_QRS)
            next_QRS = slpos(this_lead).qrs(cant_QRS) + 1.5* (slpos(this_lead).qrs(cant_QRS)-slpos(this_lead).qrs(cant_QRS-1));
        else            
            next_QRS = slpos(this_lead).qrs(ii+1);
        end

        slpos_out(this_lead).Pon(ii) = match_value(slpos(this_lead).Pon, slpos(this_lead).qrs(ii), prev_QRS, slpos(this_lead).qrs(ii));
        slpos_out(this_lead).P(ii) = match_value(slpos(this_lead).P, slpos(this_lead).qrs(ii), prev_QRS, slpos(this_lead).qrs(ii));
        slpos_out(this_lead).Poff(ii) = match_value(slpos(this_lead).Poff, slpos(this_lead).qrs(ii), prev_QRS, slpos(this_lead).qrs(ii));
        slpos_out(this_lead).Ptipo(ii) = nan;
        
        slpos_out(this_lead).QRSon(ii) = match_value(slpos(this_lead).QRSon, slpos(this_lead).qrs(ii), prev_QRS, slpos(this_lead).qrs(ii));
        slpos_out(this_lead).qrs(ii) = slpos(this_lead).qrs(ii);
        slpos_out(this_lead).Q(ii) = match_value(slpos(this_lead).Q, slpos(this_lead).qrs(ii), prev_QRS, next_QRS);
        slpos_out(this_lead).R(ii) = match_value(slpos(this_lead).R, slpos(this_lead).qrs(ii), prev_QRS, next_QRS);
        slpos_out(this_lead).S(ii) = match_value(slpos(this_lead).S, slpos(this_lead).qrs(ii), prev_QRS, next_QRS);
        slpos_out(this_lead).QRSoff(ii) = match_value(slpos(this_lead).QRSoff, slpos(this_lead).qrs(ii), slpos(this_lead).qrs(ii), next_QRS);

        slpos_out(this_lead).Ton(ii) = match_value(slpos(this_lead).Ton, slpos(this_lead).qrs(ii), slpos(this_lead).qrs(ii), next_QRS);
        slpos_out(this_lead).T(ii) = match_value(slpos(this_lead).T, slpos(this_lead).qrs(ii), slpos(this_lead).qrs(ii), next_QRS);
        slpos_out(this_lead).Tprima(ii) = match_value(slpos(this_lead).Tprima, slpos(this_lead).qrs(ii), slpos(this_lead).qrs(ii), next_QRS);
        slpos_out(this_lead).Toff(ii) = match_value(slpos(this_lead).Toff, slpos(this_lead).qrs(ii), slpos(this_lead).qrs(ii), next_QRS);
        slpos_out(this_lead).Ttipo(ii) = nan;
        
        prev_QRS = slpos(this_lead).qrs(ii);
        
    end
    
    if( cant_QRS == 0 )
        slposout_mat = [];
    else
        % group remaining waves if any
        slposout_mat = positions2matrix(slpos_out(this_lead));
    end
    
    cant_names = length(slpos_names);
    other_waves = cell(1,cant_names);
    
    for ii = 1:cant_names
        if( all(isnan(this_pos_mat(:,ii))) )
            other_waves{ii} = [];
        else
            if( isempty(slposout_mat) )
                aux_val = this_pos_mat(:,ii);
            else
                aux_val = setdiff(this_pos_mat(:,ii), slposout_mat(:,ii));
            end
            other_waves{ii} = aux_val(~isnan(aux_val));
        end
    end
    
    if( any(cellfun(@(a)(~isempty(a)), other_waves)) )
        % Group P waves and add as a new heartbeat
        [~, wave_idx] = intersect(slpos_names, {'Pon' 'P' 'Poff'});
        if( any(cellfun(@(a)(~isempty(a)), other_waves(wave_idx) )) )
            other_waves(wave_idx) = alineate_positions( other_waves(wave_idx) );
            aux_val = nan(length(other_waves{wave_idx(1)}),cant_names);
            aux_val(:,wave_idx) = cell2mat(other_waves(wave_idx)); 
            slposout_mat = [ slposout_mat; aux_val ];
        end
        
        % Group QRS waves and add as a new heartbeat
        [~, wave_idx] = intersect(slpos_names, {'QRSon' 'qrs' 'Q' 'R' 'S' 'QRSoff'});
        if( any(cellfun(@(a)(~isempty(a)), other_waves(wave_idx) )) )
            other_waves(wave_idx) = alineate_positions( other_waves(wave_idx) );
            aux_val = nan(length(other_waves{wave_idx(1)}),cant_names);
            aux_val(:,wave_idx) = cell2mat(other_waves(wave_idx)); 
            slposout_mat = [ slposout_mat; aux_val ];
        end
        
        % Group T waves and add as a new heartbeat
        [~, wave_idx] = intersect(slpos_names, {'Ton' 'T' 'Toff'});
        if( any(cellfun(@(a)(~isempty(a)), other_waves(wave_idx) )) )
            other_waves(wave_idx) = alineate_positions( other_waves(wave_idx) );
            aux_val = nan(length(other_waves{wave_idx(1)}),cant_names);
            aux_val(:,wave_idx) = cell2mat(other_waves(wave_idx)); 
            slposout_mat = [ slposout_mat; aux_val ];
        end
        
        slpos_out(this_lead) = matrix2positions(slposout_mat, slpos_names);
    end
        
    % sort heartbeats and waves and add as a new heartbeat
    [~, aux_idx] = sort((sum(~isnan(this_pos_mat))), 'descend');
    aux_val = slpos_out(this_lead).( slpos_names{aux_idx(1)} );

    [~, aux_idx] = sort(aux_val);
    
    for fname = slpos_names
        if( isfield(slpos_out(this_lead), fname{1}) )
            aux_val = slpos_out(this_lead).(fname{1});
            if( isempty(aux_val) )
                slpos_out(this_lead).(fname{1}) = nan(cant_anns,1);
            else
                slpos_out(this_lead).(fname{1}) = colvec(aux_val(aux_idx));
            end
        else
            slpos_out(this_lead).(fname{1}) = nan(cant_anns,1);
        end
    end
    
end


function close_val = get_closer(values, reference)

    close_val = values( min_index(abs(values-reference)) );

function vals_in_range = get_values_in_range(values, start_range, end_range)

    vals_in_range = values( values >= start_range & values <= end_range );

function matched_val = match_value(values, reference, start_range, end_range)

    aux_val = get_values_in_range(values, start_range, end_range);
    
    if( isempty(aux_val) )
        matched_val = nan;
    else
        matched_val = get_closer(aux_val, reference);
        if( isempty(matched_val) )
            matched_val = nan;
        end
    end
    
