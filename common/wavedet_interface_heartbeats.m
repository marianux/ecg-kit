function template_positions = wavedet_interface_heartbeats( ecg_template, header, qrs_idx )
        
    template_size = size(ecg_template,1);
    my_win = ones(template_size,1);
    N_half_win = 64;
    aux_win = blackman(2*N_half_win);
    my_win(1:N_half_win) = my_win(1:N_half_win) .* aux_win(1:N_half_win);
    my_win(end-N_half_win+1:end) = flipud(my_win(1:N_half_win));

    mean_fetal_heartbeat = bsxfun( @times, ecg_template, my_win);

    repeat_beats = 7;
    header.nsamp = repeat_beats*template_size;
    
    % line for debug
    %     plot_ecg_mosaic(mean_fetal_heartbeat, [], [], 2, 2)

    template_location = qrs_idx:template_size:repeat_beats*template_size;
    [~, all_positions] = wavedet_interface(repmat(mean_fetal_heartbeat,repeat_beats,1), header, template_location );

    % line for debug
    %     plot_ecg_strip( repmat(mean_fetal_heartbeat,repeat_beats,1) , 'ECG_header', header, 'Annotations', all_positions)

    for ii = 1:header.nsig
        for fname = rowvec(fieldnames(all_positions(ii)))
            aux_val = all_positions(ii).(fname{1});
            aux_val = aux_val(end) - (template_location(end) - qrs_idx + 1) + 1;
            if( aux_val < 1 || aux_val > template_size )
                aux_val = nan;
            end
            template_positions(ii).(fname{1}) = aux_val;
        end
    end
    
    