function this_hdl = PlotWaveMarks( user_data, this_annotation, field_names, lead, yTextOffset, this_color)

    this_hdl = [];
    aux_on = colvec(this_annotation.(field_names{1})) ;
    aux_peak = colvec(this_annotation.(field_names{2})) ;
    aux_off = colvec(this_annotation.(field_names{3})) ;
    
    nsamp = size(user_data.ECG,1);
    
    aux_on_idx = find(~isnan(aux_on) & aux_on >= user_data.start_sample & aux_on <= user_data.start_sample+user_data.start_sample+nsamp);
    this_hdl = [ this_hdl; plot(user_data.axes_hdl, repmat(rowvec(aux_on(aux_on_idx)), 2, 1 ), bsxfun( @plus, repmat([-yTextOffset; yTextOffset], 1, length(aux_on_idx)), rowvec(user_data.ECG(aux_on(aux_on_idx) - user_data.start_sample + 1, lead)) ), 'Color' , this_color, 'LineStyle', ':', 'Marker', '<' , 'MarkerSize', 2, 'LineWidth', 0.25)];
    aux_off_idx = find(~isnan(aux_off) & aux_off >= user_data.start_sample & aux_off <= user_data.start_sample+nsamp);
    this_hdl = [ this_hdl; plot(user_data.axes_hdl, repmat(rowvec(aux_off(aux_off_idx)), 2, 1 ), bsxfun( @plus, repmat([-yTextOffset; yTextOffset], 1, length(aux_off_idx)), rowvec(user_data.ECG(aux_off(aux_off_idx) - user_data.start_sample + 1, lead)) ), 'Color' , this_color, 'LineStyle', ':', 'Marker', '>', 'MarkerSize', 2, 'LineWidth', 0.25 )];
    aux_peak_idx = find(~isnan(aux_peak) & aux_peak >= user_data.start_sample & aux_peak <= user_data.start_sample+nsamp);
    this_hdl = [ this_hdl; plot(user_data.axes_hdl, repmat(rowvec(aux_peak(aux_peak_idx)), 2, 1 ), bsxfun( @plus, repmat([-yTextOffset; yTextOffset]*1.3, 1, length(aux_peak_idx)), rowvec(user_data.ECG(aux_peak(aux_peak_idx) - user_data.start_sample + 1, lead)) ), 'Color' , this_color, 'LineStyle', ':', 'Marker', '^', 'MarkerSize', 2, 'LineWidth', 0.25 )];
    aux_complete_idx = find(~isnan(aux_on) & ~isnan(aux_off) & aux_on >= user_data.start_sample & aux_on <= user_data.start_sample+nsamp & aux_off >= user_data.start_sample & aux_off <= user_data.start_sample+nsamp);
    this_hdl = [ this_hdl; plot(user_data.axes_hdl,  repmat([rowvec(aux_on(aux_complete_idx));rowvec(aux_off(aux_complete_idx))], 1, 2 ) , ...
                            [ [-yTextOffset + rowvec(user_data.ECG(aux_on(aux_complete_idx) - user_data.start_sample + 1, lead)) ; -yTextOffset + rowvec(user_data.ECG(aux_off(aux_complete_idx) - user_data.start_sample + 1, lead))] [ yTextOffset + rowvec(user_data.ECG(aux_on(aux_complete_idx) - user_data.start_sample + 1, lead)); yTextOffset + rowvec(user_data.ECG(aux_off(aux_complete_idx) - user_data.start_sample + 1, lead)) ] ], 'Color' , this_color, 'LineStyle', ':', 'Marker', 'none', 'LineWidth', 0.25 )];
