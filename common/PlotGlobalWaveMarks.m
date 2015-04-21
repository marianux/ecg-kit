%% (Internal) Internal function of plot_ecg_strip
%   
%   this_hdl = PlotGlobalWaveMarks( user_data, field_names, limits, this_color)
% 
% 
% Example:
% 
% 
% See also plot_ecg_strip
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function this_hdl = PlotGlobalWaveMarks( user_data, field_names, limits, this_color)

    this_hdl = [];
    aux_on = user_data.global_annotations.(field_names{1});
    aux_peak = user_data.global_annotations.(field_names{2});
    aux_off = user_data.global_annotations.(field_names{3});
    
    aux_on_idx = find(~isnan(aux_on));
    this_hdl = [ this_hdl; plot(user_data.axes_hdl, repmat(rowvec(aux_on(aux_on_idx)), 2, 1 ), repmat(colvec(limits), 1, length(aux_on_idx)), 'Color' , this_color, 'LineStyle', ':', 'Marker', '<' , 'MarkerSize', 4, 'LineWidth', 0.25)];
    aux_off_idx = find(~isnan(aux_off));
    this_hdl = [ this_hdl; plot(user_data.axes_hdl, repmat(rowvec(aux_off(aux_off_idx)), 2, 1 ), repmat(colvec(limits), 1, length(aux_off_idx)), 'Color' , this_color, 'LineStyle', ':', 'Marker', '>', 'MarkerSize', 4, 'LineWidth', 0.25 )];
    aux_peak_idx = find(~isnan(aux_peak));
    this_hdl = [ this_hdl; plot(user_data.axes_hdl, repmat(rowvec(aux_peak(aux_peak_idx)), 2, 1 ), repmat(colvec(limits) + [-0.02; 0.02] * diff(limits), 1, length(aux_peak_idx)), 'Color' , this_color, 'LineStyle', ':', 'Marker', '^', 'MarkerSize', 4, 'LineWidth', 0.25 )];
    aux_complete_idx = find(~isnan(aux_on) & ~isnan(aux_off));
    this_hdl = [ this_hdl; plot(user_data.axes_hdl,  repmat([rowvec(aux_on(aux_complete_idx));rowvec(aux_off(aux_complete_idx))], 1, 2 ) , ...
                            [ repmat(limits(1), 2, length(aux_complete_idx)) repmat(limits(2),2,length(aux_complete_idx)) ], 'Color' , this_color, 'LineStyle', ':', 'Marker', 'none', 'MarkerSize', 4, 'LineWidth', 0.25 )];
