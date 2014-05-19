function ECG_hdl = plot_ecg_heartbeat(ECG, QRS_locations, QRS_start_idx, cant_qrs, heasig, filtro, axes_hdl )

% obsolete, use plot_ecg_strip.m

if( nargin < 6 )
    filtro = [];
end

if( nargin < 7 || isempty(axes_hdl) )
    axes_hdl = gca;
end

ECG_hdl = [];

% cla(axes_hdl);
% axes(axes_hdl);
fig_hdl = gcf;

set(fig_hdl, 'CurrentAxes', axes_hdl);

[lECG cant_sig] = size(ECG);

if( iscell(QRS_locations) )
    other_QRS_locations = QRS_locations(2:end);
    QRS_locations = QRS_locations{1};
else
    other_QRS_locations = [];
end

cant_QRS_locations = length(QRS_locations);

QRS_start_idx = min(cant_QRS_locations, QRS_start_idx);

if( isempty(QRS_locations) ) 
    start_idx = 1;
    end_idx = 10*heasig.freq;
else
    
    if( (QRS_start_idx-cant_qrs) < 1 )
        start_idx = max(1, QRS_locations(QRS_start_idx) - 5 * heasig.freq );
    else
        
        % avoid gaps between marks
        aux_start = max(1, QRS_start_idx-cant_qrs);
        aux_marks = QRS_locations(aux_start:QRS_start_idx);
        aux_diffs = diff(aux_marks);
        median_diff = median(aux_diffs);
        gap_idx = find( aux_diffs > 5*median_diff , 1, 'last' );
        
        if( isempty(gap_idx) )
            start_idx = max(1, QRS_locations(aux_start));
        else
            start_idx = max(1, aux_marks(gap_idx+1) - 5 * heasig.freq);
        end
            
    end
    
    if( QRS_start_idx >= (cant_QRS_locations-cant_qrs) )
        end_idx = min(lECG, QRS_locations(QRS_start_idx) + 5 * heasig.freq );
    else
        
        aux_end = min(cant_QRS_locations, QRS_start_idx+cant_qrs);
        aux_marks = QRS_locations(QRS_start_idx:aux_end);
        aux_diffs = diff(aux_marks);
        median_diff = median(aux_diffs);
        gap_idx = find( 5*median_diff < aux_diffs, 1, 'last' );
        
        if( isempty(gap_idx) )
            end_idx = min(lECG, QRS_locations(aux_end));
        else
            end_idx = min(lECG, aux_marks(gap_idx) + 5 * heasig.freq);
        end
        
    end
    
    
end
aux_idx = start_idx:end_idx;


if( ~isempty(filtro) && ~isempty(aux_idx) )
    % zero phase filtering
    aux_idx2 = max(1, aux_idx(1) - 1 * heasig.freq  ):min(lECG, aux_idx(end) + 1 * heasig.freq );
    ECG(aux_idx2,:) = filter(filtro, flipud(ECG(aux_idx2,:)) );
    ECG(aux_idx2,:) = filter(filtro, flipud(ECG(aux_idx2,:)) );
end

if( cant_sig > 1 )
    max_values = max(ECG(aux_idx,:));
    min_values = min(ECG(aux_idx,:));
    aux_ranges = max_values - min_values;
    aux_sig = bsxfun( @minus, ECG(aux_idx,:), mean(ECG(aux_idx,:)));
    aux_sig = bsxfun( @times, aux_sig, max(aux_ranges)./aux_ranges );
    ecg_max = max(max(aux_sig));
    
else
    aux_sig = ECG(aux_idx,:);
    ecg_max = max(aux_sig);
end

if( isempty(aux_idx) )
    qrs_ploted =[];
else
    qrs_ploted = find(QRS_locations >= aux_idx(1) & QRS_locations <= aux_idx(end) );
end

ECG_hdl = plot(axes_hdl, aux_idx, aux_sig );
hold(axes_hdl, 'on')

ECG_hdl = [ECG_hdl; plot(axes_hdl, repmat(rowvec(QRS_locations(qrs_ploted)),2,1), [ repmat(ecg_max, 1,length(qrs_ploted)) ; zeros(1,length(qrs_ploted))] , 'r--' )];

if( ~isempty(qrs_ploted) )
    for ii = rowvec(qrs_ploted)
        aux_hdl = text( QRS_locations(ii), double(1.1*ecg_max), num2str(ii) );
    end

    aux_text_pos = get(aux_hdl, 'Extent');

    ytext_loc = aux_text_pos(2)+aux_text_pos(4);

    aux_val = unique(QRS_locations(qrs_ploted));
    set(axes_hdl, 'XTick', rowvec(aux_val) );
    set(axes_hdl, 'XTickLabel', Seconds2HMS(colvec(aux_val)*1/heasig.freq) );

end

aux_yrange = get(axes_hdl, 'Ylim');

if( isempty(other_QRS_locations) )
    
    ytext_loc = aux_yrange(2);
    
else
   
    ColorOrder = hsv(12);
    ii = 1;
    
    for this_QRS_location = other_QRS_locations
        this_QRS_locations = this_QRS_location{1};
        qrs_ploted = find(this_QRS_locations >= aux_idx(1) & this_QRS_locations <= aux_idx(end) );
        ECG_hdl = [ECG_hdl; plot(axes_hdl, repmat(rowvec(this_QRS_locations(qrs_ploted)),2,1), [ repmat(ecg_max, 1,length(qrs_ploted)) ; zeros(1,length(qrs_ploted))] , 'LineStyle', ':', 'Color', ColorOrder(ii,:) )];
        for ii = rowvec(qrs_ploted)
            text( this_QRS_locations(ii), ytext_loc, num2str(ii), 'FontSize', 7 );
        end
        
        ytext_loc = ytext_loc+aux_text_pos(4)/3;
        ii = ii + 1;
        
    end
    
end

aux_yrange = [aux_yrange(1) ytext_loc ];

set(axes_hdl, 'Ylim', aux_yrange);
set(axes_hdl, 'Box', 'off');

hold(axes_hdl, 'off');

