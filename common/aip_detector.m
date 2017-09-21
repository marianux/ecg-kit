function [payload, interproc_data ] = aip_detector( ECG_matrix, ECG_header, ECG_start_offset, progress_handle, payload_in, interproc_data)
% aip_detector Arbitrary impulsive pseudoperiodic detector
%   Detailed explanation goes here

    payload = [];
    
    payload.series_quality.AnnNames = {};
    payload.series_quality.ratios = [];
    payload.series_quality.estimated_labs = {};

    %% config parsing
    
    if( ~isfield(payload_in, 'trgt_width') )
        % Default QRS width
        payload_in.trgt_width = 0.06; % seconds
    end
    
    if( ~isfield(payload_in, 'trgt_min_pattern_separation') )
        % Default QRS width
        payload_in.trgt_min_pattern_separation = 0.3; % seconds
    end
    
    if( ~isfield(payload_in, 'trgt_max_pattern_separation') )
        % Default QRS width
        payload_in.trgt_max_pattern_separation = 2; % seconds
    end

    if( ~isfield(payload_in, 'stable_RR_time_win') )
        % Default time window to consider an stable rhythm
        payload_in.stable_RR_time_win = 2; % seconds
    end
    
    if( ~isfield(payload_in, 'final_build_time_win') )
        % Default time window to build the final detection based on each
        % pattern found in the recording
        payload_in.final_build_time_win = 20; % seconds
    end
    
    if( ~isfield(payload_in, 'powerline_interference') )
        % Default QRS width
        payload_in.powerline_interference = 50; % Hz
    end
    
    if( ~isfield(payload_in, 'max_patterns_found') )
        % Default QRS width
        payload_in.max_patterns_found = 3; % # patterns
    end
    
    %% start
    
%     lp_size = 11;
    pattern_size = 2*round(payload_in.trgt_width/2*ECG_header.freq)+1; % force odd number
    
    pattern_coeffs = diff(gausswin(pattern_size+1)) .* gausswin(pattern_size);
    
    progress_handle.Loops2Do = ECG_header.nsig;

    pwl_size = round(ECG_header.freq / payload_in.powerline_interference);

    lp_size = round(1.2*pattern_size);
    
    progress_handle.start_loop();    

    %% Initial Pattern search 
    
    progress_handle.checkpoint('Initial filtering');

%         rise_detector = filter(ones(lp_size,1)/lp_size,1, abs([0; diff(double(ECG_matrix(:,ss)))]));
%         rise_detector = [zeros(((lp_size-1)/2)-1,1); rise_detector((lp_size-1)/2:end)];
    rise_detector = filter( ones(pwl_size,1)/pwl_size, 1, flipud(double(ECG_matrix(:,1))) );
    rise_detector = filter( ones(pwl_size,1)/pwl_size, 1, flipud(rise_detector) );
    rise_detector = filter( pattern_coeffs, 1, flipud(rise_detector) );
    rise_detector = filter( pattern_coeffs, 1, flipud(rise_detector) );
    rise_detector = filter( ones(lp_size,1)/lp_size, 1, flipud(abs(rise_detector)) );
    rise_detector = filter( ones(lp_size,1)/lp_size, 1, flipud(rise_detector) );

% figure(3); plot(ECG_matrix(:,1) ); xlims = xlim(); ylims = ylim(); box_h = ylims + [0.01 -0.1] * diff(ylims); aux_sc = 0.8*diff(ylims)/(max(rise_detector)-min(rise_detector));aux_off = ylims(1) + 0.6*diff(ylims); hold on; plot(rise_detector * aux_sc + aux_off); hold off; ylim(ylims);

%     % maybe an statistical check for assimetry and kurtosis ???
%     [skewness(rise_detector) kurtosis(rise_detector)]

    progress_handle.checkpoint('Peak detection');

    initial_thr = 30; % percentil
    first_bin_idx = 2; % bin donde comenzar a buscar el min en el histograma

    actual_thr = prctile(rise_detector, initial_thr);

    [~, max_values ] = modmax(rise_detector, 1, actual_thr, 0, round(payload_in.trgt_min_pattern_separation * ECG_header.freq));
    
    prctile_grid = prctile( max_values, 1:100 );
    
    grid_step = median(diff(prctile_grid));
    
    thr_grid = actual_thr: grid_step:max(max_values);
    
    hist_max_values = histcounts(max_values, thr_grid);

    [thr_idx, thr_max ] = modmax( colvec( hist_max_values ) , first_bin_idx, 0, 0, [], 10);
    
    thr_idx_expected = floor(rowvec(thr_idx) * colvec(thr_max) *1/sum(thr_max));

    aux_seq = 1:length(thr_grid);
    
    min_hist_max_values = min( hist_max_values( aux_seq >= first_bin_idx & aux_seq < thr_idx_expected) );
    
    thr_min_idx = round(mean(find(aux_seq >= first_bin_idx & aux_seq < thr_idx_expected & [hist_max_values 0] == min_hist_max_values)));
    
    actual_thr = thr_grid(thr_min_idx );

    first_detection_idx = modmax(rise_detector, 1, actual_thr, 0, round(payload_in.trgt_min_pattern_separation * ECG_header.freq));
    
    %% look for stable segments
    
    RRserie = RR_calculation(first_detection_idx, ECG_header.freq);               
    
    RRserie_filt = MedianFiltSequence(first_detection_idx, RRserie, round(5 * ECG_header.freq));
    
    RR_scatter = abs(RRserie - RRserie_filt)./RRserie_filt;
    
    % find the most stable regions
    RR_thr = prctile(RR_scatter, 50);
    
    stable_rhythm_regions = get_segments_from_sequence(first_detection_idx, RR_scatter, round(payload_in.stable_RR_time_win * ECG_header.freq), RR_thr );

% figure(3); plot(first_detection_idx, [RRserie RRserie_filt] ); ylims = ylim(); ylims = ylims + [0.1 -0.1] * diff(ylims); hold on; plot( [ stable_rhythm_regions'; flipud(stable_rhythm_regions'); stable_rhythm_regions(:,1)' ],  repmat([ylims(1); ylims(1); ylims(2);ylims(2);ylims(1)],1, size(stable_rhythm_regions,1) ), 'b:' ); hold off;
% figure(3); plot(first_detection_idx, [RR_scatter], ':xb' ); xlims = xlim(); ylims = ylim(); ylims = ylims + [0.01 -0.1] * diff(ylims); hold on; plot(xlims, [RR_thr RR_thr], '--r'); plot( [ stable_rhythm_regions'; flipud(stable_rhythm_regions'); stable_rhythm_regions(:,1)' ],  repmat([ylims(1); ylims(1); ylims(2);ylims(2);ylims(1)],1, size(stable_rhythm_regions,1) ), 'b:' ); hold off;
% figure(3); plot(first_detection_idx, [RRserie RRserie_filt] ); xlims = xlim(); ylims = ylim(); ylims = ylims + [0.1 -0.1] * diff(ylims); aux_sc = 0.25*diff(ylims)/(max(RR_scatter)-min(RR_scatter));aux_off = ylims(1) + 0.65*diff(ylims); hold on; plot(first_detection_idx, (RR_scatter*aux_sc) +aux_off , ':xb' ); plot(xlims, ([RR_thr RR_thr]*aux_sc ) + aux_off, '--r'); plot( [ stable_rhythm_regions'; flipud(stable_rhythm_regions'); stable_rhythm_regions(:,1)' ],  repmat([ylims(1); ylims(1); ylims(2);ylims(2);ylims(1)],1, size(stable_rhythm_regions,1) ), 'b:' ); hold off;    

    if( isempty(stable_rhythm_regions) )
        cprintf('Red', '\nNo stable rhythm segments found. Consider decreasing "stable_RR_time_win" (Current %f seconds)\n\n', payload_in.stable_RR_time_win);                    
        return
    end

    lstable_rhythm_regions = size(stable_rhythm_regions,1);
    
    % get the larger or more stable segments first
    [~, aux_idx] = sort(diff(stable_rhythm_regions,1,2), 'descend');

    larger_stable_rank = arrayfun(@(a)(find(a==aux_idx)), 1:lstable_rhythm_regions);
    
    pattern_prewin = round(payload_in.trgt_width/2*ECG_header.freq); % force odd number
    
    % get the segments with less change in morphology
    pack_variance = repmat(realmax,lstable_rhythm_regions,1);
    
    % from the longest - most rhythm stable regions
    for ii = 1: min( 2*payload_in.max_patterns_found, round(lstable_rhythm_regions/4) )
        
        jj = aux_idx(ii);
        % refinamos con el metodo de woody
        aux_idx2 = find(first_detection_idx >= stable_rhythm_regions(jj,1) & first_detection_idx <= stable_rhythm_regions(jj,2));
        avg_pack = pack_signal(ECG_matrix(:,1), first_detection_idx(aux_idx2), [ pattern_prewin pattern_prewin ], true);    

        pack_variance(jj) = mean(var(squeeze(double(avg_pack))'));
        
    end
    [~, aux_idx] = sort(pack_variance);
    
    morphology_rank = arrayfun(@(a)(find(a==aux_idx)), 1:lstable_rhythm_regions);
    
    % final rank
    
    [~, aux_idx] = sort(larger_stable_rank + morphology_rank);
    
    stable_rhythm_regions = stable_rhythm_regions(aux_idx( 1:min(payload_in.max_patterns_found , size(stable_rhythm_regions,1))), : );

% figure(3); plot(first_detection_idx, [RRserie RRserie_filt] ); xlims = xlim(); ylims = ylim(); ylims = ylims + [0.1 -0.1] * diff(ylims); aux_sc = 0.25*diff(ylims)/(max(RR_scatter)-min(RR_scatter));aux_off = ylims(1) + 0.65*diff(ylims); hold on; plot(first_detection_idx, (RR_scatter*aux_sc) +aux_off , ':xb' ); plot(xlims, ([RR_thr RR_thr]*aux_sc ) + aux_off, '--r'); plot( [ stable_rhythm_regions'; flipud(stable_rhythm_regions'); stable_rhythm_regions(:,1)' ],  repmat([ylims(1); ylims(1); ylims(2);ylims(2);ylims(1)],1, size(stable_rhythm_regions,1) ), 'b:' ); arrayfun(@(a)(text( stable_rhythm_regions(a,1), ylims(2), num2str(a))), 1:min(payload_in.max_patterns_found )); hold off;
    
%     RRserie_filt = {[colvec(first_detection_idx) colvec(RRserie_filt)]};
    
    %% save first detections as a separate time serie
    
    str_aux = [ 'aip_first_guess' ];
    payload.(str_aux).time = first_detection_idx + ECG_start_offset - 1;
    payload.series_quality.AnnNames = [ payload.series_quality.AnnNames ; {str_aux} {'time'} ];
    RRserie_mean_sq_error = mean((RRserie - RRserie_filt).^2);
    payload.series_quality.ratios = [payload.series_quality.ratios; RRserie_mean_sq_error];
    payload.series_quality.estimated_labs = [ payload.series_quality.estimated_labs; {[]} ];
    

    for ii = 1:size(stable_rhythm_regions,1)
        
        % refinamos con el metodo de woody
        aux_idx2 = find(first_detection_idx >= stable_rhythm_regions(ii,1) & first_detection_idx <= stable_rhythm_regions(ii,2));
        
% figure(3); avg_pack = pack_signal(ECG_matrix(:,1), first_detection_idx(aux_idx2), 10*[ pattern_prewin pattern_prewin ], true); plot_ecg_mosaic(avg_pack)
        
        pattern_coeffs = woody_method(ECG_matrix(:,1), first_detection_idx(aux_idx2), [ pattern_prewin (pattern_prewin+1) ], 0.95);

        progress_handle.checkpoint(['Pattern match ' num2str(ii) ]);

    %         rise_detector = filter(ones(lp_size,1)/lp_size,1, abs([0; diff(double(ECG_matrix(:,ss)))]));
    %         rise_detector = [zeros(((lp_size-1)/2)-1,1); rise_detector((lp_size-1)/2:end)];
        rise_detector = filter( pattern_coeffs, 1, flipud(double(ECG_matrix(:,1))) );
        rise_detector = filter( pattern_coeffs, 1, flipud(rise_detector) );
        rise_detector = filter( ones(pwl_size,1)/pwl_size, 1, flipud(rise_detector) );
        rise_detector = filter( ones(pwl_size,1)/pwl_size, 1, flipud(rise_detector) );
        rise_detector = filter( ones(lp_size,1)/lp_size, 1, flipud(abs(rise_detector)) );
        rise_detector = filter( ones(lp_size,1)/lp_size, 1, flipud(rise_detector) );

        progress_handle.checkpoint('Peak detection');

        actual_thr = prctile(rise_detector, initial_thr);

        [~, max_values ] = modmax(rise_detector, 1, actual_thr, 0, round(payload_in.trgt_min_pattern_separation * ECG_header.freq));

        prctile_grid = prctile( max_values, 1:100 );

        grid_step = median(diff(prctile_grid));

        thr_grid = actual_thr: grid_step:max(max_values);

        hist_max_values = histcounts(max_values, thr_grid);

        [thr_idx, thr_max ] = modmax( colvec( hist_max_values ) , first_bin_idx, 0, 0, [], 10);

        thr_idx_expected = floor(rowvec(thr_idx) * colvec(thr_max) *1/sum(thr_max));

        aux_seq = 1:length(thr_grid);

        min_hist_max_values = min( hist_max_values( aux_seq >= first_bin_idx & aux_seq < thr_idx_expected) );

        thr_min_idx = round(mean(find(aux_seq >= first_bin_idx & aux_seq < thr_idx_expected & [hist_max_values 0] == min_hist_max_values)));

        actual_thr = thr_grid(thr_min_idx );

        this_idx = modmax(rise_detector, 1, actual_thr, 0, round(payload_in.trgt_min_pattern_separation * ECG_header.freq));
    
                   
        RRserie = RR_calculation(this_idx, ECG_header.freq);               
        
        this_RRserie_filt = colvec(MedianFiltSequence(this_idx, RRserie, round(5 * ECG_header.freq)));

        RRserie_mean_sq_error = mean((RRserie - this_RRserie_filt).^2);
        
%         RRserie_filt = [RRserie_filt; {[colvec(this_idx) RRserie]}];        
        
        % cargo la estructura de detecciones
        str_aux = [ 'aip_auto_' num2str(ii) ];
        payload.(str_aux).time = this_idx + ECG_start_offset - 1;
        payload.series_quality.AnnNames = [ payload.series_quality.AnnNames ; {str_aux} {'time'} ];
        payload.series_quality.ratios = [payload.series_quality.ratios; RRserie_mean_sq_error];
        payload.series_quality.estimated_labs = [ payload.series_quality.estimated_labs; {[]} ];
        
    end
    
    % como el ratio se usa para rankear mediciones, hacemos un ratio
    % ficticio para cada detección por el ranking de "parsimoniosidad" que
    % sacó  
    payload.series_quality.ratios = 1-(payload.series_quality.ratios * 1./max(payload.series_quality.ratios));
    
%     % estimate the mean RR sequence
%     [RRserie_resampled, time_resampled ] = resample_sequences(RRserie_filt);
%     
%     RRserie_med = median(RRserie_resampled,2);
% %     RRserie_var = var(RRserie_resampled,2);
%     RRserie_med_filt = colvec(MedianFiltSequence(time_resampled, RRserie_med, round(5 * ECG_header.freq)));
%     
%     RRserie_mean_sq_error = mean(bsxfun(@minus, RRserie_resampled, RRserie_med_filt).^2);
% 
%     [~, aux_idx] = sort(RRserie_mean_sq_error);
%     
%     lstable_rhythm_regions = size(stable_rhythm_regions,1);
%     
%     % como el ratio se usa para ranquear mediciones, hacemos un ratio
%     % ficticio para cada detección por el ranking de "parsimoniosidad" que
%     % sacó  
%     % contabilizo el first_guess como primera detección
%     aux_seq = linspace(1,0,lstable_rhythm_regions+1);
% 
%     % contabilizo el first_guess como primera detección
%     for ii = 1:(lstable_rhythm_regions+1)
%     
%         payload.series_quality.ratios(ii) = aux_seq(aux_idx(ii));
% 
%     end
    
% figure(3); plot(time_resampled, RRserie_resampled ); hold on; plot(time_resampled, RRserie_med, 'r', 'LineWidth', 2); hold off;

% 
%     progress_handle.checkpoint('Building final detection');
% 
%     final_idx = [];
%     WinSize = round(payload_in.final_build_time_win*ECG_header.freq);
% 
%     for ii = 1:WinSize:ECG_header.nsamp
% 
%         aux_idx = find(time_resampled >= ii & time_resampled <= (ii + WinSize));
%         
%         aux_idx2 = min_index(mean(bsxfun(@minus, RRserie_resampled(aux_idx,:), RRserie_med(aux_idx) ).^2));
% 
%         aux_val = RRserie_filt{aux_idx2};
%         final_idx = [final_idx; aux_val(aux_val(:,1) >= ii & aux_val(:,1) <= (ii + WinSize),1) ];
% 
%     end
% 
%     final_idx = unique(final_idx);
% 
%     % max local refinement in the signal
% 
%     [~, aux_idx, aux_idx2] = pack_signal(ECG_matrix(:,1), final_idx, [ pattern_prewin pattern_prewin ] );
% 
%     refined_max_idx = cellfun(@(a)( max_index( ECG_matrix(a,1) ) ), aux_idx);
% 
%     final_idx(aux_idx2) = final_idx(aux_idx2) - pattern_prewin + refined_max_idx ;
%     
%     %% save results
%     progress_handle.checkpoint('Saving results');
% 
%     % upsmple detections
%     str_aux = [ 'aip_auto' ];
%     payload.(str_aux).time = final_idx + ECG_start_offset - 1;
%     payload.series_quality.AnnNames = [ payload.series_quality.AnnNames ; {str_aux} {'time'} ];
%     payload.series_quality.ratios = [payload.series_quality.ratios; 0];
%     payload.series_quality.estimated_labs = [ payload.series_quality.estimated_labs; {[]} ];

    