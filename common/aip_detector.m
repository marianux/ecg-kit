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
        % Default auto
        payload_in.powerline_interference = nan; % Hz
    end
    
    if( ~isfield(payload_in, 'max_patterns_found') )
        % Default QRS width
        payload_in.max_patterns_found = 3; % # patterns
    end
    
    if( ~isfield(payload_in, 'sig_idx') )
        % Default QRS width
        payload_in.sig_idx = 1:ECG_header.nsig; % # patterns
    end

    %% preprocessing 
    
    % lead names desambiguation
    str_aux = regexprep(cellstr(ECG_header.desc), '\W*(\w+)\W*', '$1');
    lead_names = regexprep(str_aux, '\W', '_');

    [str_aux2, ~ , aux_idx] = unique( lead_names );
    aux_val = length(str_aux2);

    if( aux_val ~= ECG_header.nsig )
        for ii = 1:aux_val
            bAux = aux_idx==ii;
            aux_matches = sum(bAux);
            if( sum(bAux) > 1 )
                lead_names(bAux) = strcat( lead_names(bAux), repmat({'v'}, aux_matches,1), cellstr(num2str((1:aux_matches)')) );
            end
        end
    end

    lead_names = regexprep(lead_names, '\W*(\w+)\W*', '$1');
    lead_names = regexprep(lead_names, '\W', '_');
    
    
    progress_handle.checkpoint('Powerline filtering');

    ECG_matrix = double(ECG_matrix);

    if( isnan(payload_in.powerline_interference) )
        f_armonic = 1:3;
        f_armonic = f_armonic*50;
        f_armonic = f_armonic( f_armonic < (ECG_header.freq / 2) );
        f50_idx = round( f_armonic * ECG_header.nsamp / ECG_header.freq ) + 1;
        pwl50_e = abs(goertzel( ECG_matrix(:,payload_in.sig_idx), f50_idx ));
        
        f_armonic = 1:3;
        f_armonic = f_armonic*60;
        f_armonic = f_armonic( f_armonic < (ECG_header.freq / 2) );
        f60_idx = round( f_armonic * ECG_header.nsamp /ECG_header.freq ) + 1;
        pwl60_e = abs(goertzel( ECG_matrix(:,payload_in.sig_idx), f60_idx ));
        
        if( sum(sum(pwl50_e)) > sum(sum(pwl60_e)) )
            payload_in.powerline_interference = 50; % Hz
        else
            payload_in.powerline_interference = 60; % Hz
        end
    end
    
    % Resample for the closest smaller sampling freq
    target_fs = floor(ECG_header.freq/payload_in.powerline_interference) * payload_in.powerline_interference;
    
    if( target_fs == ECG_header.freq )
        ECG_matrix_res = ECG_matrix(:,payload_in.sig_idx);
    else
        [ups, downs] = rat(target_fs / ECG_header.freq);
        ECG_matrix_res = resample(ECG_matrix(:,payload_in.sig_idx), ups, downs );
    end

    pwl_size = round(target_fs / payload_in.powerline_interference);
    
    ECG_matrix_res = filter( ones(pwl_size,1)/pwl_size, 1, flipud(ECG_matrix_res) );
    ECG_matrix_res = filter( ones(pwl_size,1)/pwl_size, 1, flipud(ECG_matrix_res) );
    
    if( target_fs ~= ECG_header.freq )
        ECG_matrix_res = resample(ECG_matrix_res, downs, ups );
        aux_val = size(ECG_matrix_res,1);
        
        % Pad or trim after resampling
        if( aux_val < ECG_header.nsamp )
            ECG_matrix_res = [ECG_matrix_res; zeros( ECG_header.nsamp - aux_val, size(ECG_matrix_res,1) )];
        elseif( aux_val > ECG_header.nsamp )
            ECG_matrix_res = ECG_matrix_res(1:ECG_header.nsamp,:);
        end
    end
    
    ECG_matrix(:,payload_in.sig_idx) = ECG_matrix_res;
    
    %% Start
    pattern_size = 2*round(payload_in.trgt_width/2*ECG_header.freq)+1; % force odd number
    
    pattern_coeffs = diff(gausswin(pattern_size+1)) .* gausswin(pattern_size);
    first_pattern_coeffs = pattern_coeffs;
    
    progress_handle.Loops2Do = length(payload_in.sig_idx);

    lp_size = round(1.2*pattern_size);
    
    
    for this_sig_idx = rowvec(payload_in.sig_idx)
    
        progress_handle.start_loop();    

        pb_str_prefix = [ lead_names{this_sig_idx} ': ' ];
        %% Initial Pattern search 

        progress_handle.checkpoint([ pb_str_prefix 'First pattern guess' ]);
        
        rise_detector = filter( first_pattern_coeffs, 1, flipud(ECG_matrix(:,this_sig_idx)) );
        rise_detector = filter( first_pattern_coeffs, 1, flipud(rise_detector) );
        rise_detector = filter( ones(lp_size,1)/lp_size, 1, flipud(abs(rise_detector)) );
        rise_detector = filter( ones(lp_size,1)/lp_size, 1, flipud(rise_detector) );
        
        % figure(4); plot(rise_detector ./ max(rise_detector)); hold on; plot(ECG_matrix(:,this_sig_idx) ./ max(ECG_matrix(:,this_sig_idx)))

        % figure(3); plot(ECG_matrix(:,1) ); xlims = xlim(); ylims = ylim(); box_h = ylims + [0.01 -0.1] * diff(ylims); aux_sc = 0.8*diff(ylims)/(max(rise_detector)-min(rise_detector));aux_off = ylims(1) + 0.6*diff(ylims); hold on; plot(rise_detector * aux_sc + aux_off); hold off; ylim(ylims);

        progress_handle.checkpoint([ pb_str_prefix 'Peak detection first guess' ]);

        actual_thr = thr_calc(rise_detector, payload_in.trgt_min_pattern_separation, ECG_header.freq);
        
        if isnan(actual_thr)
            cprintf('Red', '\nRecording %s, lead %s, no patterns above threshold found.(Current %f seconds). \n\n',ECG_header.recname , lead_names{this_sig_idx}, payload_in.stable_RR_time_win);
            continue
        end
        
        first_detection_idx = modmax(rise_detector, 1, actual_thr, 1, round(payload_in.trgt_min_pattern_separation * ECG_header.freq));
         
        %% Look for stable segments
        if numel(first_detection_idx) <= 1
            cprintf('Red', '\nRecording %s, lead %s, Not enough patterns found over noise floor. \n\n',ECG_header.recname , lead_names{this_sig_idx});                    
            continue
        end
        
        RRserie = RR_calculation(first_detection_idx, ECG_header.freq);               

        RRserie_filt = MedianFiltSequence(first_detection_idx, RRserie, round(5 * ECG_header.freq));

        RR_scatter = abs(RRserie - RRserie_filt)./RRserie_filt;

        % Find the most stable regions
        RR_thr = prctile(RR_scatter, 50);

        stable_rhythm_regions = get_segments_from_sequence(first_detection_idx, RR_scatter, round(payload_in.stable_RR_time_win * ECG_header.freq), RR_thr );

%     figure(3); plot(first_detection_idx, [RRserie RRserie_filt] ); ylims = ylim(); ylims = ylims + [0.1 -0.1] * diff(ylims); hold on; plot( [ stable_rhythm_regions'; flipud(stable_rhythm_regions'); stable_rhythm_regions(:,1)' ],  repmat([ylims(1); ylims(1); ylims(2);ylims(2);ylims(1)],1, size(stable_rhythm_regions,1) ), 'b:' ); hold off;
%     figure(3); plot(first_detection_idx, [RR_scatter], ':xb' ); xlims = xlim(); ylims = ylim(); ylims = ylims + [0.01 -0.1] * diff(ylims); hold on; plot(xlims, [RR_thr RR_thr], '--r'); plot( [ stable_rhythm_regions'; flipud(stable_rhythm_regions'); stable_rhythm_regions(:,1)' ],  repmat([ylims(1); ylims(1); ylims(2);ylims(2);ylims(1)],1, size(stable_rhythm_regions,1) ), 'b:' ); hold off;
%     figure(3); plot(first_detection_idx, [RRserie RRserie_filt] ); xlims = xlim(); ylims = ylim(); ylims = ylims + [0.1 -0.1] * diff(ylims); aux_sc = 0.25*diff(ylims)/(max(RR_scatter)-min(RR_scatter));aux_off = ylims(1) + 0.65*diff(ylims); hold on; plot(first_detection_idx, (RR_scatter*aux_sc) +aux_off , ':xb' ); plot(xlims, ([RR_thr RR_thr]*aux_sc ) + aux_off, '--r'); plot( [ stable_rhythm_regions'; flipud(stable_rhythm_regions'); stable_rhythm_regions(:,1)' ],  repmat([ylims(1); ylims(1); ylims(2);ylims(2);ylims(1)],1, size(stable_rhythm_regions,1) ), 'b:' ); hold off;    

        if( isempty(stable_rhythm_regions) )
            cprintf('Red', '\nRecording %s, lead %s, no stable rhythm segments found. Consider decreasing "stable_RR_time_win" (Current %f seconds). \n\n',ECG_header.recname , lead_names{this_sig_idx}, payload_in.stable_RR_time_win);                    
            continue
        end

        lstable_rhythm_regions = size(stable_rhythm_regions,1);

        % Get the larger or more stable segments first
        [~, aux_idx] = sort(diff(stable_rhythm_regions,1,2), 'descend');

        larger_stable_rank = arrayfun(@(a)(find(a==aux_idx)), 1:lstable_rhythm_regions);

        pattern_prewin = round(payload_in.trgt_width/2*ECG_header.freq); % force odd number

        % Get the segments with less change in morphology
        pack_variance = repmat(realmax,lstable_rhythm_regions,1);

        % From the longest - most rhythm stable regions
        for ii = 1: round(lstable_rhythm_regions/4)

            jj = aux_idx(ii);
            
            % Woody method refinement
            aux_idx2 = find(first_detection_idx >= stable_rhythm_regions(jj,1) & first_detection_idx <= stable_rhythm_regions(jj,2));
            avg_pack = pack_signal(ECG_matrix(:,1), first_detection_idx(aux_idx2), [ pattern_prewin pattern_prewin ], true);    

            % Mean variance across the ensemble
            pack_variance(jj) = nanmean(nanvar(squeeze(avg_pack), [], 2));

        end
        [~, aux_idx] = sort(pack_variance);

        morphology_rank = arrayfun(@(a)(find(a==aux_idx)), 1:lstable_rhythm_regions);

        % Final rank
        [~, aux_idx] = sort(larger_stable_rank + morphology_rank);

        stable_rhythm_regions = stable_rhythm_regions(aux_idx( 1:min(payload_in.max_patterns_found , size(stable_rhythm_regions,1))), : );

        % figure(3); plot(first_detection_idx, [RRserie RRserie_filt] ); xlims = xlim(); ylims = ylim(); ylims = ylims + [0.1 -0.1] * diff(ylims); aux_sc = 0.25*diff(ylims)/(max(RR_scatter)-min(RR_scatter));aux_off = ylims(1) + 0.65*diff(ylims); hold on; plot(first_detection_idx, (RR_scatter*aux_sc) +aux_off , ':xb' ); plot(xlims, ([RR_thr RR_thr]*aux_sc ) + aux_off, '--r'); plot( [ stable_rhythm_regions'; flipud(stable_rhythm_regions'); stable_rhythm_regions(:,1)' ],  repmat([ylims(1); ylims(1); ylims(2);ylims(2);ylims(1)],1, size(stable_rhythm_regions,1) ), 'b:' ); arrayfun(@(a)(text( stable_rhythm_regions(a,1), ylims(2), num2str(a))), 1:min(payload_in.max_patterns_found )); hold off;

        %% save first detections as a separate time serie

        str_aux = [ 'aip_guess' '_' lead_names{this_sig_idx} ];
        payload.(str_aux).time = first_detection_idx + ECG_start_offset - 1;
        payload.series_quality.AnnNames = [ payload.series_quality.AnnNames ; {str_aux} {'time'} ];
        RRserie_mean_sq_error = mean((RRserie - RRserie_filt).^2);
        payload.series_quality.ratios = [payload.series_quality.ratios; RRserie_mean_sq_error];
        payload.series_quality.estimated_labs = [ payload.series_quality.estimated_labs; {[]} ];

        k_woody = 3;
        win_limits = [ pattern_prewin (pattern_prewin+1) ];
        extended_win_limits = [k_woody * win_limits(1), k_woody * win_limits(1)+1];
        aux_offset = (extended_win_limits(1) - win_limits(1) + 1); % ???
        
        aux_pat_coeff = zeros(sum(win_limits),size(stable_rhythm_regions,1));
        
        for ii = 1:size(stable_rhythm_regions,1)

            % refinamos con el metodo de woody
            aux_idx2 = find(first_detection_idx >= stable_rhythm_regions(ii,1) & first_detection_idx <= stable_rhythm_regions(ii,2));

    % figure(3); avg_pack = pack_signal(ECG_matrix(:,this_sig_idx), first_detection_idx(aux_idx2), 10*[ pattern_prewin pattern_prewin ], true); plot_ecg_mosaic(avg_pack)

            if isempty(aux_idx2)
                % Cannot found any stable regions to apply the woody_method
                cprintf('Red', '\nRecording %s, lead %s, Cannot found any stable regions. \n\n',ECG_header.recname , lead_names{this_sig_idx});
                continue
            end
    
            pattern_coeffs = woody_method(ECG_matrix(:,this_sig_idx), first_detection_idx(aux_idx2), extended_win_limits, 0.95);
            % trim the pattern 
            aux_pat_coeff(:,ii) = pattern_coeffs( aux_offset:(aux_offset+size(aux_pat_coeff,1)-1) );
            
            % figure(5); plot(first_pattern_coeffs ./ max(abs(first_pattern_coeffs))); hold on; plot(aux_pat_coeff(:,ii) ./ max(abs(aux_pat_coeff(:,ii))))

            progress_handle.checkpoint([ pb_str_prefix 'Pattern match ' num2str(ii) ]);

        %         rise_detector = filter(ones(lp_size,1)/lp_size,1, abs([0; diff(double(ECG_matrix(:,ss)))]));
        %         rise_detector = [zeros(((lp_size-1)/2)-1,1); rise_detector((lp_size-1)/2:end)];
            rise_detector = filter( pattern_coeffs, 1, flipud(ECG_matrix(:,this_sig_idx)) );
            rise_detector = filter( pattern_coeffs, 1, flipud(rise_detector) );
            rise_detector = filter( ones(lp_size,1)/lp_size, 1, flipud(abs(rise_detector)) );
            rise_detector = filter( ones(lp_size,1)/lp_size, 1, flipud(rise_detector) );
            
            % figure(6); plot(rise_detector ./ max(rise_detector)); hold on; plot(ECG_matrix(:,this_sig_idx) ./ max(ECG_matrix(:,this_sig_idx)))

            progress_handle.checkpoint([ pb_str_prefix 'Peak detection pattern ' num2str(ii) ]);

            actual_thr = thr_calc(rise_detector, payload_in.trgt_min_pattern_separation, ECG_header.freq);
            
            if isnan(actual_thr)
                cprintf('Red', '\nRecording %s, lead %s, no patterns above threshold found.(Current %f seconds). \n\n',ECG_header.recname , lead_names{this_sig_idx}, payload_in.stable_RR_time_win);
                continue
            end
            
            % Here the MIN time restriction is important
            this_idx = modmax(rise_detector, 1, actual_thr, 1, round(payload_in.trgt_min_pattern_separation * ECG_header.freq));
            
            if numel(first_detection_idx) <= 1
                cprintf('Red', '\nRecording %s, lead %s, Not enough patterns found over noise floor. \n\n',ECG_header.recname , lead_names{this_sig_idx});                    
                continue
            end
            
            RRserie = RR_calculation(this_idx, ECG_header.freq);               

            this_RRserie_filt = colvec(MedianFiltSequence(this_idx, RRserie, round(5 * ECG_header.freq)));

            RRserie_mean_sq_error = mean((RRserie - this_RRserie_filt).^2);       

            % Structure of detection loading
            str_aux = [ 'aip_patt_' num2str(ii) '_' lead_names{this_sig_idx} ];
            payload.(str_aux).time = this_idx + ECG_start_offset - 1;
            payload.series_quality.AnnNames = [ payload.series_quality.AnnNames ; {str_aux} {'time'} ];
            payload.series_quality.ratios = [payload.series_quality.ratios; RRserie_mean_sq_error];
            payload.series_quality.estimated_labs = [ payload.series_quality.estimated_labs; {[]} ];

        end

        % Como el ratio se usa para rankear mediciones, hacemos un ratio
        % ficticio para cada detección por el ranking de "parsimoniosidad" que
        % sacó  
        payload.series_quality.ratios = 1-(payload.series_quality.ratios * 1./max(payload.series_quality.ratios));

        progress_handle.end_loop();    

    end
end
    
function [actual_thr] = thr_calc(rise_detector, trgt_min_pattern_separation, freq)
    initial_thr = 30;

    actual_thr = prctile(rise_detector, initial_thr);

    [~, max_values ] = modmax(rise_detector, 1, actual_thr, 1, round(trgt_min_pattern_separation * freq));
    
    prctile_grid = prctile(max_values,1:100);

	grid_step = median(diff(prctile_grid));
	
    if grid_step == 0
		actual_thr = nan;
    else
        thr_grid = actual_thr:grid_step:max(max_values);
        
        hist_max_values = histcounts(max_values, thr_grid);

        first_bin_idx = 1;

        [thr_idx, thr_max ] = modmax( colvec( hist_max_values ) , first_bin_idx, 0, 0, [], [] );

        if( isempty(thr_max) )
            actual_thr = nan;
        else
            % Mass center of the distribution, probably the value where the
            % patterns under search are located.
            thr_idx_expected = floor(rowvec(thr_idx) * colvec(thr_max) *1/sum(thr_max));

            aux_seq = 1:length(thr_grid);

            % MPS0018 recording from Basel VIII database forced to add the
            % lesser or equal to aux_seq <= thr_idx_expected. Probably a
            % recording where no clear maximum above the noise threshold were
            % found.
			min_hist_max_values = min( hist_max_values( aux_seq >= first_bin_idx & aux_seq < thr_idx_expected) );
			
            % In case several indexes match the boolean condition, the mean
            % index is the center of all those indexes. Other criteria such as
            % min or max can be explored
            % thr_min_idx = round(mean(find(bAux & [hist_max_values 0] == min_hist_max_values)));
			thr_min_idx = round(mean(find(aux_seq >= first_bin_idx & aux_seq < thr_idx_expected & [hist_max_values 0] == min_hist_max_values)));
            
            actual_thr = thr_grid(thr_min_idx);
        end
    end
end
