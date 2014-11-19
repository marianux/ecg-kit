function all_detections = wavedet_QRS_detection_mix(struct_in, ECG_header, start_end_this_segment )

    % attemp to build a better detection from single-lead detections.

    AnnNames = struct_in.series_quality.AnnNames(:,1);
    detector_name = cellfun( @(a)( strtok(a, '_')), AnnNames, 'UniformOutput', false);
    
    aux_idx = find(strcmpi(detector_name, 'wavedet'));
    
    all_annotations = {};
    for ii = rowvec(aux_idx)
        all_annotations = [all_annotations; {struct_in.(struct_in.series_quality.AnnNames{ii,1}).time}];
    end

    [ ratios, estimated_labs ] = CalcRRserieRatio(all_annotations, ECG_header, start_end_this_segment);

    [~, best_detections_idx] = sort(ratios, 'descend');

    % generate artificial annotations combining K best annotations
    aux_idx = best_detections_idx(1:min(10, length(best_detections_idx)) );
    artificial_annotations = combine_anns(all_annotations(aux_idx), estimated_labs(aux_idx), ECG_header );

    for ii = 1:length(artificial_annotations)
        aux_str = ['wavedetMix_ECGmix' num2str(ii)];
        all_detections.(aux_str) = artificial_annotations(ii);
    end



