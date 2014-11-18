function all_detections = calculate_artificial_QRS_detections(all_detections, ECG_header, start_end_this_segment)

    if( nargin < 3 || isempty(start_end_this_segment) ) 
        start_end_this_segment = [1 ECG_header.nsamp];
    end
    % attemp to build a better detection from single-lead detections.

    [~, all_annotations] = getAnnNames(all_detections);

    [ ratios, estimated_labs ] = CalcRRserieRatio(all_annotations, ECG_header, start_end_this_segment);

    [~, best_detections_idx] = sort(ratios, 'descend');

    % generate artificial annotations combining K best annotations
    aux_idx = best_detections_idx(1:min(10, length(best_detections_idx)) );
    artificial_annotations = combine_anns(all_annotations(aux_idx), estimated_labs(aux_idx), ECG_header );

    for ii = 1:length(artificial_annotations)
        aux_str = ['mixartif_ECGmix' num2str(ii)];
        all_detections.(aux_str) = artificial_annotations(ii);
    end

    [AnnNames, all_annotations] = getAnnNames(all_detections);

    [ ratios, estimated_labs] = CalcRRserieRatio(all_annotations, ECG_header, start_end_this_segment);

    [ratios, best_detections_idx] = sort(ratios, 'descend');

    all_detections.series_quality.ratios = ratios;
    all_detections.series_quality.estimated_labs = estimated_labs;
    all_detections.series_quality.AnnNames = AnnNames(best_detections_idx,:); %#ok<STRNU>



