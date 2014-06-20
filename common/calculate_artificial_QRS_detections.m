function all_detections = calculate_artificial_QRS_detections(all_detections, ECG_header)

    % attemp to build a better detection from single-lead detections.

    [AnnNames, all_annotations] = getAnnNames(all_detections);

    cant_anns = size(AnnNames,1);

    [ ratios, estimated_labs ] = CalcRRserieRatio(all_annotations, ECG_header);

    [~, best_detections_idx] = sort(ratios, 'descend');

    % generate artificial annotations combining K best annotations
    aux_idx = best_detections_idx(1:min(10, length(best_detections_idx)) );
    artificial_annotations = combine_anns(all_annotations(aux_idx), estimated_labs(aux_idx), ECG_header );

    for ii = 1:length(artificial_annotations)
        aux_str = ['artificial_' num2str(ii)];
        all_detections.(aux_str) = artificial_annotations(ii);
    end

    [AnnNames, all_annotations] = getAnnNames(all_detections);

    [ ratios, estimated_labs] = CalcRRserieRatio(all_annotations, ECG_header);

    [ratios, best_detections_idx] = sort(ratios, 'descend');

    all_detections.series_quality.ratios = ratios;
    all_detections.series_quality.estimated_labs = estimated_labs;

    all_detections.series_quality.all_annotations = all_annotations(best_detections_idx);

    all_detections.series_quality.AnnNames = AnnNames(best_detections_idx,:); %#ok<STRNU>


function [AnnNames, all_annotations] = getAnnNames(aux_struct)

    AnnNames = [];

    for fname = rowvec(fieldnames(aux_struct))
        if( isfield(aux_struct.(fname{1}), 'time') )
            AnnNames = [AnnNames; cellstr(fname{1}) cellstr('time')];
        end
        if( isfield(aux_struct.(fname{1}), 'qrs') )
            AnnNames = [AnnNames; cellstr(fname{1}) cellstr('qrs')];
        end
    end

    cant_anns = size(AnnNames,1);

    all_annotations = cell(cant_anns,1);
    for ii = 1:cant_anns
        all_annotations{ii} = aux_struct.(AnnNames{ii,1}).(AnnNames{ii,2});
    end
