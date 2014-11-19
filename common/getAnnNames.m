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
    