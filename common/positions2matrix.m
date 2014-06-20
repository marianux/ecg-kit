function [pos_matrix, fields ] = positions2matrix(this_position, fields)

    if( nargin < 2 || isempty(fields) )
        fields = rowvec(fieldnames(this_position));
    else
        if( ~iscell(fields) )
            fields = rowvec(cell2str(fields));
        end
    end

    cant_leads = length(this_position);
    pos_matrix = cell(cant_leads,1);
    
    for ii = 1:cant_leads
        
        % get the field with max num of annotations
        aux_maxlength = 0;
        for fname = fields
            aux_maxlength = max(aux_maxlength, length(this_position(ii).(fname{1})));
        end

        %store to use rules on all annotations at the end
        aux_matrix = nan(aux_maxlength,length(fields));

        count = 1;
        for fn = fields
            fn = fn{1};
            if( isfield( this_position(ii), fn) )
                aux_matrix(1:length(this_position(ii).(fn)),count) = colvec(this_position(ii).(fn));
            end
            count = count + 1;
        end
        
        pos_matrix{ii} = aux_matrix;
        
    end
    
    if( cant_leads == 1 )
        pos_matrix = aux_matrix;
    end
    
end
