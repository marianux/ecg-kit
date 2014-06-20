function [this_position, fields ] = matrix2positions(pos_matrix, fields)

    fields = rowvec(cell2str(fields));

    cant_leads = length(pos_matrix);
    
    for ii = 1:cant_leads
    
        aux_matrix = pos_matrix{ii};
        aux_struct = [];
        
        count = 1;
        for fn = fields
            fn = fn{1};
            aux_struct.(fn) = colvec(aux_matrix(:,count));
            count = count + 1;
        end       
        
        if( isempty(this_position) )
            this_position = aux_struct; 
        else
            this_position(ii) = aux_struct; 
        end

    end
    
end
