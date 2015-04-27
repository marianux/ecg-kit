%% (Internal) Convert matrix of ECG wave annotations to a struct position format, used in wavedet algorithm
%   
%   [this_position, fields ] = matrix2positions(pos_matrix, fields)
% 
% Arguments:
% 
%      + pos_matrix: with the fields defined in fields variable
% 
%      + fields: field names in pos_matrix
%             
% Output:
% 
%      + this_position: struct with fields vectors of annotations
% 
%      + fields: field names in pos_matrix
% 
% Example:
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
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
