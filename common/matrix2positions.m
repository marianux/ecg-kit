%% (Internal) Convert matrix of ECG wave annotations to a struct position format, used in wavedet algorithm
%   
%   this_position = matrix2positions(pos_matrix, fields)
% 
% Arguments:
% 
%      + pos_matrix: cell with the fields defined in fields variable
% 
% Output:
% 
%      + this_position: struct with fields vectors of annotations
% 
%      + fields: field names in pos_matrix
% 
% Example:
% 
% See also positions2matrix
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/05/2014
% Birthdate  : 18/05/2016
% Copyright 2008-2016
% 
function this_position = matrix2positions(pos_matrix, fields)

    this_position = [];
    
    if(~iscell(pos_matrix)) 
        pos_matrix = {pos_matrix};
    end
    
    fields = rowvec(cellstr(fields));

    cant_leads = length(pos_matrix);
    cant_fields = length(fields);
    
    for ii = 1:cant_leads
    
        aux_matrix = pos_matrix{ii};
        aux_struct = [];
        
        for count = 1:cant_fields
            aux_struct.(fields{count}) = colvec(aux_matrix(:,count));
        end       
        
        if( isempty(this_position) )
            this_position = aux_struct; 
        else
            this_position(ii) = aux_struct; 
        end

    end
    
end
