function IdentField = getident_multids(dsAny, strIdent, path2ds )

if( nargin < 3 )
    path2ds = '.';
end

dataset_user_data = getuser(dsAny);

if( isfield(dataset_user_data, 'DSnames') )
    %dataset multiparte

    db_id = dataset_user_data.db_id;
    
    ds2work = cellstr(dataset_user_data.DSnames);

    IdentField = [];

    for this_ds = rowvec(ds2work)

        filename = [ path2ds filesep char(this_ds) '.mat' ];
        
        if( exist(filename, 'file') )
            %por si me quedo sin mem
            load(filename);

            dataset_user_data = getuser(dsTarget_dataset);

            if( db_id ~= dataset_user_data.db_id)
                error('MEC:NotPartOfMultipartDs', [ filename ' no es una parte válida. Rearme los datasets con DoFeatureMatricesFromFiles'] )
            end
            
            IdentField = [IdentField; getident(dsTarget_dataset, strIdent )];
            
        else
            error('MEC:PartNotFound', [ 'No encontramos la parte ' filename])
        end

    end

else
    %si no es multiparte, lo tomo como leido correctamente.
    IdentField = getident(dsAny, strIdent );
end
