function Class_indexes = getnlab_multids(dsAny, path2ds)

if( nargin < 2 )
    path2ds = '.';
end

dataset_user_data = getuser(dsAny);

if( isfield(dataset_user_data, 'DSnames') )
    %dataset multiparte

    db_id = dataset_user_data.db_id;
    
    ds2work = cellstr(dataset_user_data.DSnames);

    Class_indexes = [];
    for this_ds = rowvec(ds2work)

        filename = [ path2ds filesep char(this_ds) '.mat' ];
        
        if( exist(filename, 'file') )
            %por si me quedo sin mem
            load(filename);

            dataset_user_data = getuser(dsTarget_dataset);

            if( db_id ~= dataset_user_data.db_id)
                error('MEC:NotPartOfMultipartDs', [ filename ' no es una parte válida. Rearme los datasets con DoFeatureMatricesFromFiles'] )
            end
            
            Class_indexes = [Class_indexes; getnlab(dsTarget_dataset)];
            
        else
            error('MEC:PartNotFound', [ 'No encontramos la parte ' filename])
        end

    end

else
    %si no es multiparte, lo tomo como leido correctamente.
    Class_indexes = getnlab(dsAny);
end
