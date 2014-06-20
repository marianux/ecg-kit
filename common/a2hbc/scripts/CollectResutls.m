function CollectResutls( tmp_path, op_mode )

if( tmp_path(end) ~= filesep )
    tmp_path = [tmp_path filesep];
end

%Java user interface is started. Not started in clusters for example.
bHaveUserInterface = usejava('desktop');

%path related constants.
mylocation_path = [fileparts(mfilename('fullpath')) filesep ];

default_paths = { ...
                    [ mylocation_path 'common' filesep ';' ]; ...
                    [ mylocation_path 'common' filesep 'prtools' filesep ';' ]; ...
                    [ mylocation_path 'common' filesep 'prtools_addins' filesep ';' ]; ...
                    [ mylocation_path 'common' filesep 'kur' filesep ';' ]; ...
                    [ mylocation_path 'common' filesep 'LIBRA' filesep ';' ]; ...
                };
            
default_paths = char(default_paths)';
default_paths = (default_paths(:))';
addpath(default_paths);

Cleanup_hdl = onCleanup(@()(rmpath(default_paths)));

DB_name = find(tmp_path == filesep, 2, 'last' );
DB_name = tmp_path(DB_name(1)+1:DB_name(2)-1);

result_files = dir([tmp_path 'tmpfile_*_results_' op_mode '*.mat' ]);
result_filenames = {result_files(:).name};

tokens = regexp(result_filenames, [ 'tmpfile_(.+)_results_' op_mode '_(.+)_(\d+).mat' ], 'tokens');

if( isempty(tokens) )
    error('Results not found');
end

ltokens = length(tokens);
recordings = [];
configs = [];
iterations = [];

for ii = 1:ltokens
    aux = tokens{ii};
    if( ~isempty(aux) )
        aux = aux{1};
        recordings = [recordings ; cellstr(aux{1})];
        configs = [configs ; cellstr(aux{2})];
        iterations = [iterations ; cellstr(aux{3})];
    end
end

all_recordings = unique( recordings );
all_configs = unique( configs );
cant_iteraciones = max( str2double(iterations) );
cant_registros = length(all_recordings);

try
    
    for config = rowvec(all_configs)
        
        file_skiped = 0;

        ConfusionMatrix = [];

        fprintf(1, [ 'Processing ' DB_name ' config ' config{1} '\n\n']);
    %     bSkipConfig = false;

        for iter = 1:cant_iteraciones 

            result_files = dir([tmp_path 'tmpfile_*_results_' op_mode '_' config{1} '_' num2str(iter) '.mat' ]);

            if( length(result_files) ~= cant_registros )
                if(usejava('desktop'))
    %                 bSkipConfig = true;
                    fprintf(2, [ 'Not all files processed for ' DB_name ' config ' config{1} ' iter ' num2str(iter) '\n\n']);
                    break
                else
                    error([ 'Not all files processed for ' DB_name ' config ' config{1} ' iter ' num2str(iter) ]);
                end
            end

            %intento leer algo, lo primero que encuentre
            for ii = 1:cant_registros
                try
                    aux = load( [tmp_path result_files(ii).name]);
                    break
                catch ME
                    continue;
                end
            end
            
            %Hago una matriz de NaNs para poder ignorar las iteraciones que
            %no se hayan calculado.
            cm = nan(size(aux.LabelList,1), size(aux.LabelList,1), size(aux.ConfusionMatrix,3) , cant_registros);
            
            for ii = 1:cant_registros

                if( result_files(ii).bytes == 0)
                    fprintf(2,['Skiping file ' tmp_path result_files(ii).name ' cause it is empty !\n' ]);
                    file_skiped = file_skiped + 1;
                    continue;
                end
                
                aux = load( [tmp_path result_files(ii).name]);

                if( ii == 1)
                    %, 'Labels', 'ConfusionMatrix', 'LabelList'
                    orig_lablist = aux.LabelList;
                    orig_lablist_size = size(orig_lablist,1);

                    if( iter == 1)
%                         cm = aux.ConfusionMatrix;
                    else
                        aux_intersect = intersect(cellstr(orig_lablist), cellstr(other_iters_LabelList));
                        laux_intersect = length(aux_intersect);

                        if(laux_intersect == orig_lablist_size && laux_intersect == other_iters_LabelList_size )
                            cm = aux.ConfusionMatrix;
                        else
                            error('Different lablists. Check results.');    
                        end
                    end
                else
                    aux_intersect = intersect(cellstr(aux.LabelList), cellstr(orig_lablist));
                    laux_intersect = length(aux_intersect);
                    this_lablist_size = size(aux.LabelList,1);

                    if(laux_intersect == this_lablist_size && laux_intersect == orig_lablist_size )
%                         cm(:,:,:,ii) = aux.ConfusionMatrix;
                    else
                        error('Different lablists. Check results.');    
                    end

                end

                cm(:,:,:,ii) = aux.ConfusionMatrix;

            end

            ConfusionMatrix = cat(3, ConfusionMatrix, cm);

            other_iters_LabelList = orig_lablist;
            other_iters_LabelList_size = orig_lablist_size;

        end
        
        if( file_skiped > 0)
            fprintf(2, 'Skiping %d files for %s configuration.\n', file_skiped, [DB_name ' ' config{1}] );
        end

        Lablist = orig_lablist;

        save([tmp_path 'tmp_' DB_name '_results_' op_mode '_' config{1} '.mat' ], 'ConfusionMatrix', 'Lablist');

        % delete([tmp_path 'tmpfile_*_results_' op_mode '_' config{1} '_*.mat' ]);

        movefile([tmp_path 'tmp_' DB_name '_results_' op_mode '_' config{1} '.mat' ], [tmp_path 'tmpfile_results_' DB_name '_' op_mode '_' config{1} '.mat' ], 'f' );

    end
    
    if( ~bHaveUserInterface )
        % flag that the program ended correctly
        setenv('A2HBC_ES', '0');
    end

    

catch MException
    
    fprintf(2,'\n\n')
    fprintf(2,'###########\n')
    fprintf(2,'## ERROR ##\n')
    fprintf(2,'###########\n')

    fprintf(2,'Path: %s Opmode:%s\n', tmp_path, op_mode);

    local_host = getenv('HOSTNAME');
    computer_arch = computer();

    fprintf(2,'Computer: %s (%s) \n', local_host, computer_arch);

    report = getReport(MException);
    fprintf(2, '%s', report);
    fprintf(2,'###########\n')
    fprintf(2,'## ERROR ##\n')
    fprintf(2,'###########\n')

    rethrow(MException)
            
end
