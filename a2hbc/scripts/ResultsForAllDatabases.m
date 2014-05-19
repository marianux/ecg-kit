function ResultsForAllDatabases( tmp_path, op_modes )

if( ~iscell(op_modes) )
    op_modes = {op_modes};
end

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

tmp_path_files = dir(tmp_path);
dirs_idx = find(cell2mat({tmp_path_files(:).isdir}));
all_global_performances = [];

try
    
    for dir_idx = rowvec(dirs_idx)

        bDataFound = false;

        for op_mode_idx = 1:length(op_modes)
            op_mode = op_modes{op_mode_idx};
            DB_name = tmp_path_files(dir_idx).name;

            result_files = dir([tmp_path DB_name filesep 'tmpfile_results_' DB_name '_' op_mode '*.mat' ]);

            if( ~isempty(result_files) )

                bDataFound = true;

                result_filenames = {result_files(:).name};

                tokens = regexp(result_filenames, [ '.*' DB_name '_' op_mode '_(.+).mat' ], 'tokens');
                ltokens = length(tokens);
                configs = [];

                for ii = 1:ltokens
                    aux = tokens{ii};
                    if( ~isempty(aux) )
                        configs = [configs ; cellstr(aux{1})];
                    end
                end
                
                all_configs = unique( configs );                 

                if(strcmp(op_mode, 'auto' ))
                    all_configs = flipud( all_configs );                 
                end
                
                for config = rowvec(all_configs)

                    title_str = ['## ' DB_name ' ' op_mode '_' config{1} ' ##\n'];
                    str_length = length(title_str)-2;
                    fprintf(1, [ '\n\n' ...
                                repmat('#',1,str_length) '\n' ...
                                title_str ...
                                repmat('#',1,str_length) '\n'] );

                    load( [tmp_path DB_name filesep 'tmpfile_results_' DB_name  '_' op_mode '_' config{1} '.mat' ] );

                    fprintf(1, 'Encontramos %3.2f iteraciones.\n', mean(sum(~isnan(squeeze(ConfusionMatrix(1,1,:,:))))));
                    
                    [~,~,~, global_performances] = DisplayResults('dsResult', nansum(ConfusionMatrix,4), 'datasetName', DB_name, 'ClassLabels', Lablist);

                    aux = sprintf( [ '%s\t%s\t%s\t' repmat('%1.0fXX%1.0f\t', 1,  size(global_performances,2)/2 )], DB_name, op_mode, config{1}, global_performances(2,:));
                    all_global_performances = [all_global_performances; cellstr(aux(1:end-1)) ];
                end
                
               
            else
                if( bDataFound )

                    title_str = ['## ' DB_name ' ' op_mode ' ##\n'];
                    str_length = length(title_str)-2;
                    fprintf(1, [ '\n\n' ...
                                repmat('#',1,str_length) '\n' ...
                                title_str ...
                                repmat('#',1,str_length) '\n'] );

                    fprintf(2, '\n\n-Data not found-\n\n');

                end
            end
            
        end
        
        if( bDataFound )
            all_global_performances = [all_global_performances; {' '}  ];
        end
        
    end

    fprintf(1,'Resumen de performances globales\n');
    fprintf(1,'--------------------------------\n\n');
    for ii = 1:length(all_global_performances)
        fprintf(1,'%s\n', all_global_performances{ii} );
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

    fprintf(2,'Path: %s Opmode:%s\n', tmp_path, op_modes);

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
