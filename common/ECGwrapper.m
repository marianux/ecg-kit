%% Allow acces to ECG recordings of arbitrary format and length.
% 
% Description:
% ECG wrapper is a class to allow the access to cardiovascular signal
% recordings of several formats (MIT, ISHNE, AHA, HES, MAT) and lengths,
% from minutes to days. Also it can be plugged to an ECGtask object to
% perform several tasks, such as QRS detection and ECG delineation among
% others. 
% 
% Arguments: (specified as ECGwrapper('arg_name1', arg_val1, ... , 'arg_nameN', arg_valN) )
% 
%           + recording_name : (char) the full filename of the ECG
%                recording. 
% 
%           + recording_format : (char) the format of the ECG recording. By
%                default or if not specified, the wrapper will attemp to
%                auto-detect the format.
% 
%           + this_pid : (char) In case working in a multiprocess
%                environment, this value will identify the current process.
%                Can be a numeric value, or a string of the form 'N/M'.
%                This pid is N and the total amount of pid's to divide the
%                whole work is M.
% 
%           + tmp_path: (char) path to store the temp files. By default
%                will be the output of tempdir function.
% 
%           + output_path: (char) path to store the result files. By default
%                will be the same path of the recordings.
% 
%           + ECGtaskHandle: (char || ECGtask) The task to perform, can be
%                the name of the task, or an ECGtask object. Available
%                ECGtasks can be listed with list_all_ECGtask() command. 
% 
%           + overlapping_time: (numeric) Time in seconds of overlapp among
%                consequtive segments. This segment is useful for ensuring
%                transitory responses of systems to be finished. Default is
%                30 seconds.
% 
%           + partition_mode: (numeric) The way that this object will
%                partition lengthy signals:
%              - 'ECG_contiguous' no overlapp between segments
%              - 'ECG_overlapped' overlapp 'overlapping_time' between
%                 segments. Default.
%              - 'QRS' do the partition based on the annotations in
%                 ECG_annotations.time property. Typically but not
%                 necessary are QRS annotations.
% 
%           + cacheResults: (logical) To save intermediate results to
%             recover in case of failure. Default is true.
% 
%           + syncSlavesWithMaster: (logical) In multiprocess environments
%             sometimes it is useful to terminate all pid's together in
%             orther to do further tasks later. Default is false.
% 
%           + repetitions: (numeric) In case the ECGtask is not
%             deterministic, the repetition property allows to repeat the
%             task several times. Default value is 1.
% 
% 
% Output:
% The constructed object.
% 
% Examples:
% 
%       you can also check 'examples.m' in the same folder of this file.
% 
% 
% See also ECGtask
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 20/3/2012
% Copyright 2008-2015
% 
classdef ECGwrapper < handle
       
    properties(GetAccess = private, Constant)
        
        % all ECG formats that ecg-kit can handle
        cKnownFormats = {'MIT' 'ISHNE', 'AHA', 'HES', 'MAT', 'Mortara', 'auto'};
        % The fields required in the ECGheader property
        cHeaderFieldNamesRequired = {'freq' 'nsamp' 'nsig' 'gain' 'adczero'};
        % Possible Partitions modes 
        cPartitionModes = {'ECG_contiguous' 'ECG_overlapped' 'QRS'};
        % The fields required in the ECGannotations property
        cAnnotationsFieldNamesRequired = {'time'};
        % The methods required for an ECGtask object
        cObjMethodsRequired = {'Process' 'Start' 'Concatenate'};
        % The properties required for an ECGtask object
        cObjPropsRequired = {'progress_handle' 'name' 'tmp_path'};
        
        % maxQRSxIter: Maximum amount of heartbeats in ECG recordings of 2 leads
        maxQRSxIter = 5e3;
        % Minimum amount of heartbeats to be processed by a PID
        minQRSxIter = 1500;
        % Minimum amount of samples to be processed by a PID
        min_ECG_perPID = 432000; % 10*60*360*2 (10 minutes of 2-lead ECG@360Hz typical MITBIH arrhythmia rec)
        %gzip compression ratio
        typical_compression_ratio = 1.5; %times
        % maximum results file size
        payload_max_size = 2 * 1024^3; % bytes
        % maximum debugging email report file size
        max_report_filesize = 4 * 1024^2; % bytes
        % Maximum time to wait for other PIDs finishing their work.
        Time2WaitPIDs = 0 * 60; % seconds.
        
    end
    
    properties ( Access = private )
        % private
        rec_filename
        % private
        rec_path
        % private
        common_path
        % private
        bHaveUserInterface
        % private
        Loops2io
        % private
        MaxNodesReading
        % private
        MaxNodesWriting
        % private
        QRS_locations
        % private
        bArgChanged = true;
        % private
        bECG_rec_changed = true;
        % private
        bCreated = false;
        % maxECGxIter: Maximum samples ECG
        maxECGxIter
        %list of handles in the toolbox
        cKnownECGtasksHdl
        % private
        cKnownECGtasks
        % private
        ErrorReport = [];
        % using annotations from ECG heartbeat classification task.
        bHB_task_annotations = false;
    end
    
    properties (SetAccess = private, GetAccess = public)  
        % The filename generated as a result of the processing of the
        % ECGtask
        Result_files = []; 
        % Flag check if the task generates a file payload as a result
        doPayload = true;
    end

    properties
        % The ECGtask generated any error ?
        Error
        % Was the ECGtask processed ?
        Processed
        % Had this pid work to do ?
        NoWork2Do
        % temporary path
        tmp_path
        % output path
        output_path
        % the full filename of the ECG recording
        recording_name
        % recording file format
        recording_format
        % annotations performed in the ECG (in MIT format)
        ECG_annotations 
        % Wave delineation performed in the ECG
        ECG_delineation 
        % information about the ECG recording
        ECG_header        
        % total amount of pid's
        cant_pids
        % identification of this pid
        this_pid
        % Force a fixed amount of samples per iter. If <1 it is a fraction
        % of the total amount of samples.
        max_samples_per_iter = [];
        % number of times to repeat the ECGtask
        repetitions
        % handle to the ECGtask object
        ECGtaskHandle
        % partition mode chosen
        partition_mode
        % time of overlapp among signal segments
        overlapping_time
        % is this caching ?
        cacheResults
        % are all pid's finishing the ECGtask together ?
        syncSlavesWithMaster
        % ECGannotations label format for heartbeat types.
        class_labeling = 'AAMI';
        % user string to individualize each run
        user_string
        
    end
    
    methods 
        
        function obj = ECGwrapper(varargin)
        %%
        % object constructor: parse arguments, check if user interface is
        % available
            %% Constants and definitions
           
            %Java user interface is started. Not started in clusters for example.
            obj.bHaveUserInterface = usejava('desktop');

            if( obj.bHaveUserInterface )
                % PC tyle, ignore this.
                obj.Loops2io = 1;
                obj.MaxNodesReading = inf;
                obj.MaxNodesWriting = inf;
            else
                %Cluster settings. Ignore them if running in a PC style computer.
                %Limits of the cluster to have multiple process accesing I/O
                obj.Loops2io = 100;
                obj.MaxNodesReading = 15;
                obj.MaxNodesWriting = 15;
            end

            % discover available ECG tasks installed
            [obj.cKnownECGtasks, obj.cKnownECGtasksHdl ] = list_all_ECGtask();
            
            %% Argument parsing

            %argument definition
            p = inputParser;   % Create instance of inputParser class.
            p.addParamValue('recording_name', [], @(x)( ischar(x) || isempty(x) ));
            aux_val = obj.cKnownFormats;
            p.addParamValue('recording_format', [], @(x)( ischar(x) && any(strcmpi(x,aux_val))) || isempty(x) );
            p.addParamValue('this_pid', 1, @(x)(ischar(x) || (isnumeric(x) && all(x > 0) ) ) );
            p.addParamValue('tmp_path', [], @(x)(ischar(x) || isempty(x) ) );
            p.addParamValue('user_string', [], @(x)(ischar(x) || isempty(x)) );
            p.addParamValue('output_path', [], @(x)(ischar(x) || isempty(x)) );
            p.addParamValue('repetitions', 1, @(x)(isnumeric(x) && x > 0 ) );
            p.addParamValue('ECGtaskHandle', ECGtask_do_nothing, @(x)( isobject(x) || ischar(x) ) );
            p.addParamValue('overlapping_time', 30, @(x)(isnumeric(x) && x > 0 ) );
            aux_val = obj.cPartitionModes;
            p.addParamValue('partition_mode', 'ECG_overlapped', @(x)( ischar(x) && any(strcmpi(x,aux_val))) );
            p.addParamValue('cacheResults', true, @(x)(islogical(x)) );
            p.addParamValue('syncSlavesWithMaster', false, @(x)(islogical(x)) );

            try
                p.parse( varargin{:} );
            catch MyError
                rethrow(MyError);    
            end

            obj.recording_name = p.Results.recording_name;
            obj.recording_format = p.Results.recording_format;
            obj.this_pid = p.Results.this_pid;
            obj.tmp_path = p.Results.tmp_path;
            obj.output_path = p.Results.output_path;
            obj.repetitions = p.Results.repetitions;
            obj.ECGtaskHandle = p.Results.ECGtaskHandle;
            obj.overlapping_time = p.Results.overlapping_time;
            obj.partition_mode = p.Results.partition_mode;
            obj.cacheResults = p.Results.cacheResults;
            obj.syncSlavesWithMaster = p.Results.syncSlavesWithMaster;
            
            % Dont know why this variable uses a lot of bytes to store at disk.
            clear p
            
            obj.bCreated = true;
            obj.Processed = false;
            obj.NoWork2Do = false;
            obj.Error = false;
            
        end
        
        function result_payload = Run(obj)
        % Prepare and execute the ECGtask, creating a payload as a result      

            result_payload = [];
            result_files = {};
            repeat_idx = 1;
            obj.ErrorReport = [];
            obj.Processed = false;
            obj.NoWork2Do = false;
            obj.Error = false;
            

            if( obj.bArgChanged )
                obj = obj.CheckArguments();
                obj.bArgChanged = false;
            end

            fprintf(1, '\n');
            cprintf( 'Blue', disp_option_enumeration( 'Description of the process:', { ['Recording: ' obj.recording_name] ['Task name: ' obj.ECGtaskHandle.name] } ));
            fprintf(1, '\n');
            
            if( isempty(obj.user_string) )
                aux_user_prefix = [];
            else
                if(length(obj.rec_filename) > 30 && length(obj.user_string) > 10 )
                    % shorter version
                    aux_user_prefix = [ adjust_string(obj.user_string, 10, 'center', '_') '_'];
                else
                    aux_user_prefix = [obj.user_string '_'];
                end

            end
            
            % check for cached results
            if( obj.cacheResults )

                if(length(obj.rec_filename) > 30 && length(obj.ECGtaskHandle.name) > 10 )
                    % shorter version
                    task_name = adjust_string(obj.ECGtaskHandle.name, 10, 'center', '_');
                else
                    task_name = obj.ECGtaskHandle.name;
                end
                
                if( isprop(obj.ECGtaskHandle, 'signal_payload') && obj.ECGtaskHandle.signal_payload )
                    % The results is a signal for arbitrary tasks, so check
                    % the result signal instead

                    MIT_filename = [ obj.rec_filename '_' aux_user_prefix task_name ];
                    MIT_filename = regexprep(MIT_filename, '\W*(\w+)\W*', '$1');
                    MIT_filename = regexprep(MIT_filename, '\W', '_');

                    files_this_pid = dir([obj.output_path MIT_filename '*.dat']);
                    
                else
                    
                    files_this_pid = dir([obj.output_path obj.rec_filename '_' aux_user_prefix task_name '*.mat']);

                end

                if( ~isempty(files_this_pid) )
                    obj.Processed = true;
                    obj.NoWork2Do = true;
                    obj.Result_files = strcat(obj.output_path, {files_this_pid(:).name});
                    cellfun(@(a)(cprintf('[1,0.5,0]','Cached results found in %s.\n', a)), obj.Result_files);
                    return
                end
                
            end            
            
            overlapping_samples = obj.overlapping_time * obj.ECG_header.freq;
            
            while ( repeat_idx <= obj.repetitions )

                try

                    %% work starts here

                    %% Prepare jobs to perform.

                    % Activate the progress_struct bar.
                    pb = progress_bar([obj.rec_filename ': ' obj.ECGtaskHandle.name ]);
                    
                    % Initialization of the process.
                    obj.ECGtaskHandle.progress_handle = pb;
                    obj.ECGtaskHandle.tmp_path = obj.tmp_path;
                    obj.ECGtaskHandle.Start( obj.ECG_header, obj.ECG_annotations );
                    
                    if( obj.ECGtaskHandle.started )
                        
                        %% 
                        
                        % find out the ammount of memory in the system
                        if( ispc() )
                            user = memory();
                        else
                            % empirical value from a x64 platform.
                            user.MaxPossibleArrayBytes = 9.3784e+09;
                        end

                        % upper bound of 2 mins @ 12-leads 1000 Hz
                        obj.maxECGxIter = round(obj.ECGtaskHandle.memory_constant*(min(2*60*60*12*1000,user.MaxPossibleArrayBytes/8))); % double = 8 bytes

                        if( strcmpi(obj.partition_mode, 'QRS') )
                            %% QRS mode
                            
                            % this annotations are qrs locations provided
                            % to the task, result of previous QRS
                            % detections / correction, that were parsed in
                            % the Start method.
                            if( isprop(obj.ECGtaskHandle, 'annotations') && ~isempty(obj.ECGtaskHandle.annotations) )
                                % take precedence to QRS locations provided
                                % with the signal.
                                obj.QRS_locations = obj.ECGtaskHandle.annotations;
                                obj.bHB_task_annotations = true;
                            else
                                obj.bHB_task_annotations = false;
                            end

                            if( isprop(obj.ECGtaskHandle, 'min_heartbeats_required') )
                                aux_min_QRS_iter = obj.ECGtaskHandle.min_heartbeats_required;
                            else
                                aux_min_QRS_iter = obj.minQRSxIter;
                            end
                            
                            if( isprop(obj.ECGtaskHandle, 'max_heartbeats_per_iter')  )                            
                                aux_max_QRS_iter = obj.ECGtaskHandle.max_heartbeats_per_iter;
                            else
                                aux_max_QRS_iter = obj.maxQRSxIter;
                            end
                            
                            cant_QRS_locations = length(obj.QRS_locations);

                            %PID parsing
                            if( obj.cant_pids > 1 )

                                max_recommended_cant_pids = max(1, floor( cant_QRS_locations / aux_min_QRS_iter ));

                                if( obj.cant_pids > max_recommended_cant_pids )
                                    warning('ECGwrapper:TooMuchPIDs', 'CantPIDs too large for the work to do, consider decreasing it.');
                                    obj.cant_pids = max_recommended_cant_pids;
                                end

                                [pid_starts, pid_ends] = TaskPartition( cant_QRS_locations, obj.cant_pids);

                                if( obj.this_pid <= obj.cant_pids )
                                    QRS_start_idx = pid_starts(obj.this_pid);
                                    QRS_end_idx = pid_ends(obj.this_pid);
                                else
                                    % Extra PIDs
                                    warning('ECGwrapper:TooMuchPIDs', 'This PID has nothing to do. Exiting.');
                                    return
                                end

                            else
                                %Only one PID
                                QRS_start_idx = 1;
                                QRS_end_idx = cant_QRS_locations;
                                
                                
                            end

                            cant_QRS2do = QRS_end_idx - QRS_start_idx + 1;
                            
                            if( cant_QRS2do > aux_max_QRS_iter )
                                warning('ECGwrapper:TooFewPIDs', 'CantPIDs is too small for the work to do, consider increasing it.');
                            end
                            
                            cant_samples= obj.QRS_locations(QRS_end_idx) - obj.QRS_locations(QRS_start_idx);
                            
                            if( isempty(obj.max_samples_per_iter) )
                                cant_iter = ceil(cant_samples * obj.ECG_header.nsig / obj.maxECGxIter);
                            else
                                if( obj.max_samples_per_iter <= 1 )
                                    cant_iter = round(1/obj.max_samples_per_iter);
                                else
                                    cant_iter = ceil(cant_samples * obj.ECG_header.nsig / obj.max_samples_per_iter);
                                end
                            end
                            %calculate iters.
                            [iter_starts, iter_ends] = TaskPartition( cant_QRS2do, cant_iter);                        

                            % check that none of the iterations exceed the
                            % memory quota
                            cant_samples_each_iter = obj.QRS_locations(iter_ends) - obj.QRS_locations(iter_starts) + 1;
                            iter_exceeded_idx = find( cant_samples_each_iter > (obj.maxECGxIter/obj.ECG_header.nsig) );

                            if( ~isempty(iter_exceeded_idx ) )
                                aux_jj_start = 1;
                                aux_iter_starts = [];
                                aux_iter_ends = [];
                                for  jj= rowvec(iter_exceeded_idx)

                                    aux_iter_starts = [aux_iter_starts ; iter_starts(aux_jj_start:jj-1)];
                                    aux_iter_ends = [aux_iter_ends; iter_ends(aux_jj_start:jj-1)];

                                    aux_jj_start = jj+1;

                                    aux_cant_this_iter = ceil(cant_samples_each_iter(jj) * obj.ECG_header.nsig / obj.maxECGxIter );
                                    [aux_sample_starts, aux_sample_ends] = TaskPartition( cant_samples_each_iter(jj), aux_cant_this_iter);

                                    aux_starts = [ ];
                                    aux_ends = [ ];
                                    for kk = 1:length(aux_sample_starts)
                                        bAux = obj.QRS_locations >= aux_sample_starts(kk) + obj.QRS_locations(iter_starts(jj)) - 1 & obj.QRS_locations <= aux_sample_ends(kk) + obj.QRS_locations(iter_starts(jj)) - 1;
                                        aux_starts = [ aux_starts; find( bAux, 1, 'first' ) ];
                                        aux_ends = [ aux_ends ; find( bAux, 1, 'last' )];
                                    end

                                    aux_iter_starts = [aux_iter_starts ; ( aux_starts )];
                                    aux_iter_ends = [aux_iter_ends; ( aux_ends )];

                                end

                                aux_iter_starts = [aux_iter_starts ; iter_starts(aux_jj_start:cant_iter)];
                                aux_iter_ends = [aux_iter_ends; iter_ends(aux_jj_start:cant_iter)];

                                iter_starts = aux_iter_starts;
                                iter_ends = aux_iter_ends;
                                cant_iter = max(1, length(iter_starts));

                            end   

                        else
                            %% ECG modes

                            %PID parsing
                            if( obj.cant_pids > 1 )

                                max_recommended_cant_pids = max(1, round( obj.ECG_header.nsamp * obj.ECG_header.nsig / obj.min_ECG_perPID ));

                                if( obj.cant_pids > max_recommended_cant_pids )
                                    warning('ECGwrapper:TooMuchPIDs', 'CantPIDs too large for the work to do, consider decreasing it.');
                                    obj.cant_pids = max_recommended_cant_pids;
                                end

                                % make a partition of the ECG
                                [pid_starts, pid_ends] = TaskPartition( obj.ECG_header.nsamp, obj.cant_pids);

                                if( obj.this_pid <= obj.cant_pids )
                                    ECG_start_idx = pid_starts(obj.this_pid);
                                    ECG_end_idx = pid_ends(obj.this_pid);
                                else
                                    % Extra PIDs
                                    warning('ECGwrapper:TooMuchPIDs', 'This PID has nothing to do. Exiting.');
                                    return
                                end

                            else
                                %Only one PID
                                ECG_start_idx = 1;
                                ECG_end_idx = obj.ECG_header.nsamp;
                            end

                            % calculate iterations 
                            cant_samples2do = ECG_end_idx - ECG_start_idx + 1;
                            
                            if( isempty(obj.max_samples_per_iter) )
                                cant_iter = max(1, ceil(cant_samples2do * obj.ECG_header.nsig / obj.maxECGxIter ));
                            else
                                if( obj.max_samples_per_iter <= 1 )
                                    cant_iter = round(1/obj.max_samples_per_iter);
                                else
                                    cant_iter = ceil(cant_samples2do * obj.ECG_header.nsig / obj.max_samples_per_iter);
                                end
                            end
                            
                            [iter_starts, iter_ends] = TaskPartition( cant_samples2do, cant_iter);

                        end

                        if( obj.this_pid > obj.cant_pids )
                            %este pid no trabaja.
                            if(obj.bHaveUserInterface)
                                cprintf('*[1,0.5,0]', 'This PID %d will not work, consider decreasing obj.cant_pids to %d.\n', obj.this_pid, obj.cant_pids);
                            else
                                fprintf(2, 'This PID %d will not work, consider decreasing obj.cant_pids to %d.\n', obj.this_pid, obj.cant_pids);
                            end
                            return
                        end

                        % Update point
                        pb.checkpoint('Initializing');

                        %% Iterations over the whole ECG recording

                        start_iter = 1;
                        payload_files = [];

                        if(length(obj.rec_filename) > 30 && length(obj.ECGtaskHandle.name) > 10 )
                            % shorter version
                            task_name = adjust_string(obj.ECGtaskHandle.name, 10, 'center', '_');
                        else
                            task_name = obj.ECGtaskHandle.name;
                        end
                        
                        %check for previous iterations already done, and try to restore.
                        if( cant_iter > 1 )
                            % look for the previous cached point to continue
                            for ii = 1:cant_iter
                                
                                aux_filename =  [obj.output_path 'tmpfile_' aux_user_prefix task_name '_' obj.rec_filename '_payload_cantpids_' num2str(obj.cant_pids) '_thispid_' num2str(obj.this_pid) '_iteration_' num2str(ii) '_of_' num2str(cant_iter)  '.mat'];
                                
                                if( exist(aux_filename, 'file')  )
                                    payload_files = [ payload_files; cellstr(aux_filename)];
                                    start_iter = ii+1;
                                else
                                    break
                                end
                            end
                        else
                            aux_filename =  [obj.output_path 'tmpfile_' aux_user_prefix task_name '_' obj.rec_filename '_payload_cantpids_' num2str(obj.cant_pids) '_thispid_' num2str(obj.this_pid) '_iteration_1_of_1.mat'];
                            
                            if( exist(aux_filename, 'file')  )
                                payload_files = [ payload_files; cellstr(aux_filename)];
                                start_iter = 2;
                                obj.NoWork2Do = true;
                                cprintf('*[1,0.5,0]', 'All work done. Exiting.\n');
                            end
                        end

                        pb.Loops2Do = cant_iter;

                        for this_iter = start_iter:cant_iter
                            %% loop initialization

                            % por q deberia saberlo ?
                            %obj.ECGtaskHandle.this_iteration = this_iter;

                            %start of the progress_struct loop 0%
                            pb.start_loop();

                            if( strcmpi(obj.partition_mode, 'QRS') )
                                %% la partici�n se hizo en la cantidad de latidos

                                this_iter_QRS_start_idx = QRS_start_idx - 1 + iter_starts(this_iter);
                                this_iter_QRS_end_idx = QRS_start_idx - 1 + iter_ends(this_iter);

                                this_iter_ECG_start_idx = max(1, obj.QRS_locations(this_iter_QRS_start_idx) - overlapping_samples);
                                this_iter_ECG_end_idx = min(obj.ECG_header.nsamp, obj.QRS_locations(this_iter_QRS_end_idx) + overlapping_samples);

                                
                                if( obj.bHB_task_annotations )
                                    
                                    if( ~isempty(obj.QRS_locations) )
                                        this_ann.time = obj.QRS_locations(this_iter_QRS_start_idx:this_iter_QRS_end_idx) - this_iter_ECG_start_idx + 1;
                                    end
                                    
                                else
                                    
                                    % create an annotation struct for this
                                    % iteration.
                                    this_ann = obj.ECG_annotations;

                                    if( ~isempty(this_ann) )

                                        for field_names = rowvec(fieldnames(this_ann))
                                            if( ~isempty( this_ann.(field_names{1}) ) )
                                                this_ann.(field_names{1}) = this_ann.(field_names{1})(this_iter_QRS_start_idx:this_iter_QRS_end_idx);
                                            end
                                        end
                                        this_ann.time = this_ann.time - this_iter_ECG_start_idx + 1;
                                    end
                                    
                                end
                                
                                % in QRS mode, this is not useful.
                                this_iter_ECG_relative_start_end_idx = [1 (this_iter_ECG_end_idx - this_iter_ECG_start_idx + 1)];

                            else
                                %% la partici�n se hizo en el registro de ECG 

                                this_iter_QRS_start_idx = nan;
                                this_iter_QRS_end_idx = nan;
                                this_iter_cant_QRS = nan;

                                this_iter_ECG_start_idx = max(1, ECG_start_idx - 1 + iter_starts(this_iter) - overlapping_samples);
                                this_iter_ECG_end_idx = min(obj.ECG_header.nsamp, ECG_start_idx - 1 + iter_ends(this_iter) + overlapping_samples);

                                % create an annotation struct for this
                                % iteration.
                                this_ann = obj.ECG_annotations;
                                if( ~isempty(this_ann) )
                                    bAux = this_ann.time > this_iter_ECG_start_idx & this_ann.time < this_iter_ECG_end_idx;
                                    for field_names = rowvec(fieldnames(this_ann))
                                        if( ~isempty( this_ann.(field_names{1}) ) )
                                            aux_val = this_ann.(field_names{1});
                                            this_ann.(field_names{1}) = aux_val( bAux );
                                        end
                                    end
                                    this_ann.time = this_ann.time - this_iter_ECG_start_idx + 1;
                                
                                end
                                %Sample where the ECG starts. Useful for
                                %overlapped mode.
                                this_iter_ECG_relative_start_end_idx = [(ECG_start_idx - 1 + iter_starts(this_iter) - this_iter_ECG_start_idx + 1), (ECG_start_idx - 1 + iter_ends(this_iter) - this_iter_ECG_start_idx + 1) ];
                            end

                            % create a header struct for this iteration.
                            this_header = obj.ECG_header;
                            this_header.nsamp = this_iter_ECG_end_idx - this_iter_ECG_start_idx + 1;

                            %% ECG Recording reading
                            % Update point
                            pb.checkpoint('ECG Recording reading');

                            if(strcmpi(obj.ECGtaskHandle.target_units, 'Wrapper') )
                                % pass this same object to the target class.
                                ECG = obj;

                            else
                                %% read from file

                                ECG = read_ECG(obj.recording_name, this_iter_ECG_start_idx, this_iter_ECG_end_idx, obj.recording_format );

                                %convert ADC samples (int16) to obj.ECGtaskHandle.target_units volts.
                                if( ~isempty(obj.ECGtaskHandle.target_units) )

                                    if(strcmpi(obj.ECGtaskHandle.target_units, 'ADCu') )

                                        if( any(any((ECG - fix(ECG)) < -(2^15-1))) || any(any((ECG - fix(ECG))) > 2^15 ) || any(any((ECG - fix(ECG)) ~= 0)) )
                                            % outside int16 range -> match range
                                            % non-integer values
                                            ECG = bsxfun( @minus, ECG, rowvec(median(ECG)) );
                                            aux_val = max(abs(ECG));
                                            ECG = round(bsxfun( @times, ECG, 2^15 ./ rowvec(aux_val) ));
                                        end

                                        ECG = int16(ECG);

                                    else
                                        [ECG this_header] = ADC2units(double(ECG), this_header, obj.ECGtaskHandle.target_units);
                                    end

                                end

                            end

                            %% User defined function calculation

                            % Update point
                            pb.checkpoint([ 'User function @ ' Seconds2HMS(this_iter_ECG_start_idx/this_header.freq) ]);

    %                         if( strcmpi(obj.partition_mode, 'QRS') )
    %                             fprintf(1, 'ECG from: %d - %d\n', obj.QRS_locations(this_iter_QRS_start_idx) , obj.QRS_locations(this_iter_QRS_end_idx) );
    %                         else
    %                             fprintf(1, 'ECG from: %d - %d\n', this_iter_ECG_start_idx + this_iter_ECG_relative_start_end_idx - 1);
    %                         end

                            %user-defined execution according to 'process'
                            payload = obj.ECGtaskHandle.Process(ECG, this_iter_ECG_start_idx, this_iter_ECG_relative_start_end_idx, this_header, this_ann, [this_iter_QRS_start_idx this_iter_QRS_end_idx] );

                            if( ~isempty(payload) )

                                if( obj.ECGtaskHandle.doPayload )

                                    % Update point
                                    pb.checkpoint('Saving results');

                                    if(length(obj.rec_filename) > 30 && length(obj.ECGtaskHandle.name) > 10 )
                                        % shorter version
                                        task_name = adjust_string(obj.ECGtaskHandle.name, 10, 'center', '_');
                                    else
                                        task_name = obj.ECGtaskHandle.name;
                                    end
                                    
                                    % save results
                                    save([obj.tmp_path 'tmpfile_' aux_user_prefix task_name '_' obj.rec_filename '_payload_cantpids_' num2str(obj.cant_pids) '_' num2str(obj.this_pid) '_' num2str(this_iter) '_' num2str(cant_iter) '.mat'], '-struct', 'payload');

                                    % rename results.
                                    auxStr = [obj.tmp_path 'tmpfile_' aux_user_prefix task_name '_' obj.rec_filename '_payload_cantpids_' num2str(obj.cant_pids) '_thispid_' num2str(obj.this_pid) '_iteration_' num2str(this_iter) '_of_' num2str(cant_iter)  '.mat'];

                                    movefile(   [obj.tmp_path 'tmpfile_' aux_user_prefix task_name '_' obj.rec_filename '_payload_cantpids_' num2str(obj.cant_pids) '_' num2str(obj.this_pid) '_' num2str(this_iter) '_' num2str(cant_iter)  '.mat'], ...
                                                auxStr, 'f' );

                                    payload_files = [ payload_files; cellstr(auxStr)];

                                end

                            end

                            pb.end_loop();

                        end


                        % move to the destination folder if needed
                        if( ~strcmpi(obj.tmp_path, obj.output_path) )

                            pause_time = 1;
                            for jj = 1:length(payload_files)

                                % sometimes this file move takes long
                                % time, specially in distributed
                                % filesystems.

                                tic_id = tic;
                                
                                if( ~strcmpi(fileparts(payload_files{jj}), fileparts(obj.output_path)) )
                                    movefile( payload_files{jj}, obj.output_path, 'f' );
                                end

                                time_elapsed = toc(tic_id);

                                if( time_elapsed > pause_time )
                                    % network congestion, wait for a while
                                    pause(pause_time)
                                    pause_time = min(180, 3 * pause_time );
                                else
                                    pause_time = max(1, pause_time / 1.5 );
                                end

                            end                            
                        end    

                        %indicate end of loop.
                        pb.reset();

                        %clear some not needed variables 
                        clear *ECG this_iter* ECG_annotations

                        if( obj.ECGtaskHandle.doPayload )
                            % if PIDs generates payloads, try to collect them.
                            %% multiPIDs. Try to build the whole thing from every parts, otherwise exit
                            if( obj.this_pid == obj.cant_pids )
                            %% Master PID
                            %last pid

                                % Update point
                                pb.checkpoint('Collecting results');

                                %%Master PID last pid
                                bContinue = true;

                                %Wait for obj.Time2WaitPIDs seconds the finalization of all PIDs. Otherwise exit
                                %with error.
                                Start2Wait = tic();
                                start_pid = 1;
                                start_iter_this_pid = 1;
                                resume_iter_this_pid = 0;
                                payload_dump_counter = 1;

                                while(bContinue)

                                    try

                                        % Update point
                                        pb.checkpoint('Collecting other PIDs results');

                                        this_header = [];
                                        
                                        %feature matrices
                                        payload = [];
                                        payload_idx = 1;

                                        % to resume after the last dump to disk.
                                        if( resume_iter_this_pid > 0 )
                                            start_iter_this_pid = resume_iter_this_pid;
                                        end

                                        if(length(obj.rec_filename) > 30 && length(obj.ECGtaskHandle.name) > 10 )
                                            % shorter version
                                            task_name = adjust_string(obj.ECGtaskHandle.name, 10, 'center', '_');
                                        else
                                            task_name = obj.ECGtaskHandle.name;
                                        end
                                        
                                        for ii = start_pid:obj.cant_pids

                                            files_this_pid = dir([obj.output_path 'tmpfile_' aux_user_prefix task_name '_' obj.rec_filename '_payload_cantpids_' num2str(obj.cant_pids) '_thispid_' num2str(ii) '_*.mat' ]);

                                            if( isempty(files_this_pid) )

                                                if(obj.cant_pids == 1)
                                                    % no payload generated
                                                    bContinue = false;
                                                    break
                                                else
%                                                     fprintf(2,'%s not found\n', [obj.output_path 'tmpfile_' aux_user_prefix obj.ECGtaskHandle.name '_' obj.rec_filename '_payload_cantpids_' num2str(obj.cant_pids) '_thispid_' num2str(ii) '_*.mat' ]);
                                                    error('ECGwrapper:PIDnotFinished', 'Handled error');
                                                end
                                            end

                                            cant_iteraciones_this_pid = length(files_this_pid);

                                            % por q deberia saberlo ?
                                            %obj.ECGtaskHandle.cant_iteration = cant_iteraciones_this_pid;
                                            
                                            if(length(obj.rec_filename) > 30 && length(obj.ECGtaskHandle.name) > 10 )
                                                % shorter version
                                                task_name = adjust_string(obj.ECGtaskHandle.name, 10, 'center', '_');
                                            else
                                                task_name = obj.ECGtaskHandle.name;
                                            end
                                            
                                            for jj = start_iter_this_pid:cant_iteraciones_this_pid 

                                                % por q deberia saberlo ?
                                                %obj.ECGtaskHandle.this_iteration = jj;

                                                aux_FM_filename = [obj.output_path 'tmpfile_' aux_user_prefix task_name '_' obj.rec_filename '_payload_cantpids_' num2str(obj.cant_pids) '_thispid_' num2str(ii) '_iteration_' num2str(jj) '_of_' num2str(cant_iteraciones_this_pid) '.mat' ];
                                                if( ~exist(aux_FM_filename, 'file'))
                                                    %Probably this PID not finished yet.
                                                    error('ECGwrapper:PIDnotFinished', 'Handled error');
                                                end
                                                aux = load(aux_FM_filename);

                                                if( isprop(obj.ECGtaskHandle, 'signal_payload') && obj.ECGtaskHandle.signal_payload  )
                                                % Process the results as a
                                                % signal, dump sequentialy
                                                % to disk

                                                    % gain offset calculation to fit in destination class int16
                                                    % offset/gain transform
%                                                         range_conversion_offset = mean(obj.ECGtaskHandle.range_min_max_tracking, 2);
%                                                         range_conversion_gain = (2^15-1)./(diff(obj.ECGtaskHandle.range_min_max_tracking,1, 2)./2); 

                                                    % gain transform
                                                    range_conversion_offset = 0;
                                                    range_conversion_gain = (2^15-1)./(max(max(abs(obj.ECGtaskHandle.range_min_max_tracking)))); 

                                                    result_signal = cast( round( bsxfun( @times, bsxfun( @minus, aux.result_signal, range_conversion_offset), range_conversion_gain) ) , 'int16');
                                                
                                                    if( isempty(this_header) )

                                                        
                                                        this_header = obj.ECG_header;

                                                        % in ADCu / physical units
                                                        this_header.gain = repmat(range_conversion_gain, this_header.nsig, 1);
                                                        this_header.offset = repmat(range_conversion_offset, this_header.nsig, 1) ;
                                                        [this_header.nsamp, this_header.nsig] = size(result_signal);

                                                        if( obj.repetitions > 1)
                                                            aux_sufix = [ '_repetition_' num2str(repeat_idx)];
                                                        else
                                                            aux_sufix = [];                                                       
                                                        end

                                                        MIT_filename = [ obj.rec_filename '_' aux_user_prefix task_name aux_sufix ];
                                                        MIT_filename = regexprep(MIT_filename, '\W*(\w+)\W*', '$1');
                                                        MIT_filename = regexprep(MIT_filename, '\W', '_');

                                                        this_header.recname = MIT_filename;
                                                        
                                                        this_header.freq = round(obj.ECGtaskHandle.sampling_rate_out);
                                                        
                                                        result_files = [ result_files; cellstr([obj.output_path MIT_filename '.dat'])];
                                                        
                                                        fidECG = fopen([obj.output_path MIT_filename '.dat'], 'w');
                                                        
                                                    else
                                                        this_header.nsamp = this_header.nsamp  + size(result_signal, 1);
                                                        
                                                        fidECG = fopen([obj.output_path MIT_filename '.dat'], 'a');
                                                    end
                                                    
                                                    
                                                    try
                                                        fwrite(fidECG, result_signal', 'int16', 0 );
                                                        fclose(fidECG);
                                                    catch MEE
                                                        fclose(fidECG);
                                                        rethrow(MEE);
                                                    end

                                                else
                                                % Process the results according to the
                                                % 'Concatenate' method
                                                    payload = obj.ECGtaskHandle.Concatenate(payload, aux);

                                                    %check variable growth
                                                    payload_var_size = whos('payload');
                                                    payload_var_size = payload_var_size.bytes/obj.typical_compression_ratio;

                                                    if( payload_var_size > obj.payload_max_size )
                                                        %force dump to disk.
                                                        payload.payload_idx = payload_idx;
                                                        save([obj.output_path 'payloads_' aux_user_prefix obj.ECGtaskHandle.name '_' obj.rec_filename '_' num2str(payload_dump_counter) '.mat'], '-struct', 'payload');
                                                        %update counter and reset payload.
                                                        payload_idx = payload_idx + size(payload,1) + 1;
                                                        payload = [];
                                                        payload_dump_counter = payload_dump_counter + 1;
                                                        %to continue from this pid/iteration.
                                                        start_pid = ii;
                                                        resume_iter_this_pid = jj+1;
                                                    end
                                                end

                                            end

                                            %the first time in this loop we resume after the last
                                            %file dump, and never again.
                                            start_iter_this_pid = 1;

                                        end

                                        % write finally the header info.
                                        if( isprop(obj.ECGtaskHandle, 'signal_payload') && obj.ECGtaskHandle.signal_payload  )
                                            writeheader(obj.output_path, this_header);   
                                        end
                                        
                                        % Update point
                                        pb.checkpoint('Generating final results');

                                        if( ~isempty(payload) )
                                            % data remain in memory -> dump 2 disk
                                            if( payload_idx > 1 )
                                                payload.payload_idx = payload_idx;
                                            end
                                            save([obj.output_path 'payloads_' aux_user_prefix obj.ECGtaskHandle.name '_' obj.rec_filename '_' num2str(payload_dump_counter) '.mat'], '-struct', 'payload');
                                            %update counter and reset payload.
                                            payload = [];
                                        end

                                        %rename results to put some ordering.
                                        files_this_pid = dir([obj.output_path 'payloads_' aux_user_prefix obj.ECGtaskHandle.name '_' obj.rec_filename '_*.mat']);
                                        cant_iteraciones_this_pid = length(files_this_pid);

                                        if( isempty(obj.user_string) )
                                            aux_user_prefix = [];
                                        else
                                            if(length(obj.rec_filename) > 30 && length(obj.user_string) > 10 )
                                                % shorter version
                                                aux_user_prefix = [ adjust_string(obj.user_string, 10, 'center', '_') '_'];
                                            else
                                                aux_user_prefix = [obj.user_string '_'];
                                            end
                                        end

                                        if(length(obj.rec_filename) > 30 && length(obj.ECGtaskHandle.name) > 10 )
                                            % shorter version
                                            task_name = adjust_string(obj.ECGtaskHandle.name, 10, 'center', '_');
                                        else
                                            task_name = obj.ECGtaskHandle.name;
                                        end

                                        for jj = 1:cant_iteraciones_this_pid

                                            if( cant_iteraciones_this_pid > 1 )
                                                aux_sufix = ['_part_' num2str(jj) '_of_' num2str(cant_iteraciones_this_pid)];
                                            else
                                                aux_sufix = [];
                                            end

                                            if( obj.repetitions > 1)
                                                aux_sufix = [aux_sufix '_repetition_' num2str(repeat_idx)];
                                            end

                                            auxStr = [obj.output_path obj.rec_filename '_' aux_user_prefix task_name aux_sufix '.mat'];

                                            movefile(   [obj.output_path 'payloads_' aux_user_prefix obj.ECGtaskHandle.name '_' obj.rec_filename '_' num2str(jj) '.mat'], ...
                                                        auxStr, 'f' );

                                            result_files = [ result_files; cellstr(auxStr)];

                                        end

                                        % Update point
                                        pb.checkpoint('Deleting temporal files');

                                        %clean temporal files
                                        delete([obj.output_path 'tmpfile_' aux_user_prefix task_name '_' obj.rec_filename '_payload_cantpids_' num2str(obj.cant_pids) '*.mat' ]);

                                        bContinue = false;


                                    catch ME
                                        % Error handling

                                        if( strfind(ME.identifier, 'ECGwrapper') )
                                            if( obj.cant_pids == 1 || obj.bHaveUserInterface || toc(Start2Wait) > obj.Time2WaitPIDs )
%                                                 fprintf(2, 'Timeout. Master give up waitng Slaves.\n');
                                                error('ECGwrapper:PIDnotFinished', 'Master give up waiting for PID %d iteration %d\n', ii, jj);
                                            end
                                            pause(30);
                                        else
                                            rethrow(ME)
                                        end
                                    end

                                end

                                if( isprop(obj.ECGtaskHandle, 'signal_payload') )
                                    
                                    if( ~obj.ECGtaskHandle.signal_payload )
                                        %Dump results as standard payload.
                                        for result_fn = rowvec(result_files)
                                            % Master PID can operate over the global
                                            % payload.
                                            payload = load(result_fn{1});
                                            payload = obj.ECGtaskHandle.Finish( payload, obj.ECG_header );
                                            save(result_fn{1}, '-struct', 'payload');
                                        end
                                    end
                                    
                                else
                                    
                                    %Dump results as standard payload.
                                    for result_fn = rowvec(result_files)
                                        % Master PID can operate over the global
                                        % payload.
                                        payload = load(result_fn{1});
                                        payload = obj.ECGtaskHandle.Finish( payload, obj.ECG_header );
                                        save(result_fn{1}, '-struct', 'payload');
                                    end

                                end
                            else
                            %% Slave PIDs

                                if( obj.syncSlavesWithMaster )
                                % other PIDS sync here after Master build the
                                % results file/s.

                                    if(length(obj.rec_filename) > 30 && length(obj.ECGtaskHandle.name) > 10 )
                                        % shorter version
                                        task_name = adjust_string(obj.ECGtaskHandle.name, 10, 'center', '_');
                                    else
                                        task_name = obj.ECGtaskHandle.name;
                                    end
                                
                                    bContinue = true;
                                    %Wait for Time2WaitPIDs seconds the finalization of all PIDs. Otherwise exit
                                    %with error.
                                    Start2Wait = tic();

                                    while(bContinue)

                                        try

                                            files_this_pid = dir([obj.tmp_path 'tmpfile_' aux_user_prefix task_name '_' obj.rec_filename '_payload_cantpids_' num2str(obj.cant_pids) '_thispid_' num2str(obj.this_pid) '_*.mat' ]);

                                            if( ~isempty(files_this_pid) )
                                                error('ECGwrapper:MasterNotFinished', 'Handled error');
                                            end


                                            bContinue = false;

                                        catch ME

                                            if( strfind(ME.identifier, 'ECGwrapper') )
                                                if( obj.bHaveUserInterface || toc(Start2Wait) > 1.5*obj.Time2WaitPIDs )
                                                    error('ECGwrapper:MasterNotFinished', 'Timeout. Slave give up waitng Master.');
                                                end
                                                pause(30);
                                            else
                                                rethrow(ME)
                                            end
                                        end


                                    end                        

                                end

                            end

                        else
                            obj.ECGtaskHandle.Finish( payload, obj.ECG_header );
                        end

                    else
                        cprintf('[1,0.5,0]', 'Requirements not satisfied in %s for task %s.\n', obj.recording_name, obj.ECGtaskHandle.name);
                    end
                    
                catch MException
                    %% Error handling

                    obj.Error = true;
                    
                    if( obj.bHaveUserInterface )
                        %% with UI
                        if( ~isempty(strfind(MException.identifier, 'ECGwrapper')) || isempty(MException.identifier) )
                            %Our errors
                            if( strfind(MException.identifier, 'ArgCheck') )
                                %Argument error, try other settings if user interface is
                                %available
                                fprintf(2, '%s', MException.message);
                                fprintf(2, '\nTry a different set of arguments.\n');
                            else
                                if( isempty(MException.identifier) )
                                    %User break with CTRL-C ??
                                    fprintf(2, '\nUser interruption.\n');
                                else
                                    obj.ErrorReport = [obj.ErrorReport {getReport(MException)}];
                                    %other home-made errors. Make an educated exit ...
                                    rethrow(MException)
                                end
                            end

                        else
                            
                            obj.ErrorReport = getReport(MException);
                            
                            %% Other unknown errors
                            rethrow(MException);
                            
                        end
                        
                    else
                        %% No User interface report clearly to log file
                        
                        str_aux = disp_string_framed(0, [ 'Error in ' obj.recording_name ] );
                        
                        report = getReport(MException);
                        
                        fprintf(2, '\n\n%s\n%s\n%s', str_aux, report, str_aux);

                        rethrow(MException)

                    end


                end

                repeat_idx = repeat_idx + 1;
                
                % Update point
                pb.checkpoint('Work done.');

            end

            % if execution did not report errors, mark as processed.
            obj.Processed = true;
            obj.Result_files = result_files;
            
            if( isempty(result_files) )
                
                if( obj.ECGtaskHandle.started )
                    if( obj.cant_pids == 1 || obj.this_pid == obj.cant_pids )
                        obj.Error = true;
                        obj.ErrorReport = [obj.ErrorReport {sprintf('No payload generated for recording %s\n', obj.recording_name )}];
                        disp_string_framed(2, 'No payload generated');
                    else
                        disp_string_framed('*Blue', 'Work done!');
                    end
                else
                    obj.NoWork2Do = true;
                    disp_string_framed('[1,0.5,0]', 'Nothing to do here');
                end
                
            else
                
                disp_string_framed('*Blue', 'Work done!');

                if( nargout > 0 )
                    % result required

                    if( cant_iteraciones_this_pid > 1 )
                        result_payload = obj.Result_files;
                    elseif( cant_iteraciones_this_pid == 1 )
                        result_payload = load(obj.Result_files{1});
                    else
                        clear result_payload
                    end
                else
                    clear result_payload
                    fprintf(1, '\n');
                    cprintf( 'Blue', disp_option_enumeration( 'Results saved in', obj.Result_files));
                    fprintf(1, '\n');
                end
            end
            
            % destroy the progress bar.
            pb.delete;
            
        end
        
        function ECG = read_signal(obj, ECG_start_idx, ECG_end_idx)
        % in case using the this object just as an I/O interface, this
        % method can read the samples of a recording.
        
            if( obj.bArgChanged )
                obj = obj.CheckArguments();
                obj.bArgChanged = false;
            end
            
            if( nargin < 2 )
                ECG_start_idx = 1;
            end
            
            if( nargin < 3 )
                ECG_end_idx = ECG_start_idx + 10 * obj.ECG_header.freq;
            end
            
            ECG = read_ECG(obj.recording_name, max(1,ECG_start_idx), min(obj.ECG_header.nsamp,ECG_end_idx), obj.recording_format );
            
        end
        
        function result_files = GetCahchedFileName(obj, task_names)
        % this method is useful to check if there are cached work already
        % created from an ECGtask
            
            if( nargin < 2 || isempty(task_names) )
                task_names = {obj.ECGtaskHandle};
            end

            if( ~iscell(task_names) && ~ischar(task_names))
                task_names = {obj.ECGtaskHandle};
            elseif( ischar(task_names) )
                task_names = {task_names};
            end
            
            result_files = {};
            
            prev_ECGtaskHandle = obj.ECGtaskHandle;
            
            for this_name = rowvec(task_names)
                
                obj.ECGtaskHandle = this_name{1};

                if( obj.bArgChanged )
                    obj = obj.CheckArguments();
                    obj.bArgChanged = false;
                end
                
                if( isempty(obj.user_string) )
                    aux_user_prefix = [];
                else
                    if(length(obj.rec_filename) > 30 && length(obj.user_string) > 10 )
                        % shorter version
                        aux_user_prefix = [ adjust_string(obj.user_string, 10, 'center', '_') '_'];
                    else
                        aux_user_prefix = [obj.user_string '_'];
                    end
                end

                % shorter version
                task_name = adjust_string(obj.ECGtaskHandle.name, 10, 'center', '_');

                files_this_pid = dir([obj.output_path obj.rec_filename '_' aux_user_prefix task_name '.mat']);
                
                task_name = obj.ECGtaskHandle.name;
                
                files_this_pid = [files_this_pid; dir([obj.output_path obj.rec_filename '_' aux_user_prefix task_name '.mat'])];

                if( ~isempty(files_this_pid) )
                    result_files = [ result_files; strcat(obj.output_path, colvec({files_this_pid(:).name})) ];
                end

            end

            % leave things as we found.
            obj.ECGtaskHandle = prev_ECGtaskHandle;
            
            obj.CheckTaskHandle();
            
        end
   
        function ReportErrors(obj)
        % error reporting method.
        
            if( obj.Processed )
                if( obj.Error )
                    
                    fprintf(2, 'Some error happened in %s\n', obj.ECGtaskHandle.name );
                    
                    cant_errors = length(obj.ErrorReport);
                    for ii = 1:cant_errors
                        if( cant_errors > 1)
                            disp_string_framed(2, sprintf('Error %d', ii));
                        end
                        fprintf(2, '\n%s\n', obj.ErrorReport{ii} );
                    end
                    
                end
            else
                fprintf(1, 'Task not processed yet. Execute ''Run'' method first\n');
            end            
        end
        
        function ECGw_out = horzcat(obj, varargin)
        % this method allows signal concatenation .

            ECGw_out = obj;
            
            this_header = ECGw_out.ECG_header;
            
            for ii = 1:length(varargin)
            
                this_header = concat_name( this_header, varargin{ii}, false );
                
            end
            
            if( isempty(ECGw_out.output_path) )
                ECGw_out.output_path = ECGw_out.rec_path;
            end

            if(ECGw_out.output_path(end) ~= filesep )
                ECGw_out.output_path = [ECGw_out.output_path filesep];
            end
            
            cached_filename = [obj.output_path this_header.recname '.dat'];
            
            if( obj.cacheResults && exist( cached_filename, 'file') )
                
                cprintf('[1,0.5,0]','Cached results found in %s.\n', cached_filename);
                
                ECGw_out = ECGwrapper( 'recording_name', [obj.output_path this_header.recname '.hea'] );
                
            else
                
                for ii = 1:length(varargin)

                    ECGw_out = horzcat_internal(ECGw_out, varargin{ii}, this_header.recname);

                end
            end
        
        end
            
        function ECGw_out = vertcat(obj, varargin)
        % this method allows signal concatenation .
            
            ECGw_out = obj;

            this_header = ECGw_out.ECG_header;
            
            for ii = 1:length(varargin)
            
                this_header = concat_name( this_header, varargin{ii}, true );
                
            end
            
            if( isempty(ECGw_out.output_path) )
                ECGw_out.output_path = ECGw_out.rec_path;
            end

            if(ECGw_out.output_path(end) ~= filesep )
                ECGw_out.output_path = [ECGw_out.output_path filesep];
            end
            
            cached_filename = [obj.output_path this_header.recname '.dat'];
            
            if( obj.cacheResults && exist( cached_filename, 'file') )
                
                cprintf('[1,0.5,0]','Cached results found in %s.\n', cached_filename);
                
                ECGw_out = ECGwrapper( 'recording_name', [obj.output_path this_header.recname '.hea'] );
                
            else
                
                for ii = 1:length(varargin)

                    ECGw_out = vertcat_internal(ECGw_out, varargin{ii});

                end
            end
        
        end
            
        function disp(obj)
        % this method produces a pretty-printed description about the
        % information stored in the object
            
            for ii = 1:length(obj)
                
                this_obj = obj(ii);
                
                if( this_obj.bCreated )

                    disp_string_framed( '*Blue', 'ECGwrapper object config' );

                    fprintf(1, '+ECG recording: ');
                    
                    if( isempty(this_obj.recording_name) )
                        cprintf('*Red', 'None selected\n');                    
                    else
                        cprintf('*Blue', '%s (%s)\n', this_obj.recording_name, this_obj.recording_format);                    
                    end
                    
                    fprintf(1,[ ...
                                '+PID: %d/%d\n' ... 
                                '+Repetitions: %d\n' ...
                                '+Partition mode: %s\n'], ... 
                                this_obj.this_pid, ...
                                this_obj.cant_pids, ...
                                this_obj.repetitions, ...
                                this_obj.partition_mode);

                    fprintf(1, '+Function name: ');
                    if( isobject(this_obj.ECGtaskHandle) )
                        
                        if( strcmpi( this_obj.ECGtaskHandle.name, 'Null task' ) )
                            cprintf('*Red', '%s\n', this_obj.ECGtaskHandle.name);                    
                        else
                            fprintf(1, '%s\n', this_obj.ECGtaskHandle.name);
                        end
                    else
                        cprintf('*Red', 'Task not defined\n');                    
                    end
                            
                    fprintf(1, '+Processed: ');                    
                    if( this_obj.Processed)
                        cprintf('*Magenta', 'true\n');                    
                    else
                        cprintf('*Red', 'false\n');                    
                    end

                    if( ~isempty(this_obj.tmp_path) )
                        fprintf(1, '+TMP: %s\n', adjust_string(this_obj.tmp_path, 20) );
                    end

                    if( this_obj.Processed )
                        fprintf(1, disp_option_enumeration( '+Result files:', this_obj.Result_files));
                    end

                    fprintf(1, '\n');
                    
                else
                    fprintf(1, 'Object not created yet.');
                end
                
            end
        end

        %% property access methods.
        
        function set.tmp_path(obj,value)
            
            if( isempty(value) )
                obj.tmp_path = value;
            elseif( (ischar(value) && exist(value, 'dir') ) )
                if( value(end) == filesep )
                    obj.tmp_path = value;
                else
                    obj.tmp_path = [value filesep];
                end
                obj.bArgChanged = true;
            else
                warning('ECGwrapper:BadArg', 'tmp_path must be a string.');
            end
            
        end
        
        function set.user_string(obj,value)
            if( isempty(value) )
                obj.user_string = value;
                
            elseif(ischar(value) ) 

                obj.user_string = value;
                
            else
                warning('ECGwrapper:BadArg', 'user_string must be a string.');
            end
        end
        
        function set.output_path(obj,value)
            if( isempty(value) )
                obj.output_path = value;
                
            elseif(ischar(value) ) 
                
                if( ~exist(value, 'dir') )
                    if( ~mkdir(value) )
                        warning('ECGwrapper:BadArg', 'output_path must exist, or privileges should be granted to this script.');
                        return
                    end
                end
                
                if( value(end) == filesep )
                    obj.output_path = value;
                else
                    obj.output_path = [value filesep];
                end
                obj.bArgChanged = true;
                
            else
                warning('ECGwrapper:BadArg', 'output_path must be a string.');
            end
        end
        
        function set.recording_name(obj,value)
            if( ischar(value) || isempty(value) )
                obj.recording_name = value;
                obj.bArgChanged = true;
                obj.bECG_rec_changed = true;
            else
                warning('ECGwrapper:BadArg', 'recording_name must be a string.');
            end
        end

        function set.recording_format(obj,value)
            if( isempty(value) )
                obj.recording_format = 'auto';
                obj.bArgChanged = true;
            elseif( ischar(value) )
                obj.recording_format = value;
                obj.bArgChanged = true;
                obj.bECG_rec_changed = true;
            else
                warning('ECGwrapper:BadArg', 'recording_format must be a string.');
            end
        end

        function set.ECG_annotations(obj,value)
            
            if( isempty(value) )
                obj.ECG_annotations = [];
                obj.QRS_locations = [];                
            else
                if( isstruct(value) )

                    if( all(isfield(value, obj.cAnnotationsFieldNamesRequired )) )

                        obj.ECG_annotations = value;
                        obj.QRS_locations = value.time;

                    else
                        warning( 'ECGwrapper:BadArg', disp_option_enumeration( 'Please provide the following fields in the annotations struct:',  obj.cAnnotationsFieldNamesRequired) );
                    end

                else
                    warning('ECGwrapper:BadArg', 'ECG_annotations must be a struct.');
                end
            end
        end

        function value = get.ECG_annotations(obj)
            
            if( obj.bECG_rec_changed )
                obj = obj.CheckECGrecording();
            end
            
            value = obj.ECG_annotations;
        end
        
        function value = get.ECG_delineation(obj)
            
            if( obj.bECG_rec_changed )
                obj = obj.CheckECGrecording();
            end
            
            value = obj.ECG_delineation;
        end
                
        function value = get.ECG_header(obj)
            
            if( obj.bECG_rec_changed )
                obj = obj.CheckECGrecording();
            end
            
            value = obj.ECG_header;
            
        end
        
        function set.cant_pids(obj,value)
            if( value > 0 )
                obj.cant_pids = value;
                obj.bArgChanged = true;
            else
                warning('ECGwrapper:BadArg', 'cant_pids must be > 0.');
            end
        end

        function set.this_pid(obj,value)

            [tp, cp] = parse_pids( value );

            if( all(~isnan([tp, cp])) && tp > 0 && tp <= cp)
                obj.this_pid = tp;
                obj.cant_pids = cp; %#ok<MCSUP>
                obj.bArgChanged = true; %#ok<MCSUP>
            else
                warning('ECGwrapper:BadArg', 'Incorrect format, use ''1/3'' or [1 3]');
            end
        end

        function set.max_samples_per_iter(obj,value)

            if( value > 0 )
                obj.max_samples_per_iter = value;
                obj.bArgChanged = true;
            else
                warning('ECGwrapper:BadArg', 'max_samples_per_iter must be > 0 and <= 1.');
            end
            
        end

        function set.repetitions(obj,value)
            if( value > 0 )
                obj.repetitions = value;
                obj.bArgChanged = true;
            else
                warning('ECGwrapper:BadArg', 'repetitions must be > 0.');
            end
        end

        function set.ECGtaskHandle(obj,value)
            if( isa(value, 'ECGtask') )
                obj.ECGtaskHandle = value;
                obj.CheckTaskHandle();
                
            elseif( ischar(value) )
                
                aux_idx = find(strcmpi(obj.cKnownECGtasks, value));
                if( isempty(aux_idx) )
                    cprintf( 'Red', disp_option_enumeration( 'Invalid ECGtask name. ECGtaskHandle must be one of these strings:', obj.cKnownECGtasks ));
                    obj.ECGtaskHandle = ECGtask_do_nothing();
                else
                    obj.ECGtaskHandle = obj.cKnownECGtasksHdl{aux_idx};
                    obj.CheckTaskHandle();
                end
            elseif( isempty(value) )
                
                obj.ECGtaskHandle = ECGtask_do_nothing();
                
            else
                warning('ECGwrapper:BadArg', 'ECGtaskHandle must be a ECGtask object.');
            end
        end
        
        function set.cacheResults(obj,value)
            if( islogical(value) )
                obj.cacheResults = value;
            else
                warning('ECGwrapper:BadArg', 'cacheResults must be boolean.');
            end
        end
        
        function set.syncSlavesWithMaster(obj,value)
            if( islogical(value) )
                obj.syncSlavesWithMaster = value;
            else
                warning('ECGwrapper:BadArg', 'syncSlavesWithMaster must be boolean.');
            end
        end
              
    end
    
    methods (Access = private)  
        
        function obj = CheckTaskHandle(obj)
        %%
        % 
        
            %Object parsing
            if( isobject(obj.ECGtaskHandle) )
                for ii = 1:length(obj.cObjMethodsRequired)
                    if( ~ismethod(obj.ECGtaskHandle, obj.cObjMethodsRequired{ii}) )
                       error( 'ECGwrapper:ArgCheck:UserObjHdl', ['Method ' obj.cObjMethodsRequired{ii} ' not implemented in UserObjHdl.\n\n'] );
                    end        
                end

                for ii = 1:length(obj.cObjPropsRequired)
                    if( ~isprop(obj.ECGtaskHandle, obj.cObjPropsRequired{ii}) )
                       error( 'ECGwrapper:ArgCheck:UserObjHdl', ['Property ' obj.cObjPropsRequired{ii} ' not present in UserObjHdl.\n\n'] );
                    end        
                end
            else
               error( 'ECGwrapper:ArgCheck:UserObjHdl', 'ECGtaskHandle is not a valid ECGtask handle.' );
            end    
        end
        
        
        function obj = CheckArguments(obj)

            %Object parsing
            obj.CheckTaskHandle();

            if( isprop(obj.ECGtaskHandle, 'min_heartbeats_required') || isprop(obj.ECGtaskHandle, 'max_heartbeats_per_iter')  )
                obj.partition_mode = 'QRS';
            elseif( isprop(obj.ECGtaskHandle, 'min_length_required') || isprop(obj.ECGtaskHandle, 'max_length_per_iter') )
                obj.partition_mode = 'ECG_overlapped';
            end

            %ECG parsing
            if( obj.bECG_rec_changed )
                obj = obj.CheckECGrecording();
            end
            
            if( isempty(obj.tmp_path) )
                
                obj.tmp_path = tempdir ;
                
            %         obj.tmp_path = [fileparts(mfilename('fullpath')) filesep 'tmp' filesep ];
            %     obj.tmp_path = [fileparts(obj.recording_name) filesep ];
            end
            
            if( isempty(obj.output_path) )
                obj.output_path = obj.rec_path;
            end

            %check path integrity.
            if(~exist(obj.tmp_path, 'dir'))
                %try to create it
                if( ~mkdir(obj.tmp_path) )
                    error('ECGwrapper:ArgCheck:InvalidPath', 'Invalid obj.tmp_path. Please provide a valid path.\n' );
                end
            end

            if( strcmpi(obj.partition_mode, 'QRS') )
                %partition using QRS detections
                if( obj.overlapping_time < 5 ) %seconds
                    warning('ECGwrapper:ArgCheck:Overlapp_too_low', 'The overlapping time between iterations is too low, consider increasing.\n' );
                end
            else
                if( strcmpi(obj.partition_mode, 'ECG_contiguous') )
                    %One segment after the other.
                    obj.overlapping_time = 0;            
                end
            end

        end

        
        function obj = CheckECGrecording(obj)
        % internal method to check the correctness of the user input.
        
            single_lead_positions = [];
            
            if( isempty(obj.recording_name) )
                error( 'ECGwrapper:ArgCheck:InvalidECGarg', 'Please provide an ECG recording as described in the documentation, help(''ECGwrapper'') maybe could help you.\n' );
            else
                
                % ECG to be read
                [obj.rec_path, obj.rec_filename, caca] = fileparts(obj.recording_name);
                
                if( isempty(obj.rec_path) )
                    obj.rec_path = ['.' filesep];
                end
                
                aux_filename = fullfile(obj.rec_path, obj.rec_filename);

                if( exist(obj.recording_name, 'file') )
                    
                    obj.rec_path = [obj.rec_path filesep];
                    
                else
                    
                    obj.rec_path = [obj.rec_path filesep];
                    
                    if( ~exist(obj.rec_path, 'dir') )
                        obj.rec_path = [pwd filesep obj.rec_path];
                    end
                    
                    aux_files = dir([obj.rec_path obj.rec_filename '.*' ]);
                    
                    if( isempty(aux_files) )
                        error( 'ECGwrapper:ArgCheck:RecNotFound', 'Can not find recording : %s', obj.recording_name );
                    else
                        % the bigger probably would have the ECG data
                        [~, aux_idx] = max(cell2mat({aux_files(:).bytes}));
                        obj.recording_name = [obj.rec_path aux_files(aux_idx).name];
                    end
                    
                    [~, obj.rec_filename] = fileparts(obj.recording_name);
                    
                end

                
                if( strcmpi(obj.recording_format, 'auto') )
                    %Try guessing the ECG format
                    aux_fmt = ECGformat(obj.recording_name);

                    if( isempty(aux_fmt) )
                       str_aux = disp_option_enumeration( sprintf('Could not guess the format of:\n\n%s\n\nChoose one of the following:', strrep(obj.recording_name, '\', '\\')), obj.cKnownFormats);
                       error( 'ECGwrapper:ArgCheck:InvalidFormat', str_aux);
                    else
                        obj.recording_format = aux_fmt;
                    end

                end

                if( strcmpi(obj.recording_format, 'MIT') )
                    strAnnExtension = {'ari' 'atr' 'ecg'};
                    annFileName = {};
                    bAnnotationFound = false;
                    for ii = 1:length(strAnnExtension)
                        aux_str = [aux_filename '.' strAnnExtension{ii}];
                        if( exist(aux_str, 'file') )
                            annFileName = [ annFileName aux_str ];
                            bAnnotationFound = true;
                        end
                    end
                    if(~bAnnotationFound)
                        ann_aux = [];
                    else
                        for aux_str = annFileName
                            ann_aux = readannot(aux_str{1});
                            if( ~isempty(ann_aux) )
                                break;                                
                            end
                        end
                    end

                    header_aux = readheader([aux_filename '.hea']);

                elseif( strcmpi(obj.recording_format, 'ISHNE') )
                    annFileName = [aux_filename '.ann'];
                    if( exist(annFileName, 'file') )
                        ann_aux = read_ishne_ann(annFileName);
                    else
                        ann_aux = [];
                    end        
                    header_aux = read_ishne_header(obj.recording_name);

                elseif( strcmpi(obj.recording_format, 'HL7a') )

                    [~, obj.ECG_header, ann_aux, single_lead_positions ] = read_ECG(obj.recording_name, [], [], obj.recording_format);
                    
                elseif( strcmpi(obj.recording_format, 'Mortara') )
                    ann_aux = [];
                    header_aux = read_Mortara_header(obj.recording_name);
                    [obj.rec_path, obj.rec_filename] = fileparts(obj.recording_name);
                    obj.rec_path = [obj.rec_path filesep];
                    

                elseif( strcmpi(obj.recording_format, 'HES') )
                    ann_aux = read_HES_ann([aux_filename '.lst']);
                    header_aux = read_HES_header(obj.recording_name);
                    ann_aux.time = round(ann_aux.time * header_aux.freq);

                elseif( strcmpi(obj.recording_format, 'AHA') )
                    ann_aux = read_AHA_ann(obj.recording_name);
                    header_aux = read_AHA_header(obj.recording_name);

                elseif( strcmpi(obj.recording_format, 'MAT') )

                    [~, header_aux, ann_aux, single_lead_positions ] = read_ECG(obj.recording_name, [], [], obj.recording_format);
                    
                end

            end

            missing_fnames = setdiff(obj.cHeaderFieldNamesRequired, fieldnames(header_aux) );

            if( ~isempty(missing_fnames) )
                str_aux = disp_option_enumeration( '\nMissing fields in the header struct:', missing_fnames);
                str_aux = sprintf('%s\nCheck the %s for details.\n', str_aux, '<a href = "matlab: web(''http://ecg-kit.readthedocs.io/en/master/header_format.html#id1'', ''-browser'' )">ecg-kit header documentation</a>');
                error( 'ECGwrapper:ArgCheck:InvalidHeader', str_aux);
            end
            
            obj.ECG_header = header_aux;
            
            obj.bECG_rec_changed = false;
            
            obj.ECG_delineation = single_lead_positions;
            
            if(~isempty(ann_aux) && isempty(obj.ECG_annotations) )
                % discard non-beats and finish annotations parsing.
                if( isfield(ann_aux, 'anntyp') )
                    ann_aux = AnnotationFilterConvert(ann_aux, obj.recording_format, obj.class_labeling);
                end
                obj.QRS_locations = ann_aux.time;
                obj.ECG_annotations = ann_aux;
            end
            
        end

        function bOk = AreWrappersHorzcatCompatible( ECGwA, ECGwB )
        % internal method to dump signal to file in MIT format
        
            bOk = ECGwA.ECG_header.freq == ECGwB.ECG_header.freq & ...
                  ECGwA.ECG_header.nsamp == ECGwB.ECG_header.nsamp;

        end
        
        function bOk = AreWrappersVertcatCompatible( ECGwA, ECGwB )
        % internal method to dump signal to file in MIT format
        
            bOk = ECGwA.ECG_header.freq == ECGwB.ECG_header.freq & ...
                  ECGwA.ECG_header.nsig == ECGwB.ECG_header.nsig & ...
                  all(ECGwA.ECG_header.gain == ECGwB.ECG_header.gain) & ...
                  all(ECGwA.ECG_header.offset == ECGwB.ECG_header.offset);

        end
        
        function dump_wrapper(fidECG, ECGwA )
        % internal method to dump signal to file in MIT format
        
            MaxPossibleArrayBytes = 5e9;
            
            max_samples_iter = MaxPossibleArrayBytes/8/2; % double = 8 bytes
        
            aux_samples_process = ECGwA.ECG_header.nsamp * ECGwA.ECG_header.nsig;

            cant_iter =  max(1, aux_samples_process / max_samples_iter);

            cant_samples_iter = round(ECGwA.ECG_header.nsamp / cant_iter);

            all_start_idx = 1:cant_samples_iter:ECGwA.ECG_header.nsamp;

            for start_idx = all_start_idx

                % gain transform
                aux_sig = ECGwA.read_signal(start_idx, start_idx + cant_samples_iter - 1);

                fwrite(fidECG, aux_sig', 'int16', 0 );

            end
            
        end
        
        function this_header = concat_name( this_header, ECGwB, bVertConcat )
        % internal method to dump signal to file in MIT format
        
            if(bVertConcat)
                str_aux = 'vconcat_';
            else
                str_aux = 'hconcat_';
            end

            if( isempty( strfind(str_aux, this_header.recname ) ) )
                MIT_filename = [ str_aux adjust_string(this_header.recname, 20, 'center', '_')  '__' adjust_string(ECGwB.ECG_header.recname, 20, 'center', '_') ];
            else
                MIT_filename = [ this_header.recname '__' adjust_string(ECGwB.ECG_header.recname, 20, 'center', '_') ];
            end

            MIT_filename = regexprep(MIT_filename, '\W*(\w+)\W*', '$1');
            MIT_filename = regexprep(MIT_filename, '\W', '_');

            this_header.recname = MIT_filename;
        
        end
        
        function ECGwC = vertcat_internal( ECGwA, ECGwB, MIT_filename )
        % internal method to dump signal to file in MIT format
            
            ECGwC = [];
            
            if( ~AreWrappersVertcatCompatible( ECGwA, ECGwB) )
                cprintf('[1,0.5,0]', 'Wrappers are not compatible for vertical concatenation. Same gain/offset and sampling freq is necessary.\n');
            else
                % Allow signal concatenation

                this_header = ECGwA.ECG_header;
                
                % vertical concatenation
                this_header.nsamp = this_header.nsamp + ECGwB.ECG_header.nsamp;
                
                fidECG = fopen([ECGwA.output_path MIT_filename '.dat'], 'w+');

                try

                    dump_wrapper(fidECG, ECGwA );

                    dump_wrapper(fidECG, ECGwB );

                    fclose(fidECG);

                catch MEE
                    fclose(fidECG);
                    rethrow(MEE);
                end                             
                
                this_header.recname = MIT_filename;
                writeheader(ECGwA.output_path, this_header);   
                
                ECGwC = ECGwrapper( 'recording_name', [ECGwA.output_path MIT_filename '.hea'] );
                
            end
        
            
        end
        
        function ECGwC = horzcat_internal( ECGwA, ECGwB, MIT_filename )
        % internal method to dump signal to file in MIT format
            
            ECGwC = [];
            
            if( ~AreWrappersHorzcatCompatible( ECGwA, ECGwB) )
                cprintf('[1,0.5,0]', 'Wrappers are not compatible for vertical concatenation. Same gain/offset and sampling freq is necessary.\n');
            else
                % Allow signal concatenation

                this_header = ECGwA.ECG_header;
                
                aux_fnames = setdiff(fieldnames(this_header), obj.cHeaderFieldNamesRequired );
                
                for fname = rowvec(aux_fnames)
                    fname = fname{1};
                    if( ischar( this_header.(fname) ) )
                        this_header.(fname) = char([ cellstr(this_header.(fname)); cellstr(ECGwB.ECG_header.(fname))]);
                    else
                        this_header.(fname) = [this_header.(fname); ECGwB.ECG_header.(fname)];
                    end
                end

                % horizontal concatenation
                this_header.nsig = this_header.nsig + ECGwB.ECG_header.nsig;

                MaxPossibleArrayBytes = 5e9;

                max_samples_iter = MaxPossibleArrayBytes/8/2; % double = 8 bytes

                aux_samples_process = ECGwA.ECG_header.nsamp * ECGwA.ECG_header.nsig + ECGwB.ECG_header.nsamp * ECGwB.ECG_header.nsig;

                cant_iter =  max(1, aux_samples_process / max_samples_iter);

                cant_samples_iter = round(ECGwA.ECG_header.nsamp / cant_iter);

                all_start_idx = 1:cant_samples_iter:ECGwA.ECG_header.nsamp;
                
                fidECG = fopen([ECGwA.output_path MIT_filename '.dat'], 'w+');

                try

                    for start_idx = all_start_idx

                        aux_sig = ECGwA.read_signal(start_idx, start_idx + cant_samples_iter - 1);
                        
                        aux_sig = [ aux_sig ECGwB.read_signal(start_idx, start_idx + cant_samples_iter - 1) ];

                        fwrite(fidECG, aux_sig', 'int16', 0 );

                    end

                    fclose(fidECG);

                catch MEE
                    fclose(fidECG);
                    rethrow(MEE);
                end                             
                
                this_header.recname = MIT_filename;
                this_header.fname = repmat( [ECGwA.output_path MIT_filename '.dat'], this_header.nsig, 1);
                
                writeheader(ECGwA.output_path, this_header);   
                
                ECGwC = ECGwrapper( 'recording_name', [ECGwA.output_path MIT_filename '.hea'] );
                
            end
        
            
        end
        
        
        
    end
   
    
end


