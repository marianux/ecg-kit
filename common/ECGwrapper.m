classdef ECGwrapper < handle

% ECGwrapper for Matlab
% ---------------------
% 
% Description:
% 
% 
% Arguments: (specified as ECGwrapper('arg_name1', arg_val1, ... , 'arg_nameN', arg_valN) )
% 
%     a. ECG specification:
% 
%         a.1 Recording filename where ECG is stored in a valid format.
% 
%           + recording_name: ECG recording to be classified.
%           + recording_format : Valid ECG format. (MIT, ISHNE, AHA, HES, MAT)
% 
%     b. Operating modes
% 
%       
%     c. Modifiers
% 
% 
% Output:
% 
% Examples:
% 
%       Lazy users can start with:
%       
%       a2hbc
% 
%       The GUI will appear and ask you for the mandatory information.
% 
%       In case you want to try the command line, here is an example:
% 
%         a2hbc( ... 
%                 'recording_name', [ '.' filesep 'example recordings' filesep '208.dat'], ...
%                 'recording_format', 'MIT', ... 
%                 'op_mode', 'auto');
% 
%       you can also check 'examples.m' in the same folder of this file.
% 
% 
% Limits and Known bugs:
%   Probably a lot :( ... but dont panic! send me feedback if you need
%   help.
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate  : 20/3/2012
% Last update: 20/3/2012
       
    properties(GetAccess = private, Constant)
        
        cKnownFormats = {'MIT' 'ISHNE', 'AHA', 'HES', 'MAT'};
        cHeaderFieldNamesRequired = {'freq' 'nsamp' 'nsig' 'gain' 'adczero' };
        cPartitionModes = {'ECG_contiguous' 'ECG_overlapped' 'QRS'};
        cAnnotationsFieldNamesRequired = {'time' };
        cObjMethodsRequired = {'Process' 'Start' 'Concatenate'};
        cObjPropsRequired = {'progress_handle' 'name'};
        
        % maxQRSxIter: Maximum amount of heartbeats in ECG recordings of 2 leads
        maxQRSxIter = 5e3;
        % Minimum amount of heartbeats to be processed by a PID
        minQRSxIter = 1500;
        % Minimum amount of samples to be processed by a PID
        min_ECG_perPID = 432000; % 10*60*360*2 (10 minutes of 2-lead ECG@360Hz typical MITBIH arrhythmia rec)
        %gzip compression ratio
        typical_compression_ratio = 1.5; %times
        % maximum results file size
        payload_max_size = 10 * 1024^2;
        % maximum debugging email report file size
        max_report_filesize = 4 * 1024^2; % bytes
        % Maximum time to wait for other PIDs finishing their work.
        Time2WaitPIDs = 15 * 60; % seconds.
        
    end
    
    properties ( Access = private )
        rec_filename
        bHaveUserInterface
        Loops2io
        MaxNodesReading
        MaxNodesWriting
        common_path
        bQRSlocations
        QRS_locations
        bArgChanged = true;
        bECG_rec_changed = true;
        bCreated = false;
        % maxECGxIter: Maximum samples ECG
        maxECGxIter
        
    end
    
    properties (SetAccess = private, GetAccess = public)  
        %read-only 
        Processed = false;
        Result_files = []; 
        doPayload = true;
    end

    properties
        %read-write
        tmp_path
        recording_name
        recording_format
        ECG_annotations 
        ECG_header        
        cant_pids
        this_pid
        repetitions
        ECGtaskHdl
        partition_mode
        overlapping_time
        cacheResults
        syncSlavesWithMaster
        class_labeling = 'AAMI'
    end
    

    methods 
        
        function obj = ECGwrapper(varargin)

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

            %% Argument parsing

            %argument definition
            p = inputParser;   % Create instance of inputParser class.
            p.addParamValue('recording_name', [], @(x)( ischar(x)));
            p.addParamValue('recording_format', [], @(x)( ischar(x) && any(strcmpi(x,obj.cKnownFormats))) || isempty(x) );
            p.addParamValue('cant_pids', 1, @(x)(isnumeric(x) && x > 0 ) );
            p.addParamValue('this_pid', 1, @(x)(isnumeric(x) && x > 0 ) );
            p.addParamValue('tmp_path', tempdir, @(x)(ischar(x)) );
            p.addParamValue('repetitions', 1, @(x)(isnumeric(x) && x > 0 ) );
            p.addParamValue('ECGtaskHandle', ECGtask_do_nothing, @(x)( isobject(x) ) );
            p.addParamValue('overlapping_time', 30, @(x)(isnumeric(x) && x > 0 ) );
            p.addParamValue('partition_mode', 'ECG_overlapped', @(x)( ischar(x) && any(strcmpi(x,obj.cPartitionModes))) );
            p.addParamValue('cacheResults', true, @(x)(islogical(x)) );
            p.addParamValue('syncSlavesWithMaster', true, @(x)(islogical(x)) );

            try
                p.parse( varargin{:} );
            catch MyError
                rethrow(MyError);    
            end

            obj.recording_name = p.Results.recording_name;
            obj.recording_format = p.Results.recording_format;
            obj.cant_pids = p.Results.cant_pids;
            obj.this_pid = p.Results.this_pid;
            obj.tmp_path = p.Results.tmp_path;
            obj.repetitions = p.Results.repetitions;
            obj.ECGtaskHdl = p.Results.ECGtaskHandle;
            obj.overlapping_time = p.Results.overlapping_time;
            obj.partition_mode = p.Results.partition_mode;
            obj.cacheResults = p.Results.cacheResults;
            obj.syncSlavesWithMaster = p.Results.syncSlavesWithMaster;

            
            % Dont know why this variable uses a lot of bytes to store at disk.
            clear p
            
            obj.bCreated = true;
            
        end
        
        function Run(obj)

            result_files = {};
            repeat_idx = 1;

            if( obj.bArgChanged )
                obj = CheckArguments(obj);
                obj.bArgChanged = false;
            end
            
            overlapping_samples = obj.overlapping_time * obj.ECG_header.freq;
            
            while ( repeat_idx <= obj.repetitions )

                try

                    %% work starts here

                    %% Prepare jobs to perform.

                    % Activate the progress_struct bar.
                    pb = progress_bar([obj.rec_filename ': ' obj.ECGtaskHdl.name ]);
                    
                    % Initialization of the process.
                    obj.ECGtaskHdl.progress_handle = pb;
                    obj.ECGtaskHdl.Start( obj.ECG_header, obj.ECG_annotations );
                    
                    % find out the ammount of memory in the system
                    user = memory();
                    obj.maxECGxIter = round(obj.ECGtaskHdl.memory_constant*((user.MaxPossibleArrayBytes)/8)); % double = 8 bytes
                    
                    if( strcmpi(obj.partition_mode, 'QRS') )
                        %% QRS mode

                        cant_QRS_locations = length(obj.QRS_locations);

                        %PID parsing
                        if( obj.cant_pids > 1 )

                            max_recommended_cant_pids = max(1, round( cant_QRS_locations / obj.minQRSxIter ));

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
                        cant_samples= obj.QRS_locations(QRS_end_idx) - obj.QRS_locations(QRS_start_idx);
                        cant_iter = ceil(cant_samples * obj.ECG_header.nsig / obj.maxECGxIter);
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
                            cant_iter = length(iter_starts);
                                                        
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
                        cant_iter = ceil(cant_samples2do * obj.ECG_header.nsig / obj.maxECGxIter );
                        [iter_starts, iter_ends] = TaskPartition( cant_samples2do, cant_iter);

                    end

                    if( obj.this_pid > obj.cant_pids )
                        %este pid no trabaja.
                        error('ECGwrapper:TooMuchPIDs', [ 'This PID ' num2str(obj.this_pid) ' will not work, consider decreasing obj.cant_pids to ' num2str(obj.cant_pids) '.\n' ] );
                    end

                    % Update point
                    pb.checkpoint('Initializing');

                    %% Iterations over the whole ECG recording

                    start_iter = 1;

                    %check for previous iterations already done, and try to restore.
                    for ii = 1:cant_iter
                        aux_filename =  [obj.tmp_path 'tmpfile_' obj.ECGtaskHdl.name '_' obj.rec_filename '_payload_cantpids_' num2str(obj.cant_pids) '_thispid_' num2str(obj.this_pid) '_iteration_' num2str(ii) '_of_' num2str(cant_iter)  '.mat'];
                        if( exist(aux_filename, 'file') )
                            start_iter = ii+1;
                        else
                            break
                        end
                    end

                    pb.Loop2do = cant_iter;
                    
                    for this_iter = start_iter:cant_iter
                        %% loop initialization

                        % por q deberia saberlo ?
                        %obj.ECGtaskHdl.this_iteration = this_iter;

                        %start of the progress_struct loop 0%
                        pb.start_loop();

                        if( strcmpi(obj.partition_mode, 'QRS') )
                            % la partición se hizo en la cantidad de latidos

                            this_iter_QRS_start_idx = QRS_start_idx - 1 + iter_starts(this_iter);
                            this_iter_QRS_end_idx = QRS_start_idx - 1 + iter_ends(this_iter);

                            this_iter_ECG_start_idx = max(1, obj.QRS_locations(this_iter_QRS_start_idx) - overlapping_samples);
                            this_iter_ECG_end_idx = min(obj.ECG_header.nsamp, obj.QRS_locations(this_iter_QRS_end_idx) + overlapping_samples);
                            
                            % create an annotation struct for this
                            % iteration.
                            this_ann = obj.ECG_annotations;
                            for field_names = rowvec(fieldnames(this_ann))
                                if( ~isempty( this_ann.(field_names{1}) ) )
                                    this_ann.(field_names{1}) = this_ann.(field_names{1})(this_iter_QRS_start_idx:this_iter_QRS_end_idx);
                                end
                            end
                            this_ann.time = this_ann.time - this_iter_ECG_start_idx + 1;

                            % in QRS mode, this is not useful.
                            this_iter_ECG_relative_start_end_idx = [1 (this_iter_ECG_end_idx - this_iter_ECG_start_idx + 1)];
                            
                        else
                            % la partición se hizo en el registro de ECG 

                            this_iter_QRS_start_idx = nan;
                            this_iter_QRS_end_idx = nan;
                            this_iter_cant_QRS = nan;

                            this_iter_ECG_start_idx = max(1, ECG_start_idx - 1 + iter_starts(this_iter) - overlapping_samples);
                            this_iter_ECG_end_idx = min(obj.ECG_header.nsamp, ECG_start_idx - 1 + iter_ends(this_iter) + overlapping_samples);

                            %No annotations provided in this mode
                            this_ann = [];
                            
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

                        %% read from file

                        ECG = read_ECG(obj.recording_name, obj.recording_format, this_iter_ECG_start_idx, this_iter_ECG_end_idx );

                        %convert ADC samples (int16) to obj.ECGtaskHdl.target_units volts.
                        [ECG this_header] = ADC2units(double(ECG), this_header, obj.ECGtaskHdl.target_units);

                        %% User defined function calculation

                        % Update point
                        pb.checkpoint('User function');

%                         if( strcmpi(obj.partition_mode, 'QRS') )
%                             fprintf(1, 'ECG from: %d - %d\n', obj.QRS_locations(this_iter_QRS_start_idx) , obj.QRS_locations(this_iter_QRS_end_idx) );
%                         else
%                             fprintf(1, 'ECG from: %d - %d\n', this_iter_ECG_start_idx + this_iter_ECG_relative_start_end_idx - 1);
%                         end

                        %user-defined execution according to 'process'
                        payload = obj.ECGtaskHdl.Process(ECG, this_iter_ECG_relative_start_end_idx, this_header, this_ann, [this_iter_QRS_start_idx this_iter_QRS_end_idx]);
                        
                        if( obj.ECGtaskHdl.doPayload )
                            
                            % Update point
                            pb.checkpoint('Saving results');

                            % save results
                            save([obj.tmp_path 'tmpfile_' obj.ECGtaskHdl.name '_' obj.rec_filename '_payload_cantpids_' num2str(obj.cant_pids) '_' num2str(obj.this_pid) '_' num2str(this_iter) '_' num2str(cant_iter) '.mat'], '-struct', 'payload');


                            % rename results.
                            movefile(   [obj.tmp_path 'tmpfile_' obj.ECGtaskHdl.name '_' obj.rec_filename '_payload_cantpids_' num2str(obj.cant_pids) '_' num2str(obj.this_pid) '_' num2str(this_iter) '_' num2str(cant_iter)  '.mat'], ...
                                        [obj.tmp_path 'tmpfile_' obj.ECGtaskHdl.name '_' obj.rec_filename '_payload_cantpids_' num2str(obj.cant_pids) '_thispid_' num2str(obj.this_pid) '_iteration_' num2str(this_iter) '_of_' num2str(cant_iter)  '.mat'], 'f' );

                        end
                        
                        pb.end_loop();
                                
                    end
                    
                    %indicate end of loop.
                    pb.reset();

                    %clear some not needed variables 
                    clear *ECG this_iter* ECG_annotations

                    if( obj.ECGtaskHdl.doPayload )
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

                                    %feature matrices
                                    payload = [];
                                    payload_idx = 1;

                                    % to resume after the last dump to disk.
                                    if( resume_iter_this_pid > 0 )
                                        start_iter_this_pid = resume_iter_this_pid;
                                    end

                                    for ii = start_pid:obj.cant_pids

                                        files_this_pid = dir([obj.tmp_path 'tmpfile_' obj.ECGtaskHdl.name '_' obj.rec_filename '_payload_cantpids_' num2str(obj.cant_pids) '_thispid_' num2str(ii) '_*.mat' ]);

                                        if( isempty(files_this_pid) )
                                            error('ECGwrapper:PIDnotFinished', 'Handled error');
                                        end

                                        cant_iteraciones_this_pid = length(files_this_pid);

                                        % por q deberia saberlo ?
                                        %obj.ECGtaskHdl.cant_iteration = cant_iteraciones_this_pid;

                                        for jj = start_iter_this_pid:cant_iteraciones_this_pid 

                                            % por q deberia saberlo ?
                                            %obj.ECGtaskHdl.this_iteration = jj;

                                            aux_FM_filename = [obj.tmp_path 'tmpfile_' obj.ECGtaskHdl.name '_' obj.rec_filename '_payload_cantpids_' num2str(obj.cant_pids) '_thispid_' num2str(ii) '_iteration_' num2str(jj) '_of_' num2str(cant_iteraciones_this_pid) '.mat' ];
                                            if( ~exist(aux_FM_filename, 'file'))
                                                %Probably this PID not finished yet.
                                                error('ECGwrapper:PIDnotFinished', 'Handled error');
                                            end
                                            aux = load(aux_FM_filename);

                                            % Process the results according to the
                                            % 'Concatenate' method
                                            payload = obj.ECGtaskHdl.Concatenate(payload, aux);

                                            %check variable growth
                                            payload_var_size = whos('payload');
                                            payload_var_size = payload_var_size.bytes/obj.typical_compression_ratio;

                                            if( payload_var_size > obj.payload_max_size )
                                                %force dump to disk.
                                                payload.payload_idx = payload_idx;
                                                save([obj.tmp_path 'payloads_' obj.ECGtaskHdl.name '_' obj.rec_filename '_' num2str(payload_dump_counter) '.mat'], '-struct', 'payload');
                                                %update counter and reset payload.
                                                payload_idx = payload_idx + size(payload,1) + 1;
                                                payload = [];
                                                payload_dump_counter = payload_dump_counter + 1;
                                                %to continue from this pid/iteration.
                                                start_pid = ii;
                                                resume_iter_this_pid = jj+1;
                                            end

                                        end

                                        %the first time in this loop we resume after the last
                                        %file dump, and never again.
                                        start_iter_this_pid = 1;

                                    end

                                    if( ~isempty(payload) )
                                        % data remain in memory -> dump 2 disk
                                        payload.payload_idx = payload_idx;
                                        save([obj.tmp_path 'payloads_' obj.ECGtaskHdl.name '_' obj.rec_filename '_' num2str(payload_dump_counter) '.mat'], '-struct', 'payload');
                                        %update counter and reset payload.
                                        payload = [];
                                        payload_dump_counter = payload_dump_counter + 1;
                                    end

                                    % Update point
                                    pb.checkpoint('Generating final results');

                                    %rename results to put some ordering.
                                    files_this_pid = dir([obj.tmp_path 'payloads_' obj.ECGtaskHdl.name '_' obj.rec_filename '_*.mat']);
                                    cant_iteraciones_this_pid = length(files_this_pid);

                                    for jj = 1:cant_iteraciones_this_pid

                                        if( obj.repetitions > 1)
                                            auxStr = [obj.tmp_path obj.ECGtaskHdl.user_string '_' obj.ECGtaskHdl.name '_' obj.rec_filename '_part_' num2str(jj) '_of_' num2str(cant_iteraciones_this_pid) '_repetition_' num2str(repeat_idx) '.mat'];
                                        else
                                            auxStr = [obj.tmp_path obj.ECGtaskHdl.user_string '_' obj.ECGtaskHdl.name '_' obj.rec_filename '_part_' num2str(jj) '_of_' num2str(cant_iteraciones_this_pid) '.mat'];
                                        end

                                        movefile(   [obj.tmp_path 'payloads_' obj.ECGtaskHdl.name '_' obj.rec_filename '_' num2str(jj) '.mat'], ...
                                                    auxStr, 'f' );

                                        result_files = [ result_files; cellstr(auxStr)];

                                    end

                                    % Update point
                                    pb.checkpoint('Deleting temporal files');

                                    %clean temporal files
                                    delete([obj.tmp_path 'tmpfile_' obj.ECGtaskHdl.name '_' obj.rec_filename '_payload_cantpids_' num2str(obj.cant_pids) '*.mat' ]);

                                    bContinue = false;


                                catch ME
                                    % Error handling

                                    if( strfind(ME.identifier, 'ECGwrapper') )
                                        if( obj.cant_pids == 1 || toc(Start2Wait) > obj.Time2WaitPIDs )
                                            error('ECGwrapper:PIDnotFinished', 'Timeout. Master give up waitng Slaves.');
                                        end
                                        pause(30);
                                    else
                                        rethrow(ME)
                                    end
                                end


                            end

                        else
                        %% Slave PIDs

                            if( obj.syncSlavesWithMaster )
                            % other PIDS sync here after Master build the
                            % results file/s.

                                bContinue = true;
                                %Wait for Time2WaitPIDs seconds the finalization of all PIDs. Otherwise exit
                                %with error.
                                Start2Wait = tic();

                                while(bContinue)

                                    try

                                        files_this_pid = dir([obj.tmp_path 'tmpfile_' obj.ECGtaskHdl.name '_' obj.rec_filename '_payload_cantpids_' num2str(obj.cant_pids) '_thispid_' num2str(obj.this_pid) '_*.mat' ]);

                                        if( ~isempty(files_this_pid) )
                                            error('ECGwrapper:MasterNotFinished', 'Handled error');
                                        end


                                        bContinue = false;

                                    catch ME

                                        if( strfind(ME.identifier, 'ECGwrapper') )
                                            if( toc(Start2Wait) > 1.5*obj.Time2WaitPIDs )
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
                        obj.ECGtaskHdl.Finish( obj.ECG_header );
                    end

                catch MException
                    %% Error handling

                    if( obj.bHaveUserInterface)
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
                                    %other home-made errors. Make an educated exit ...
                                    rethrow(MException)
                                end
                            end

                        else
                            
                            %% Other unknown errors
                            rethrow(MException);
                            
                        end
                        
                    else
                        %% No User interface
                        fprintf(2,'\n\n')
                        fprintf(2,'###########\n')
                        fprintf(2,'## ERROR ##\n')
                        fprintf(2,'###########\n')

                        report = getReport(MException);
                        
                        fprintf(2, '%s', report);
                        fprintf(2,'###########\n')
                        fprintf(2,'## ERROR ##\n')
                        fprintf(2,'###########\n')

                        rethrow(MException)

                    end


                end

                repeat_idx = repeat_idx + 1;
                
                % Update point
                pb.checkpoint('Work done.');

            end

            obj.Processed = true;
            obj.Result_files = result_files;
            
            % destroy the progress bar.
            pb.delete;
            
        end
        
        function obj = CheckECGrecording(obj)

            if( isempty(obj.recording_name) )
                error( 'ECGwrapper:ArgCheck:InvalidECGarg', 'Please provide an ECG recording as described in the documentation, help(''ECGwrapper'') maybe could help you.\n' );
            else

                % ECG to be read

                [~, obj.rec_filename] = fileparts(obj.recording_name);

                if( ~exist(obj.recording_name, 'file') )
                   error( 'ECGwrapper:ArgCheck:RecNotFound', 'Can not find recording : %s', obj.recording_name );
                end

                if( strcmp(obj.recording_format, 'auto') )
                    %Try guessing the ECG format
                    obj.recording_format = ECGformat(obj.recording_name);

                    if( isempty(obj.recording_format) )
                       strAux = [ repmat(' + ', length(obj.cKnownFormats), 1) char(obj.cKnownFormats) repmat('\n', length(obj.cKnownFormats), 1 ) ];
                       error( 'ECGwrapper:ArgCheck:InvalidFormat', ['Unknown format. Choose one of the following:\n' rowvec(strAux') ] );
                    end

                end

                if( strcmp(obj.recording_format, 'MIT') )
                    strAnnExtension = {'ari' 'atr' 'ecg'};
                    bAnnotationFound = false;
                    for ii = 1:length(strAnnExtension)
                        annFileName = [obj.recording_name(1:end-3) strAnnExtension{ii}];
                        if( exist(annFileName, 'file') )
                            bAnnotationFound = true;
                            break
                        end
                    end
                    if(~bAnnotationFound)
                        ann_aux = [];
                    else
                        ann_aux = readannot(annFileName);
                    end

                    obj.ECG_header = readheader([obj.recording_name(1:end-3) 'hea']);

                elseif( strcmp(obj.recording_format, 'ISHNE') )
                    annFileName = [obj.recording_name(1:end-3) 'ann'];
                    if( exist(annFileName, 'file') )
                        ann_aux = read_ishne_ann(annFileName);
                    else
                        ann_aux = [];
                    end        
                    obj.ECG_header = read_ishne_header(obj.recording_name);

                elseif( strcmp(obj.recording_format, 'HES') )
                    ann_aux = read_HES_ann([obj.recording_name(1:end-3) 'lst']);
                    obj.ECG_header = read_HES_header(obj.recording_name);
                    ann_aux.time = round(ann_aux.time * obj.ECG_header.freq);

                elseif( strcmp(obj.recording_format, 'AHA') )
                    ann_aux = read_AHA_ann(obj.recording_name);
                    obj.ECG_header = read_AHA_header(obj.recording_name);

                elseif( strcmp(obj.recording_format, 'MAT') )

                    aux_load = load(obj.recording_name, { 'heasig', 'ann' });

                    if( ~isfield( aux_load, 'heasig') )
                        error( 'ECGwrapper:ArgCheck:InvalidECGarg', 'MAT file does not include the "heasig" variable.\n' );
                    end

                    obj.ECG_header = aux_load.obj.ECG_header;

                    if( ~isfield( aux_load, 'ann') )
                        ann_aux = aux_load.ECG_annotations;
                    else
                        ann_aux = [];
                    end

                    clear aux_load

                end
                
                obj.ECG_header.ECG_format = obj.recording_format;

            end

            if(isempty(ann_aux))
                obj.bQRSlocations = false;
                obj.ECG_annotations = [];
                obj.QRS_locations = [];
            else
                % discard non-beats and finish annotations parsing.
                obj.bQRSlocations = true;
                ann_aux = AnnotationFilterConvert(ann_aux, obj.ECG_header.ECG_format, obj.class_labeling);
                obj.QRS_locations = ann_aux.time;
                obj.ECG_annotations = ann_aux;
            end
            
            obj.bECG_rec_changed = false;
            
        end

        
        function disp(obj)
            tab_size = 60;
            if( obj.bCreated )
                if( isobject(obj.ECGtaskHdl) )
                    task_name = obj.ECGtaskHdl.name;
                else
                    task_name = 'Undef';
                end
                
                fprintf(1,[ '+--------------------------+\n' ...
                            '| ECGwrapper object config |\n' ...
                            '+--------------------------+\n' ...
                            '+ECG recording: %s (%s) \t +PID: %d/%d\n' ... 
                            '+Function name: %s '   '\t +repetitions: %d\n' ...
                            '+Partition mode: %s\n'], ... 
                            adjust_string(obj.recording_name, tab_size-23), ...
                            obj.recording_format, ...
                            obj.this_pid, ...
                            obj.cant_pids, ...
                            adjust_string(task_name, tab_size-16), ...
                            obj.repetitions, ...
                            adjust_string(obj.partition_mode, tab_size));

                if( obj.Processed)
                    strAux = '+Processed: true';
                else
                    strAux = '+Processed: false';
                end
                fprintf(1, [ strAux repmat(' ', 1, tab_size-length(strAux)) '\t +TMP: %s\n'], adjust_string(obj.tmp_path, 20));
                
                if( obj.Processed)
                    fprintf(1, ['\n' ... 
                                'Result files\n' ...
                                '------------\n'] );
                    disp(obj.Result_files);
                end
                
            else
                fprintf(1, 'Object not created yet.');
            end
        end
        
        function set.tmp_path(obj,value)
            if( ischar(value) || isempty(value) )
                obj.tmp_path = value;
                obj.bArgChanged = true;
            else
                warning('ECGwrapper:BadArg', 'tmp_path must be a string.');
            end
        end
        
        function set.recording_name(obj,value)
            if( ischar(value) )
                obj.recording_name = value;
                obj.bArgChanged = true;
                obj.bECG_rec_changed = true;
            else
                warning('ECGwrapper:BadArg', 'recording_name must be a string.');
            end
        end

        function set.recording_format(obj,value)
            if( ischar(value) )
                obj.recording_format = value;
                obj.bArgChanged = true;
                obj.bECG_rec_changed = true;
            elseif( isempty(value) )
                obj.recording_format = 'auto';
                obj.bArgChanged = true;
            else
                warning('ECGwrapper:BadArg', 'recording_format must be a string.');
            end
        end

        function set.ECG_annotations(obj,value)
            
            if( isstruct(value) )
                
                if( all(isfield(value, obj.cAnnotationsFieldNamesRequired )) )
                    
                    obj.ECG_annotations = value;
                    obj.bQRSlocations = true;
                    obj.QRS_locations = value.time;
                
                else
                    
                    strAux = [ repmat(' + ', length(obj.cAnnotationsFieldNamesRequired), 1) char(obj.cAnnotationsFieldNamesRequired) repmat('\n', length(obj.cAnnotationsFieldNamesRequired), 1 ) ];
                    warning( 'ECGwrapper:BadArg', ['Please provide the following fields in the annotations struct:\n ' rowvec(strAux') ] );
                end
                
            else
                warning('ECGwrapper:BadArg', 'ECG_annotations must be a struct.');
            end
        end

        function value = get.ECG_annotations(obj)
            
            if( obj.bECG_rec_changed )
                obj = obj.CheckECGrecording();
            end
            
            value = obj.ECG_annotations;
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
            if( value > 0 && value <= obj.cant_pids)
                obj.this_pid = value;
                obj.bArgChanged = true;
            else
                warning('ECGwrapper:BadArg', 'this_pid must be <= cant_pids.');
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

        function set.ECGtaskHdl(obj,value)
            if( isobject(value) )
                obj.ECGtaskHdl = value;
                obj.bArgChanged = true;
            else
                warning('ECGwrapper:BadArg', 'ECGtaskHdl must be a ECGtask object.');
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
   
    
end

function obj = CheckArguments(obj)

    %Object parsing
    if( isobject(obj.ECGtaskHdl) )
        for ii = 1:length(obj.cObjMethodsRequired)
            if( ~ismethod(obj.ECGtaskHdl, obj.cObjMethodsRequired{ii}) )
               error( 'ECGwrapper:ArgCheck:UserObjHdl', ['Method ' obj.cObjMethodsRequired{ii} ' not implemented in UserObjHdl.\n\n'] );
            end        
        end

        for ii = 1:length(obj.cObjPropsRequired)
            if( ~isprop(obj.ECGtaskHdl, obj.cObjPropsRequired{ii}) )
               error( 'ECGwrapper:ArgCheck:UserObjHdl', ['Property ' obj.cObjPropsRequired{ii} ' not present in UserObjHdl.\n\n'] );
            end        
        end
    else
       error( 'ECGwrapper:ArgCheck:UserObjHdl', 'ECGtaskHdl is not a valid ECGtask handle.' );
    end

    %ECG parsing
    obj = CheckECGrecording(obj);
    

    if( isempty(obj.tmp_path) )
        if(ispc())
            obj.tmp_path = tempdir ;
        else
            obj.tmp_path = './tmp/';
        end
    %         obj.tmp_path = [fileparts(mfilename('fullpath')) filesep 'tmp' filesep ];
    %     obj.tmp_path = [fileparts(obj.recording_name) filesep ];
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
        if( obj.bQRSlocations )
            if( obj.overlapping_time < 5 ) %seconds
                warning('ECGwrapper:ArgCheck:Overlapp_too_low', 'The overlapping time between iterations is too low, consider increasing.\n' );
            end
        else
            error('ECGwrapper:ArgCheck:QRS_det_not_available', 'Please provide valid QRS detections to use this mode.\n' );
        end

    else
        if( strcmpi(obj.partition_mode, 'ECG_contiguous') )
            %One segment after the other.
            obj.overlapping_time = 0;            
        end
    end

end

