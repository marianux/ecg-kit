function examples(pid_str, examples_path, user_str)

% Example of how to use the ECGkit
% ---------------------------------
% 
% Description:
% 
% This script exemplifies the use of the most important features of the
% kit, which are:
% 
%  + QRS detection algorithms: (Wavedet - Pan&Tompkins - gqrs)
%  + QRS manual corrector graphical user interface (GUI)
%  + ECG waves delineator algorithm (Wavedet)
%  + Time-series manual corrector graphical user interface (GUI)
%  + Heartbeat classifier (a2hbc)
%  + ABP/PPG pulse detection/delineation ( WavePPG - wabp )
%  + Signal report generator
% 
%  This script is prepared to run locally without arguments, as well as in
%  a cluster environment by using "pid_str" argument. The pid_str
%  argument is a char with format 'N/M', being N <= M with default value
%  '1/1'. You can partition a big job into M pieces in cluster
%  architecture, by starting M processes with N ranging from 1 to M.
% 
%  You can watch a typical run of this script for small, local ECG
%  recording here:
% 
% 
%  or for several big recordings here:
% 
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate  : 11/2/2014
% Last update: 18/7/2014


%% Argument parsing

    if( nargin < 1 || ~ischar(pid_str) )
        % single PID run
        pid_str = '1/1';
    end
    
    if( nargin < 2 || ~exist(examples_path, 'dir') )
        % inspect ECG files in rootpath\example recordings\ folder
        root_path = fileparts(mfilename('fullpath'));
        % default folder to look at
        examples_path = [root_path filesep 'example recordings' filesep ];
        if(~exist(examples_path, 'dir'))
            disp_string_framed(2, 'Please provide a valid path with ECG recordings');
            return
        end
    else
        if( examples_path(end) ~= filesep )
            examples_path = [examples_path filesep];
        end
    end
    
    if( nargin < 3  )
        user_str = '';
    end
    
    %% Start
    
    % This script exemplifies the use of the toolbox in several typical
    % tasks
    
    filenames = dir(examples_path);
    recnames = {filenames(:).name};
    [~,recnames] = cellfun(@(a)(fileparts(a)), recnames, 'UniformOutput', false);
    recnames = unique(recnames);
    recnames = setdiff(recnames, {'' '.' '..' 'results' 'condor' });
    recnames = {'ex_ABP_PPG_Registro_01M'};
%     recnames = recnames(1)
    lrecnames = length(recnames);

    % In case of running in a user-assisted fashion.
    bUseDesktop = usejava('desktop');

    if( bUseDesktop )
        tmp_path = tempdir;
        output_path = [ examples_path 'results' filesep ];
    else
        % For cluster or distributed environment processing.
        
        InstallECGkit();
        
        % this is a local path, usually faster to reach than output_path
        tmp_path = '/scratch/'; 
        
        % distributed or cluster-wide accesible path
        output_path = [ examples_path 'results' filesep ];
    end
    
% just for debugging, keep it commented.      
%     bUseDesktop = false
    
    %% QRS automatic detection
    
    % go through all files
    
    ECG_all_wrappers = [];
    jj = 1;
    
    for ii = 1:lrecnames
        
        rec_filename = [examples_path recnames{ii}];
        
        % task name, 
%         ECGt_QRSd = 'QRS_detection';
        % or create an specific handle to have more control
        ECGt_QRSd = ECGtask_QRS_detection(); 
%         % select an specific algorithm. Default: Run all detectors
%         ECGt_QRSd.detectors = 'wavedet'; % Wavedet algorithm based on
%         ECGt_QRSd.detectors = 'pantom';  % Pan-Tompkins alg.
%         ECGt_QRSd.detectors = 'gqrs';    % WFDB gqrs algorithm.
%         ECGt_QRSd.detectors = 'user:example_worst_ever_QRS_detector';    % Example of how you can add your own QRS detector.
%         ECGt_QRSd.detectors = 'user:your_QRS_detector_func_name';    %
%         "your_QRS_detector_func_name" can be your own detector.
        ECGt_QRSd.detectors = {'wavedet' 'gqrs' 'user:example_worst_ever_QRS_detector'};    

        % you can individualize each run of the QRS detector with an
        % external string
        ECGt_QRSd.user_string = user_str;
        % or group by the config used
%         ECGt_QRSd.user_string = ECGt_QRSd.detectors;
        
        
%         ECGt_QRSd.only_ECG_leads = false;    % consider all signals ECG
        ECGt_QRSd.only_ECG_leads = true;    % Identify ECG signals based on their header description.
        
        ECG_w = ECGwrapper( 'recording_name', rec_filename, ...
                            'this_pid', pid_str, ...
                            'tmp_path', tmp_path, ...
                            'output_path', output_path, ...
                            'ECGtaskHandle', ECGt_QRSd);
        
        try
            
            % process the task
            ECG_w.Run;
            
            % collect object if were recognized as ECG recordings.
            if( jj == 1)
                ECG_all_wrappers = ECG_w;
            else
                ECG_all_wrappers(jj) = ECG_w;
            end
            jj = jj + 1;
            
        catch MException
            
            if( strfind(MException.identifier, 'ECGwrapper:ArgCheck:InvalidFormat') )
                disp_string_framed('*Red', sprintf( 'Could not guess the format of %s', ECG_w.recording_name) );
            else
                % report just in case 
                report = getReport(MException);
                fprintf(2, '\n%s\n', report);
            end
        end
    
    end
    
    % recognized recordings
    lrecnames = length(ECG_all_wrappers);
    
    % at the end, report problems if happened.
    for ii = 1:lrecnames
        ECG_all_wrappers(ii).ReportErrors;
    end
    
    %% QRS visual inspection and correction

    if( bUseDesktop )
        
        % other task can be performed on the same objects
        for ii = 1:lrecnames

            % last worker is the responsible of the visual correction.
            if( ECG_all_wrappers(ii).this_pid == ECG_all_wrappers(ii).cant_pids)

                % if there are not any previous error.
                if( ECG_all_wrappers(ii).Processed && ~ECG_all_wrappers(ii).Error ) 

                    % this is to use previous saved results as starting point,
                    % if any available
                    cached_filenames = ECG_all_wrappers(ii).GetCahchedFileName({'QRS_corrector' 'QRS_detection'});

                    % if no previous correction work, try the automatic
                    % detection task
                    
                    % if any, do the correction
                    if( ~isempty(cached_filenames) )

                        % this is to use previous saved results as starting point,
                        % if any available
                        ECG_all_wrappers(ii).ECGtaskHandle = 'QRS_corrector';
                        
                        % This task is supposed to be supervised, so only one pid is enough.
                        ECG_all_wrappers(ii).this_pid = '1/1';

                        % user provided name to individualize each run
                        ECG_all_wrappers(ii).ECGtaskHandle.user_string = user_str;
                        
                        % to avoid loading cached results and exit, this flag
                        % allows the re-editing of the current state of the
                        % detections.
                        ECG_all_wrappers(ii).cacheResults = false;

                        % maybe in your application you should run this for
                        % all files.
                        ECG_all_wrappers(ii).ECGtaskHandle.payload = load(cached_filenames{1});

                        % process the task
                        ECG_all_wrappers(ii).Run;

                        % restore the original pids configuration
                        ECG_all_wrappers(ii).this_pid = pid_str;

                        % As we changed for "QRS correction" task, we have to enable this
                        % value again in order to avoid performing the following tasks every time.
                        % If you want to recalculate any task, change it to false
                        ECG_all_wrappers(ii).cacheResults = true;
                        
                    end


                end

            end

        end

        % at the end, report problems if happened.
        for ii = 1:lrecnames
            ECG_all_wrappers(ii).ReportErrors;
        end

    end
    
    %% PPG/ABP pulse detection
    
    % other task can be performed on the same objects
    for ii = 1:lrecnames

        % set the delineator task name and run again.
        ECG_all_wrappers(ii).ECGtaskHandle = 'PPG_ABP_detector';

        % user provided name to individualize each run
        ECG_all_wrappers(ii).ECGtaskHandle.user_string = user_str;

        % process the task
        ECG_all_wrappers(ii).Run;
        
    end

    % at the end, report problems if happened.
    for ii = 1:lrecnames
        ECG_all_wrappers(ii).ReportErrors;
    end
    
    %% PPG/ABP waves visual inspection and correction

    if( bUseDesktop )
        
        % other task can be performed on the same objects
        for ii = 1:lrecnames

            % last worker is the responsible of the visual correction.
            if( ECG_all_wrappers(ii).this_pid == ECG_all_wrappers(ii).cant_pids)

                % if there are not any previous error.
                if( ECG_all_wrappers(ii).Processed && ~ECG_all_wrappers(ii).Error ) 

                    % this is to use previous saved results as starting point,
                    % if any available
                    cached_filenames = ECG_all_wrappers(ii).GetCahchedFileName({'PPG_ABP_corrector' 'PPG_ABP_detector'});

                    % if no previous correction work, try the automatic
                    % detection task
                    
                    % if any, do the correction
                    if( ~isempty(cached_filenames) )

                        % this is to use previous saved results as starting point,
                        % if any available
                        ECG_all_wrappers(ii).ECGtaskHandle = 'PPG_ABP_corrector';
                        
                        % This task is supposed to be supervised, so only one pid is enough.
                        ECG_all_wrappers(ii).this_pid = '1/1';

                        % user provided name to individualize each run
                        ECG_all_wrappers(ii).ECGtaskHandle.user_string = user_str;
                        
                        % to avoid loading cached results and exit, this flag
                        % allows the re-editing of the current state of the
                        % detections.
                        ECG_all_wrappers(ii).cacheResults = false;

                        % maybe in your application you should run this for
                        % all files.
                        ECG_all_wrappers(ii).ECGtaskHandle.payload = load(cached_filenames{1});

                        % process the task
                        ECG_all_wrappers(ii).Run;

                        % restore the original pids configuration
                        ECG_all_wrappers(ii).this_pid = pid_str;

                        % As we changed for "QRS correction" task, we have to enable this
                        % value again in order to avoid performing the following tasks every time.
                        % If you want to recalculate any task, change it to false
                        ECG_all_wrappers(ii).cacheResults = true;
                        
                    end


                end

            end

        end

        % at the end, report problems if happened.
        for ii = 1:lrecnames
            ECG_all_wrappers(ii).ReportErrors;
        end

    end
        
    %% ECG automatic delineation
    
    % other task can be performed on the same objects
    for ii = 1:lrecnames
        
        % this is to use previous cached results as starting point
        cached_filenames = ECG_all_wrappers(ii).GetCahchedFileName('QRS_corrector');
        
        % if corrected QRS detections are not available, wavedet
        % performs automatic QRS detection.
        if( ~isempty(cached_filenames) )
            % this is to use previous result from the automatic QRS
            % detection
            ECG_all_wrappers(ii).ECGtaskHandle.payload = load(cached_filenames{1});

        end

        % set the delineator task name and run again.
        ECG_all_wrappers(ii).ECGtaskHandle = 'ECG_delineation';

        % user provided name to individualize each run
        ECG_all_wrappers(ii).ECGtaskHandle.user_string = user_str;

%         ECGt_QRSd.detectors = 'wavedet'; % Wavedet algorithm based on
%         ECGt_QRSd.detectors = 'user:example_worst_ever_ECG_delineator';
%         % Example of how you can add your own ECG delineator. 
%         ECGt_QRSd.detectors = 'user:your_ECG_delineator_func_name';    
%         "your_ECG_delineator_func_name" can be your own delineator.
        ECG_all_wrappers(ii).ECGtaskHandle.delineators = {'wavedet' 'user:example_worst_ever_ECG_delineator'};
        
        % process the task
        ECG_all_wrappers(ii).Run;
        
    end

    % at the end, report problems if happened.
    for ii = 1:lrecnames
        ECG_all_wrappers(ii).ReportErrors;
    end
    
    %% Visual inspection of the detection/delineation
    
    if( bUseDesktop )

        % other task can be performed on the same objects
        for ii = 1:lrecnames

            % last worker is the responsible of the visual correction.
            if( ECG_all_wrappers(ii).this_pid == ECG_all_wrappers(ii).cant_pids)

                % if there are not any previous error.
                if( ECG_all_wrappers(ii).Processed && ~ECG_all_wrappers(ii).Error ) 

                    % this is to use previous saved results as starting point,
                    % if any available
                    cached_filenames = ECG_all_wrappers(ii).GetCahchedFileName( {'Delineation_corrector' 'ECG_delineation'} );

                    % if no previous correction work, try the automatic
                    % detection task
                    
                    % if any, do the correction
                    if( ~isempty(cached_filenames) )

                        % this is to use previous saved results as starting point,
                        % if any available
                        ECG_all_wrappers(ii).ECGtaskHandle = 'Delineation_corrector';
                        
                        % This task is supposed to be supervised, so only one pid is enough.
                        ECG_all_wrappers(ii).this_pid = '1/1';

                        % user provided name to individualize each run
                        ECG_all_wrappers(ii).ECGtaskHandle.user_string = user_str;
                        
                        % to avoid loading cached results and exit, this flag
                        % allows the re-editing of the current state of the
                        % detections.
                        ECG_all_wrappers(ii).cacheResults = false;

                        % maybe in your application you should run this for
                        % all files.
                        ECG_all_wrappers(ii).ECGtaskHandle.payload = load(cached_filenames{1});

                        % process the task
                        ECG_all_wrappers(ii).Run;

                        % restore the original pids configuration
                        ECG_all_wrappers(ii).this_pid = pid_str;

                        % As we changed for "QRS correction" task, we have to enable this
                        % value again in order to avoid performing the following tasks every time.
                        % If you want to recalculate any task, change it to false
                        ECG_all_wrappers(ii).cacheResults = true;
                        
                    end

                end

            end

        end

        % at the end, report problems if happened.
        for ii = 1:lrecnames
            ECG_all_wrappers(ii).ReportErrors;
        end

    end
    
    %% Automatic Heartbeat classification
    
    % other task can be performed on the same objects
    for ii = 1:lrecnames
        
        % this is to use previous cached results as starting point
        cached_filenames = ECG_all_wrappers(ii).GetCahchedFileName('QRS_corrector');
        % if corrected QRS detections are not available, wavedet
        % performs automatic QRS detection.
        
        if( ~isempty(cached_filenames) )
            
            ECG_all_wrappers(ii).ECGtaskHandle = 'ECG_heartbeat_classifier';

            ECG_all_wrappers(ii).ECGtaskHandle.payload = load(cached_filenames{1});
            
            % user provided name to individualize each run
            ECG_all_wrappers(ii).ECGtaskHandle.user_string = user_str;

            % process the task
            ECG_all_wrappers(ii).Run;

        end

    end

    % at the end, report problems if happened.
    for ii = 1:lrecnames
        ECG_all_wrappers(ii).ReportErrors;
    end
    
    %% Visual inspection of the signal
    

    filename = []; % default setting. Let the report function decide.
%     filename = 'container_filename'; % to put everything in one big file.
    
    % other winlengths can be added to the array in order to further
    % explore the recordings, and the algorithm results.
%     winlengths = []; % default setting
    winlengths = [ 7 ]; %seconds

    % go through all files
    for ii = 1:lrecnames

        if( ECG_all_wrappers(ii).this_pid == ECG_all_wrappers(ii).cant_pids)

            % last worker is the responsible of the reporting.
            if( ECG_all_wrappers(ii).this_pid == ECG_all_wrappers(ii).cant_pids)

                try
                    
                    reportECG(ECG_all_wrappers(ii), 'LowDetail', 'full', winlengths, 'pdf', [] );
            %         reportECG(ECG_all_wrappers(ii), 'LowDetail', 'full', winlengths, 'png', filename);
                catch MException

                    report = getReport(MException);

                    fprintf(2, '\n%s\n', report);

                end

            end

        end
        
    end
    
    %% other user-defined tasks ...
        
    if( ~bUseDesktop )
    
        UnInstallECGkit();
        
    end
    
    
    