%% Very simple example of how to use the ECGkit in a cluster
% 
% This script is the same as "first_simple_example" but it is prepared to
% run locally without arguments, or in a cluster environment by
% using "pid_str" argument. The pid_str argument is a char with format
% 'N/M', being N <= M with default value '1/1'. You can partition a big job
% into M pieces in a cluster architecture, by starting M processes with N
% ranging from 1 to M.  
% 
% See also ECGwrapper, ECGtask, UnInstallECGkit, first_simple_example
% 
% *Author*: Mariano Llamedo Soria
% <matlab:web('mailto:llamedom@electron.frba.utn.edu.ar','-browser') (email)> 
% *Version*: 0.1 beta
% *Birthdate*: 11/13/2014
% *Last update*: 18/14/2014
% *Copyright* 2008-2014
% 
function second_simple_example(pid_str)

    if( nargin < 1 || ~ischar(pid_str) )
        % single PID run
        pid_str = '1/1';
    end

    root_path = fileparts(mfilename('fullpath'));
    % get the parent or ECGkit root folder.
    root_path = fileparts(root_path);
    % default folder to look at
    examples_path = [root_path filesep 'recordings' filesep ];

    % the filename of the recording. In this example, be sure to choose a
    % longer ECG recording. Not included with the source code, get one from
    % Physionet and name this variable accordingly.
    example_filename = 'longer_example_recording.mat';

    % user-defined string to customize the experiment
    my_user = 'my_experiment_name';

    % In case of running in a user-assisted fashion.
    bUseDesktop = usejava('desktop');

    if( ~bUseDesktop )
        % For cluster or distributed environment processing.
        InstallECGkit();
    end
    
    
    % this is an ECGwrapper object, and is the interface to get samples
    % from the recordings to the processing algorithms.
    ECGw = ECGwrapper('recording_name', [examples_path example_filename]);
    % this is the to indicate an ECGwrapper object the way it shoud deal
    % with multi-processing
    ECGw.this_pid = pid_str;
    % in this example we are not interested in caching results, since the
    % recording is very short, so we will re-process every call.
    ECGw.cacheResults = false;
    % add a user-defined prefix to the result filenames
    ECGw.user_string = my_user;
    
%     % QRS detection in ALL the signals present in the recording, with ALL
%     % the QRS detection algorithms available.
%     ECGw.ECGtaskHandle = 'QRS_detection';
%     ECGw.Run    
    
    % QRS detection, but only in those signals that seems to be ECG. That
    % mean, that the description of the signal give some clue about ECG. If
    % not sure, process all signals as ECG.
    ECGw.ECGtaskHandle = 'QRS_detection';
    ECGw.ECGtaskHandle.only_ECG_leads = true;
    % just restrict the run to the three algorithms
    ECGw.ECGtaskHandle.detectors = { 'wavedet', 'gqrs', 'wqrs'};
    ECGw.Run    
    
    % Pulse detection in pulsatile signals
    ECGw.ECGtaskHandle = 'PPG_ABP_detector';
    % only in pulsatile signals. This is donde as before, based on the
    % description of the signal, if unsure, process all signals.
    ECGw.ECGtaskHandle.lead_config = 'PPG-ABP-only';
    ECGw.Run

    % here we perform ECG delineation based on the previous QRS detections
    % (the corrected ones if available) with ALL the algorithms present in
    % the toolbox (wavedet only at the moment of writing this)
    ECGw.ECGtaskHandle = 'ECG_delineation';
    cached_filenames = ECGw.GetCahchedFileName({'QRS_corrector' 'QRS_detection'});
    if( ~isempty(cached_filenames) )
        ECGw.ECGtaskHandle.payload = load(cached_filenames{1});
        ECGw.Run    
    end
    
    % finally you can classify the heartbeats according to the AAMI
    % classes, also based in the previous QRS detections (again, the
    % corrected ones if available) 
    ECGw.ECGtaskHandle = 'ECG_heartbeat_classifier';
    ECGw.ECGtaskHandle.mode = 'auto';
    cached_filenames = ECGw.GetCahchedFileName({'QRS_corrector' 'QRS_detection'});
    if( ~isempty(cached_filenames) )
        ECGw.ECGtaskHandle.payload = load(cached_filenames{1});
        ECGw.Run    

        % only the master PID will generate the report
        if( ECGw.this_pid == ECGw.cant_pids )
            % at the end, generate a report pretty-printed report with the results
            reportECG(ECGw, 'LowDetail', 'full')
        end
        
    end
        
    
    if( ~bUseDesktop )
        UnInstallECGkit();
    end
    