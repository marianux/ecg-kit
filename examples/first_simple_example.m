%% Very simple example of how to use the ECGkit
% 
% This script exemplifies the use of the ECGkit in a very short multimodal
% cardiovascular recording which includes arterial blood pressure (ABP),
% plethysmographic (PPG) and electrocardiogram signals. The following tasks
% will be performed in this example:  
% 
% * Heartbeat/QRS detection
% * ABP/PPG pulse detection
% * ECG wave delineation
% * Heartbeat classification
% * Report generation
% 
% Each automatic step is followed by a manual verification step in order to
% verify the algorithm's results.
% 
% You can watch a typical run of this script online
%  
%   http://youtu.be/8lJtkGhrqFw?list=PLlD2eDv5CIe9sA2atmnb-DX48FIRG46z7
% 
% 
% See also ECGwrapper, ECGtask, UnInstallECGkit, ECGkit_examples
% 
% *Author*: Mariano Llamedo Soria
% <matlab:web('mailto:llamedom@electron.frba.utn.edu.ar','-browser') (email)> 
% *Version*: 0.1 beta
% *Birthdate*: 11/13/2014
% *Last update*: 15/4/2015
% *Copyright* 2008-2015
% 
function first_simple_example()

    root_path = fileparts(mfilename('fullpath'));
    % get the parent or ECGkit root folder.
    root_path = fileparts(root_path);
    % default folder to look at
    examples_path = [root_path filesep 'recordings' filesep ];

    % the filename of the recording
    example_filename = '208.hea';
%     example_filename = '800.hea';
%     example_filename = 'example_recording.mat';

    % user-defined string to customize the experiment
    my_user = 'my_experiment_name';
    % Set to true, if you want to test the GUI for auditing and correcting
    % the automatic algorithm results.
    bGUICorrection = false;
    
    
    % this is an ECGwrapper object, and is the interface to get samples
    % from the recordings to the processing algorithms.
    ECGw = ECGwrapper('recording_name', [examples_path example_filename]);
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

    if( bGUICorrection )
        % this step if for correcting the automatic detections performed before
        ECGw.ECGtaskHandle = 'QRS_corrector';
        % with the following two lines, you can use as starting point the
        % manual corrected detections (if available) or the automatic
        % detections.
        cached_filenames = ECGw.GetCahchedFileName({'QRS_corrector' 'QRS_detection'});
        ECGw.ECGtaskHandle.payload = load(cached_filenames{1});
        ECGw.Run
    
        % same as before, but for the pulsatile signals.
        ECGw.ECGtaskHandle = 'PPG_ABP_corrector';
        cached_filenames = ECGw.GetCahchedFileName({'PPG_ABP_corrector' 'PPG_ABP_detector'});
        ECGw.ECGtaskHandle.payload = load(cached_filenames{1});    
    end
    
    % here we perform ECG delineation based on the previous QRS detections
    % (the corrected ones if available) with ALL the algorithms present in
    % the toolbox (wavedet only at the moment of writing this)
    ECGw.ECGtaskHandle = 'ECG_delineation';
    cached_filenames = ECGw.GetCahchedFileName({'QRS_corrector' 'QRS_detection'});
    ECGw.ECGtaskHandle.payload = load(cached_filenames{1});
    ECGw.Run    
    
    if( bGUICorrection )
        % in case you would like to correct the delineation marks, you can do
        % it with this task.
        ECGw.ECGtaskHandle = 'ECG_delineation_corrector';
        cached_filenames = ECGw.GetCahchedFileName({'ECG_delineation_corrector' 'ECG_delineation'});
        ECGw.ECGtaskHandle.payload = load(cached_filenames{1});
        ECGw.Run
    end
    
    % finally you can classify the heartbeats according to the AAMI
    % classes, also based in the previous QRS detections (again, the
    % corrected ones if available) 
    ECGw.ECGtaskHandle = 'ECG_heartbeat_classifier';
    ECGw.ECGtaskHandle.mode = 'auto';
    cached_filenames = ECGw.GetCahchedFileName({'QRS_corrector' 'QRS_detection'});
    ECGw.ECGtaskHandle.payload = load(cached_filenames{1});
    ECGw.Run    

    % at the end, generate a report pretty-printed report with the results
    reportECG(ECGw, 'LowDetail', 'full')
    