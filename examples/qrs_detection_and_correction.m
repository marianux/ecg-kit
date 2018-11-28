function qrs_detection_and_correction()

root_path = fileparts(mfilename('fullpath'));
% get the parent or ecg-kit root folder.
root_path = fileparts(root_path);
% default recordings path
database_path = [root_path filesep 'recordings' filesep ];

output_path = [ database_path 'qrs_test/' ];

rec_names = list_recordings(database_path);

rec_names = cellfun( @(a)(a(1:end-4) ), rec_names, 'UniformOutput', false );

% Arbitrary Impulsive Pseudoperiodic (AIP) configuration
payload_in.trgt_width = 0.06; % seconds
payload_in.trgt_min_pattern_separation = 0.3; % seconds
payload_in.trgt_max_pattern_separation = 2; % seconds
payload_in.max_patterns_found = 2; % # de morfologías o latidos a buscar

% Use previous results flag. Set to false to recalculate always.
bCached = true;


for this_rec = rec_names
    
    % For each recording *aa* found, this script will generate a result
    % file for each step in the processing:
    %
    % aa_QRS_detection.mat
    % aa_AIP_det_arbitrary_function.mat
    % aa_QRS_detections_post_process.mat
    % aa_reviewed_annotations.mat
    %
    % the last file is the *best* QRS detections perfomed, and is the result
    % of 1) automatic detection, 2) post-processing and merge and 3) visual
    % inspection and correction.
    % the first step is checking if this file has already been processed.
    
    rev_anns = [output_path this_rec{1} '_reviewed_annotations.mat' ];

    if(exist(rev_anns, 'file'))

        continue

    else

%% QRS detection
%%%%%%%%%%%%%%%%%

        % Wavedet and GQRS configuration

        fname = [database_path this_rec{1} ];
        
        ECGw = ECGwrapper( 'recording_name', fname);
        
        ECGw.output_path = output_path;
        ECGw.ECGtaskHandle = 'QRS_detection';
        ECGw.ECGtaskHandle.only_ECG_leads = true;
        ECGw.ECGtaskHandle.detectors = { 'wavedet', 'gqrs'};

        ECGw.cacheResults = bCached; 

        ECGw.Run

        cached_filenames = ECGw.Result_files;
        QRS_detectors_struct = load(cached_filenames{1});

        % Arbitrary Impulsive Pseudoperiodic (AIP) detector: Another
        % unpublished detector

        ECGw.ECGtaskHandle = 'arbitrary_function';

        ECGw.ECGtaskHandle.payload = payload_in;

        % Add a user-string to identify the run
        ECGw.user_string = 'AIP_det';

        % add your function pointer
        ECGw.ECGtaskHandle.function_pointer = @aip_detector;
        ECGw.ECGtaskHandle.concate_func_pointer = @aip_detector_concatenate;

        ECGw.cacheResults = bCached; 

        ECGw.Run

        cached_filenames = ECGw.Result_files;
        QRS_aip_struct = load(cached_filenames{1});

%% Post process detections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Merge together detections generated with:
        % 
        %   + Wavedet
        %   + GQRS
        %   + AIP
        % 

        ECGw.ECGtaskHandle = 'QRS_detections_post_process';

        % Mix QRS detections strategy function

        ECGw.ECGtaskHandle.post_proc_func = 'mixartif';

        QRS_detections = concat_QRS_detections( QRS_detectors_struct, QRS_aip_struct);

        ECGw.ECGtaskHandle.payload = QRS_detections;

        ECGw.cacheResults = bCached; 

        ECGw.Run

        cached_filenames = ECGw.Result_files;

        post_proc_struct = load(cached_filenames{1});

%% Detections review and correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Commonly the 

        ECGw.ECGtaskHandle = 'QRS_corrector';

        ECGw.ECGtaskHandle.payload = concat_QRS_detections( post_proc_struct, QRS_detections) ;

        ECGw.ECGtaskHandle.leads = {'II', 'V2', 'V5'};

        ECGw.cacheResults = false; 

        ECGw.Run

    end
end

