%% Use of ecg-kit for the detection-validation-correction of heartbeats in ECG
% 
%  Brief description:
%  ------------------
% 
% Example or recipe for the detection of heartbeats in a dataset or folder
% containing several ECG recordings. Several features of the toolbox were
% described in previous examples. See other scripts available in this
% folder for further details.
% This example or recipe follows these steps:
% 
%   1- Automatic QRS detection
%       a- wavedet
%       b- gqrs
%       c- AIP
%   2- Merge and selection of the best detection.
%   3- Visual inspection/correction of the results.
%   4- Export the validated detections to a file.
% 
% The goal of this example is showing a method for creating high quality
% detections for research or other purposes.
% 
% Author: Mariano Llamedo Soria (llamedom at frba.utn.edu.ar)
% Version: 0.1 beta
% Birthdate: 12/05/2018
% Last update: 12/05/2018
% Copyright 2008-2018
% 

function qrs_detection_and_correction()

root_path = fileparts(mfilename('fullpath'));
% get the parent or ecg-kit root folder.
root_path = fileparts(root_path);
% default recordings path
database_path = [root_path filesep 'recordings' filesep ];

output_path = [ database_path 'qrs_test/' ];

% recordings definition:
% List the recognized formats in *database_path* 
rec_names = list_recordings(database_path);
rec_names = cellfun( @(a)(a(1:end-4) ), rec_names, 'UniformOutput', false );

% % OR ... enumerate the recordings you want to work with
% rec_names = { ...
%               '208' ...
%               '800' ...
%             };
        
        

% Arbitrary Impulsive Pseudoperiodic (AIP) configuration:
% This is an unpublished (at the moment of writing this) detector for
% arbitrary patterns, in this case, we briefly describe the pattern to
% match in terms of:
payload_in.trgt_width = 0.06; % 60 milliseconds width. A narrow and normal QRS complex
% The minimum separation among heartbeats, 300 millisecond allow most of
% the anticipated heartbeats to be detected. MITDB has even more
% anticipated heartbeats, but 300 milliseconds is a decent starting point.
payload_in.trgt_min_pattern_separation = 0.3; % seconds. 
% The longest pause allowed among patterns. By the moment this property is
% unused, but future releases can detect pauses in order to restart the
% pattern search, or relax certain restrictions.
payload_in.trgt_max_pattern_separation = 2; % seconds
% amount of patterns to find, matching the above criteria. Could be related
% to the heartbeats morphologies in QRS detection context.
payload_in.max_patterns_found = 2; % patterns

% Use previous results flag. Set to false to recalculate always.
bCached = true;


for this_rec = rec_names
    
    % For each recording *aa* found, this script will generate a result
    % file for each step in the processing:
    %
    % aa_QRS_detection.mat
    % aa_AIP_det_arbitrary_function.mat
    % aa_QRS_detections_post_process.mat
    % aa_reviewed_annotations.mat         <-- this is the final result 
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
        % *mixartif* use the output of *all* leads and detectors to perform
        % a multi-lead-and-algorithm composition. The result of this
        % algorithm is a new set of detections based on the concatenation
        % of the ""best"" detections found for each 20-seconds window in a
        % recording. So, this algorithm generates *new* QRS detection
        % sieres.
        ECGw.ECGtaskHandle.post_proc_func = 'mixartif';
        
        % Add a user-string to identify the run
        ECGw.user_string = 'mixartif';

        QRS_detections = concat_QRS_detections( QRS_detectors_struct, QRS_aip_struct);

        ECGw.ECGtaskHandle.payload = QRS_detections;

        ECGw.cacheResults = bCached; 

        ECGw.Run

        cached_filenames = ECGw.Result_files;

        mixartif_struct = load(cached_filenames{1});
        
        % Lead selection strategy based on the best (quality metric) *m*.
        % This is a lead selection algorithm, so this algorithm filter the
        % first ranked lead of all the algorithms outputs included in the
        % payload (QRS_detections struct)
        ECGw.ECGtaskHandle.post_proc_func = 'best_m_lead';

        % Add a user-string to identify the run
        ECGw.user_string = 'best_m_lead';
        
        ECGw.ECGtaskHandle.payload = QRS_detections;

        ECGw.cacheResults = bCached; 

        ECGw.Run

        cached_filenames = ECGw.Result_files;

        best_m_struct = load(cached_filenames{1});

%% Detections review and correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Commonly the 

        ECGw.ECGtaskHandle = 'QRS_corrector';

        ECGw.ECGtaskHandle.payload = concat_QRS_detections( mixartif_struct, best_m_struct, QRS_detections) ;

        % For several applications, these leads provide a good perspective
        % of the electrical activity. The corrector will use these leads if
        % available.
        ECGw.ECGtaskHandle.leads = {'II', 'V2', 'V5'};

        ECGw.cacheResults = false; 

        ECGw.Run

    end
end

