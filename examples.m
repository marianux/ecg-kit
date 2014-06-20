function examples()

    %% Start
    
    % This script exemplifies the use of the toolbox in several typical
    % tasks
    
    % inspect ECG files in rootpath\example recordings\ folder
    root_path = fileparts(mfilename('fullpath'));
    examples_path = [root_path filesep 'example recordings' ];
    filenames = dir(examples_path);
    recnames = {filenames(:).name};
    [~,recnames] = cellfun(@(a)(fileparts(a)), recnames, 'UniformOutput', false);
    recnames = unique(recnames);
%     recnames = recnames(3:end);
    recnames = recnames(3);
    lrecnames = length(recnames);

    
    %% QRS automatic detection
    
    % go through all files
    
    ECG_all_wrappers = [];
    
    for ii = 1:lrecnames
        
        rec_filename = [examples_path filesep recnames{ii}];
        
        % task name, 
%         ECGt_QRSd = 'QRS_detection';
        % or create an specific handle to have more control
        ECGt_QRSd = ECGtask_QRS_detection(); 
%         % select an specific algorithm. Default: Run all detectors
%         ECGt_QRSd.detectors = 'wavedet'; % Wavedet algorithm based on
%         ECGt_QRSd.detectors = 'pantom';  % Pan-Tompkins alg.
        ECGt_QRSd.detectors = 'gqrs';    % WFDB gqrs algorithm.

        ECG_w = ECGwrapper( 'recording_name', rec_filename, ...
                           'ECGtaskHandle', ECGt_QRSd);
        
        % process the task
        ECG_w.Run;
        
        % collect objects
        ECG_all_wrappers = [ECG_all_wrappers; ECG_w];
    
    end

    % at the end, report problems if happened.
    for ii = 1:lrecnames
        ECG_all_wrappers(ii).ReportErrors;
    end
    
    %% QRS visual inspection and correction
    
    % other task can be performed on the same objects
    for ii = 1:lrecnames
        
        if( ECG_all_wrappers(ii).Processed && ~ECG_all_wrappers(ii).Error ) 
            
            % to avoid loading previous results and exit
            ECG_all_wrappers(ii).cacheResults = false;
            ECG_all_wrappers(ii).ECGtaskHandle = 'QRS_corrector';
            
            % this is to use previous cached results as starting point
            cached_filenames = ECG_all_wrappers(ii).GetCahchedFileName();
            
            if( isempty(cached_filenames) )
                % this is to use previous result from the automatic QRS
                % detection
                cached_filenames = ECG_all_wrappers(ii).Result_files;
            end
            
            ECG_all_wrappers(ii).Task_payload = load(cached_filenames{1});

            % process the task
            ECG_all_wrappers(ii).Run;

        end
        
    end
    

    % at the end, report problems if happened.
    for ii = 1:lrecnames
        ECG_all_wrappers(ii).ReportErrors;
    end
    
    %% ECG automatic delineation
    
    % other task can be performed on the same objects
    for ii = 1:lrecnames
        
        ECG_all_wrappers(ii).cacheResults = true;
        
        % obtain the filename of the cached results from the previous
        % QRS_corrector task.
        ECG_all_wrappers(ii).ECGtaskHandle = 'QRS_corrector';
        % this is to use previous cached results as starting point
        cached_filenames = ECG_all_wrappers(ii).GetCahchedFileName();
        % if corrected QRS detections are not available, wavedet
        % performs automatic QRS detection.
        if( ~isempty(cached_filenames) )
            ECG_all_wrappers(ii).Task_payload = load(cached_filenames{1});
        end

        % set the delineator task name and run again.
        ECG_all_wrappers(ii).ECGtaskHandle = 'ECG_delineation';

        % process the task
        ECG_all_wrappers(ii).Run;
            
    end

    % at the end, report problems if happened.
    for ii = 1:lrecnames
        ECG_all_wrappers(ii).ReportErrors;
    end
    
    
    %% Visual inspection of the detection/delineation
    
    
    %% Automatic Heartbeat classification
    
    % other task can be performed on the same objects
    for ii = 1:lrecnames
        
        if( ECG_all_wrappers(ii).Processed && ~ECG_all_wrappers(ii).Error ) 

            ECG_all_wrappers(ii).cacheResults = true;

            % obtain the filename of the cached results from the previous
            % QRS_corrector task.
            ECG_all_wrappers(ii).ECGtaskHandle = 'QRS_corrector';
            % this is to use previous cached results as starting point
            cached_filenames = ECG_all_wrappers(ii).GetCahchedFileName();
            % if corrected QRS detections are not available, wavedet
            % performs automatic QRS detection.
            if( ~isempty(cached_filenames) )
                ECG_all_wrappers(ii).Task_payload = load(cached_filenames{1});
            end
            
            ECG_all_wrappers(ii).ECGtaskHandle = 'ECG_heartbeat_classifier';

            % process the task
            ECG_all_wrappers(ii).Run;

        end
        
    end

    % at the end, report problems if happened.
    for ii = 1:lrecnames
        ECG_all_wrappers(ii).ReportErrors;
    end
    
    %% Visual inspection of the signal
    
    % go through all files
    
    filename = []; % default setting. Let the report function decide.
%     filename = 'container_filename'; % to put everything in one big file.
%     winlengths = []; % default setting
    winlengths = [7 15]; %seconds
    
    for ii = 1:lrecnames
        
        rec_filename = [examples_path filesep recnames{ii}];

        ECG_w = ECGwrapper( 'recording_name', rec_filename);
        
        reportECG(ECG_w, 'LowDetail', 'full', winlengths, 'pdf', filename);
%         reportECG(ECG_w, 'LowDetail', 'full', winlengths, 'png', filename);
    
    end
        
    
    %% other user-defined tasks ...
    
    