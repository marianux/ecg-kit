maxNoAnswers = 3;

cKnownFormats = {'MIT' 'ISHNE', 'AHA', 'HES', 'MAT'};
cHeaderFieldNamesRequired = {'freq' 'nsamp' 'nsig' 'gain' 'adczero' };
cAnnotationsFieldNamesRequired = {'time' };
cKnownLabelings = {'AAMI', 'AAMI2'};

typical_lablists =  { ...
                    { 'Normal'; 'Supraventricular'; 'Ventricular'; 'Fusion'; 'Unclass'} ; ...   AAMI
                    { 'Normal'; 'Supraventricular'; 'Ventricular'; 'Unclass'} ; ...             AAMI2
                    };


% maxQRSxIter: Maximum amount of heartbeats in ECG recordings of 2 leads
maxQRSxIter = 5e3;

% Minimum amount of heartbeats to be processed by a PID
min_QRS_perPID = 1500;

overlapping_time = 30; % seconds
cKnownModesOfOperation = {'auto' 'slightly-assisted', 'assisted'};
maxOperMode = length(cKnownModesOfOperation); % maxOperMode modes of operation

% Working sampling rate.
sampling_rate = 360; % Hz


%gzip compression ratio
typical_compression_ratio = 1.5; %times

% target maximum report file size
max_report_filesize = 4 * 1024^2; % bytes

%Java user interface is started. Not started in clusters for example.
bHaveUserInterface = usejava('desktop');

if( ~bHaveUserInterface )
    %Cluster settings. Ignore them if running in a PC style computer.
    
    %Limits of the cluster to have multiple process accesing I/O
    Loops2io = 100;
    MaxNodesReading = 15;
    MaxNodesWriting = 15;
    
end

% Veces que intentar� ver si un archivo existe en filesystems distribuidos,
% como los del cluster.
%Delayed mounting of filesystems in CLusters. Check many times befor
%asserting.
if( bHaveUserInterface )
   Retries2findfiles = 1;
else
   Retries2findfiles = 3;
end

%
Time2WaitPIDs = 15 * 60; % seconds.

%path related constants.
mylocation_path = [fileparts(mfilename('fullpath')) filesep ];

%Check compilation of source MEX files
common_path = [mylocation_path 'common' filesep];
source_files = dir([ common_path '*.c'] );
lsource_files = length(source_files);
for ii = 1:lsource_files
    [~, source_file_name] = fileparts( source_files(ii).name);
    mex_file = dir([common_path  source_file_name '.' mexext ]);
    if( isempty(mex_file) || mex_file.datenum <= source_files(ii).datenum  )
        eval(['mex -outdir ''' common_path ''' ''' [common_path source_files(ii).name] '''']);
    end
end

%the trained linear classifier.
% load('ldc_classifier.mat');
load('global_classifier.mat');
