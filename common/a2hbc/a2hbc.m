function [ Labels, ConfusionMatrix, LabelList ] = a2hbc(varargin)

% Argentino-Aragonés heartbeat classifier (a2hbc) for Matlab
% ----------------------------------------------------------
% 
% Description:
% 
% Performs ECG heartbeat classification in the input ECG signal in one of the following modes:
% 1. Automatic or Unassisted (default).
% 2. Slightly assisted. Ask minimum user assistance for cluster decision.
% 3. Fully assisted. Ask a decision for each cluster detected.
% 
% The classification is performed among the four AAMI classes: normal
% (N), supraventricular (S) and Ventricular (V) heartbeats. The 'Unknown'
% (U) class should in case of doubt or imposibility to label certainly a
% heartbeat.
% 
% 
% Arguments: (specified as a2hbc('arg_name1', arg_val1, ... , 'arg_nameN', arg_valN) )
% 
%     a. ECG specification:
% 
%         a.1 Recording filename where ECG is stored in a valid format.
% 
%           + recording_name: ECG recording to be classified.
%           + recording_format : Valid ECG format. (MIT, ISHNE, AHA, HES, MAT)
% 
%         a.2 ECG signal/s ADC samples, header and QRS locations. Useful if 
%             you plan to hook this software in your code. In this case you
%             are responsible of coding the I/O to the ECG recording.
% 
%           + ECG: ECG signal matrix. Columns are leads.
% 
%           + ECG_header: Description of the ECG typically available in the
%                         header. Structure with fields:
%                           -freq: Sampling rate in Hz.
%                           -nsig: Number of ECG leads.
%                           -nsamp: Number of ECG samples.
%                           -adczero: ADC offset (typically 2^(adc_res-1) when using unsigned integers).
%                           -adcgain: ADC gain in units/adc_sample (typically uV/adc_unit).
% 
%           + QRS_locations: QRS locations of the ECG. Typically available in the
%                            "time" field of the annotations in MIT format. 
% 
%     b. Operating modes
% 
%         + op_mode: Operating mode 'auto' (default), 'slightly-assisted'
%                   or 'assisted'. Each described above.
% 
%         + NumOfClusters: Number of cluster to find in the data.
%       
%         + ClusteringRepetitions: Repetitions in order to increase the
%               resolution of the clustering. The higher the repetitions,
%               the higher the amount of clusters to find. See the
%               documentations for details. Default: 1 (no repetitions).
%       
%         + ClusterPresence: Threshold for the qualified majority described
%               in the articles. Default: 75%
%       
%     c. Modifiers
% 
%         + InteractiveMode: Boolean for interacting with A2HBC with a GUI.
%                            Default: false.
% 
%         + tmp_path: Path to store temp files. Default: $A2HBC_PATH$\tmp
% 
%       c.1 Multiprocess modifiers:
% 
%         + cant_pids: Number of processes to divide the work. Each process
%                      will work in a 1/cant_pids part of the ECG.
%         + this_pid: Identifies each process with an integer between
%                     1 - cant_pids.
% 
% 
%       c.2 Other modifiers:
% 
%         + CacheData: Boolean for allowing A2HBC to cache data in order to
%                      speedup future execution of the same recording.
%                      Default: true.
% 
%         + SimulateExpert: Boolean for using expert annotations (field 'anntyp' 
%                           of annotation structure) as gold standard and
%                           evaluate the classification performance of
%                           A2HBC. Default: false
% 
%         + Repetitions: Repetitions of the classification process. Useful
%                        when interested in obtaining a center and
%                        dispersion of the performance distribution.
% 
% Output:
%   + Labels: Classification Labels for each QRS location provided.
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
% Note: If you found this software useful, you may  be interested in reading these articles:
% 
% References:
% 
% [0] Llamedo Soria, Mariano. Automatic Processing and Classification of
%     Electrocardiogram for the Detection of Risk Indexes. PhD thesis.
% 
% [1] Llamedo, M. Martï¿½nez, J. Heartbeat Classification Using Feature
%     Selection driven by Database Generalization Criteria. IEEE
%     Transactions on Biomedical Engineering, 2011, 58, 616-625
% [2] Llamedo, M. & Martï¿½nez, J.P. Cross-Database Evaluation of a Multilead  
%     Heartbeat Classifier. IEEE Transactions on Information Technology in 
%     Biomedicine, 2012 expected. Under review, with minor revision.
% [3] Llamedo, M. and Martï¿½nez, J. P. Multidatabase Evaluation of an
%     Automatic Algorithm for ECG Heartbeat Classification with Adjustable
%     Patient-Adaptation Capability. Under review. 
% 
% @ARTICLE{Llamedo11, author = {Llamedo, M. and Martï¿½nez, J.P.}, title = {Heartbeat Classification Using Feature Selection driven by Database Generalization Criteria}, journal = {IEEE Transactions on Biomedical Engineering}, year = {2011}, volume = {58}, pages = {616-625} }
% @ARTICLE{Llamedo11Multilead, author = {Llamedo, M. and Martï¿½nez, J.P.}, title = {Cross-Database Evaluation of a Multilead Heartbeat Classifier},   journal = {IEEE Transactions on Information Technology in Biomedicine},  year = {2012},  volume = {accepted} }
% @ARTICLE{}
% 
% Limits and Known bugs:
%   Probably a lot :( ... but dont panic! send me feedback if you need
%   help.
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate  : 16/8/2011
% Last update: 9/2/2012

%% Prepare paths, check some context, and start working ...

%path related constants.
% mylocation_path = [fileparts(mfilename('fullpath')) filesep ];
% 
% default_paths = { ...
%                     [ mylocation_path 'scripts' filesep ';' ]; ...
%                 };
%             
% default_paths = char(default_paths)';
% default_paths = (default_paths(:))';
% addpath(default_paths);

% Do path installation during the installation of the ECG kit.
default_paths = '';

%to keep all windows and figures clean and tidy when closing a2hbc
intial_state_open_handles = findall(0);

dbstop if caught error 

% uncomment for debugging.
% bDebug = false;
bDebug = true;

db_status = dbstatus();
if( length(db_status) > 1 && strcmpi(db_status(end).cond, 'caught error')  )
    % already stopping on error caught
    bRestoreErrorStatus = false;
else
    if( bDebug )
        % debug errors inside A2HBC, clear on exit
        dbstop if caught error
        bRestoreErrorStatus = true;
    end
end

Cleanup_hdl = onCleanup(@()DoHouseKeeping(default_paths, intial_state_open_handles, bRestoreErrorStatus));

%% Thats all folks !

[ Labels, ConfusionMatrix, LabelList ] = a2hbc_main(varargin{:});


