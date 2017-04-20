%% Plots and interact with ECG signal
% This function plot ECG signals, eventually with annotation marks such as
% QRS complex locations, or even P, Q, R, S and T wave locations and
% boundaries. 
%   
% Some of the relevant features:
% 
%  + User can interact using mouse shortcuts with several aspects of the
%    visualization, such as zoom, pan and measurements.
%  + Information of the multilead wave boundaries can be added to the ECG,
%    for example the delineation obtained with wavedet.
%  + It can "pretty" present the ECG charts for printing to pdf
%    documents
%  + It can be easily added to your project for debug or result
%    presentation through its versatile interface.
% 
% The mouse interaction was adapted from the Dragzoom function, by Evgeny
% Pr, this can be found in: http://www.mathworks.com/matlabcentral/fileexchange/29276-dragzoom-drag-and-zoom-tool
% 
% Prototype
% 
%              function ECG_hdl = plot_ecg_strip( ECG, varargin )
% 
% Arguments:
%     
%     +ECG: [numeric | char | ECGwrapper] REQUIRED. 
%           Signal matrix of dimension [nsamp nsig] where:
% 
%             - nsamp: time length in samples
%             - nsig: number of ECG leads or number of signals.
% 
%           Recording filename or ECGwrapper object of the ECG recording
%           are also accepted.
% 
%     +ECG_header: [struct] OPTIONAL. 
% 
%             Description of the ECG typically available in the
%             header. Structure with fields:
% 
%               -freq: Sampling rate in Hz. (1)
% 
%               -nsig: Number of ECG leads. (size(ECG,2))
% 
%               -nsamp: Number of ECG samples. (size(ECG,1))
% 
%               -adczero: ADC offset (e.g. 2^(adc_res-1) when
%               using unsigned integers). ( repmat(0, ECG_header.nsig , 1) ) 
% 
%               -adcgain: ADC gain in units/adc_sample
%               (typically uV/adc_unit). ( repmat(1, ECG_header.nsig , 1) )
% 
%               -units: Measurement units. ( repmat('uV', ECG_header.nsig , 1) )
% 
%               -desc: Signal description. ( num2str(colvec(1:ECG_header.nsig)) )
% 
%     +Start_time: [numeric] OPTIONAL. Start time in seconds.
%     
%     +End_time: [numeric] OPTIONAL. Start time in seconds.
%     
%     +QRS_locations: [numeric] OPTIONAL. Default values enclosed in ()
%                 Synchronization sample. In ECG context, this values are
%                 the QRS fiducial point. (empty)
% 
%     +QRS_start_index: [numeric] OPTIONAL. Default values enclosed in ()
%                                 Start at the i-th QRS_start_index heartbeat
%                                 in QRS_locations, or QRS_locations(QRS_start_index). (empty) 
%     
%     +QRS_complexes: [numeric] OPTIONAL. Default values enclosed in () 
%                                     Display the amount of QRS_start_complexes heartbeats from the QRS_start_index. (empty)
%     
%     +Lead_offset: [numeric] OPTIONAL. Default values enclosed in () 
%                             A DC value [nsig 1] to be added to each lead. (0)
%     
%     +Lead_gain: [numeric] OPTIONAL. Default values enclosed in () 
%                           A value [nsig 1] to be multiplied by each lead. (1)
%     
%     +ECG_delineation_single_lead: 
%                         [struct array size nsig] OPTIONAL. Default values enclosed in ()
%                         Annotation struct of size [nsig 1] with fields:
%                         Pon P Poff QRSon qrs QRSoff Ton T Toff. Each
%                         field of size [1 nhb], being nhb the amount of
%                         heartbeats. (empty)   
% 
%     +ECG_delineation_multilead: [struct] OPTIONAL. Default values enclosed in ()
%                         Annotation struct with the same fields of and
%                         characteristics of ECG_delineation_single_lead. (empty)   
% 
%     +Title: [string] OPTIONAL. Default values enclosed in ()
%                      Description title. (recname - time interval)
%             
%     +DetailLevel: [string] OPTIONAL. Default values enclosed in ()
%                            The details included in the ECG plot depends
%                            on the zoom level and the data provided.
%                            Possible values: 'all', 'single-lead',
%                            'multilead' and 'none'. (none)    
% 
%     +PrettyPrint: [bool] OPTIONAL. Default values enclosed in ()
%                          Prepare the plot for printing as a PDF. (false)
% 
%     +Figure_handle: [axes handle] OPTIONAL. Choose the figure to display the
%                                 plot. (handle produced by figure)
% 
% Limits and Known bugs:
%   Probably a lot :( ... but dont panic! send me feedback if you need help.
% 
% Example:
% 
%   plot_ecg_strip([ ECGkitrootpath '\recordings\example_recording.mat'])
% 
%   plot_ecg_strip(ECG)
% 
%   plot_ecg_strip(ECG, 'ECG_header', heasig, 'ECG_delineation_single_lead', positions_single_lead);
% 
%   plot_ecg_strip(ECG, 'ECG_header', heasig, 'ECG_delineation_single_lead', positions_single_lead, 'Start_time', 10*60, 'End_time', 20*60);
% 
% 
% Mouse actions:
% 
%   Normal mode:
%       single-click and holding LB : Activation Drag mode
%       single-click and holding RB : Activation rubber band for region zooming
%       single-click MB             : Activation measuring rubber band mode
%       scroll wheel MB             : Activation Zoom mode
%       double-click LB, RB, MB     : Reset to Original View
% 
%   Magnifier mode ('m' key):
%       single-click LB             : Not Used
%       single-click RB             : Not Used
%       single-click MB             : Reset Magnifier to Original View
%       scroll MB                   : Change Magnifier Zoom
%       double-click LB             : Increase Magnifier Size
%       double-click RB             : Decrease Magnifier Size
% 
%   Hotkeys:
%       'h'                         : Show help
%       '+'                         : Zoom plus
%       '-'                         : Zoom minus
%       'd'                         : Toggle the detail level of the annotations
%       'a'                         : Toggle the annotations graph mode ... 
%       '0'                         : Set default axes (reset to original view)
%       'c'                         : On/Off pointer in crosshair mode
%       'g'                         : If pressed and holding, change lead gain with scroll
%       'o'                         : If pressed and holding, change lead offset with scroll
%       'x'                         : If pressed and holding, zoom and drag works only for X axis
%       'y'                         : If pressed and holding, zoom and drag works only for Y axis
%       'm'                         : If pressed and holding, Magnifier mode on
%       'p'                         : On/Off paper mode
%       'r'                         : Export format (PDF/PNG)
%       's'                         : Export current view
% 
% 
% See also plot_ecg_strip
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Last update: 19/11/2014
% Birthdate  : 15/8/2012
% Copyright 2008-2015
function returned_handles = plot_ecg_strip( ECG, varargin )

returned_handles = [];

%% Constants

% Avoid IO greater than MaxIOread
MaxIOread = 10; %megabytes


Default_time_window = 10; % seconds

%multilead
color_Pwave_global =        [210 255 189]/255;
color_QRScomplex_global =   [255 200 255]/255;
color_Twave_global =        [255 240 170]/255;

%single-lead
PwaveColor =   [210 255 189]/255;
QRScplxColor = [255 200 255]/255;
TwaveColor =   [255 240 170]/255;

cBeatLabels          = { 'N'             'S'             'V'            'F'             'U'         };
cBeatLabelsColorCode = { [0 0 255]/255   [0 255 0]/255   [255 0 0]/255  [0 255 255]/255 [0 0 0]/255 };
cAnnotationsFieldNamesRequired = {'time' 'anntyp'};
cAnnotationFieldNames = { 'Pon' 'P' 'Poff' 'QRSon' 'qrs' 'QRSoff' 'Ton' 'T' 'Tprima' 'Toff' };
cDetailLevels = {'all', 'single-lead', 'multilead', 'none' };
cUnitsVoltages = {'NV' 'UV' 'MV' 'V' 'UV' 'VOLT'};
cKnownReportFormats = {'pdf', 'png'};

% k for reports
target_res = 300; % samples/inch
target_res = target_res / 2.54; % samples/cm

min_cant_samp_seconds = 0.1; % seconds

%% argument definition
p = inputParser;   % Create instance of inputParser class.
addRequired(p,  'ECG', @(x)(~isempty(x) && (isnumeric(x) || ischar(x) || isa(x, 'ECGwrapper') ) ) );
p.addParamValue('ECG_header', [], @(x)(isstruct(x)) );
p.addParamValue('QRS_locations', [], @(x)( isnumeric(x) && all(x > 0) || iscell(x) || isa(x, 'ECGwrapper') ));
p.addParamValue('Start_time', [], @(x)( isnumeric(x) && x >= 0 ) );
p.addParamValue('End_time', [], @(x)( isnumeric(x) && x > 0 ) );
p.addParamValue('QRS_start_index', [], @(x)( isnumeric(x) && x > 0 ) );
p.addParamValue('QRS_complexes', [], @(x)( isnumeric(x) && all(x > 0) ) );
p.addParamValue('Figure_handle', [], @(x)(ishandle(x)) );
p.addParamValue('Lead_offset', [], @(x)( isnumeric(x) ) );
p.addParamValue('Lead_gain', [], @(x)( isnumeric(x) ) );
p.addParamValue('Heartbeat_classification', [], @(x)(isstruct(x) || isa(x, 'ECGwrapper')) );
p.addParamValue('ECG_delineation_single_lead', [], @(x)(isstruct(x) || isa(x, 'ECGwrapper')) );
p.addParamValue('ECG_delineation_multilead', [], @(x)(isstruct(x) || isa(x, 'ECGwrapper')) );
p.addParamValue('Title', [], @(x)(ischar(x)) );
p.addParamValue('LinkedHandle', [], @(x)(ishandle(x)) );
p.addParamValue('DetailLevel', 'dont_care', @(x)( ischar(x) && any(strcmpi(x,cDetailLevels))  ) );
p.addParamValue('PrettyPrint', false, @(x)(islogical(x)) );
p.addParamValue('OnlyECG', false, @(x)(islogical(x)) );
p.addParamValue('FilterECG', false, @(x)(islogical(x)) );
p.addParamValue('ReportFilename', [], @(x)( isempty(x) ) );

try
    p.parse( ECG, varargin{:} );
catch MyError
    fprintf(2, disp_option_enumeration('Incorrect argval/argvalue, use:', p.Parameters ));
    fprintf(2, '\nOr run ''doc plot_ecg_strip'' for details.\n\n' );
    rethrow(MyError);    
end

ECG = p.Results.ECG;
heasig = p.Results.ECG_header;
QRS_locations = p.Results.QRS_locations;
start_time = p.Results.Start_time;
end_time = p.Results.End_time;
QRS_start_idx = p.Results.QRS_start_index;
cant_qrs = p.Results.QRS_complexes;
fig_hdl = p.Results.Figure_handle;
lead_offset = p.Results.Lead_offset;
gains = p.Results.Lead_gain;
hb_labels = p.Results.Heartbeat_classification;
annotations = p.Results.ECG_delineation_single_lead; %multilead
global_annotations = p.Results.ECG_delineation_multilead; %single-lead
strTitle = p.Results.Title;
linked_hdl = p.Results.LinkedHandle;
strDetail = p.Results.DetailLevel;
bPrettyPrint = p.Results.PrettyPrint;
bOnlyECG = p.Results.OnlyECG;
bFilterECG = p.Results.FilterECG;
report_filename = p.Results.ReportFilename;

clear p

aux_val = dbstack;
bCalledFromWorkspace = length(aux_val) == 1;
cLinespecs = [];
this_path = fileparts(mfilename('fullpath'));
cLinespecs = load([ this_path filesep 'clinespecs.mat']);
cLinespecsNone = cLinespecs.cLinespecsNone;
cLinespecs = cLinespecs.cLinespecs;

ECG_w = [];
bWrapper_provided = false;

if( ischar(ECG) )
    
    if( exist(ECG, 'file') )
        ECG_w = ECGwrapper( 'recording_name', ECG);
        bWrapper_provided = true;

        if( isempty(heasig) )
            heasig = ECG_w.ECG_header;
        end

        ECG = [];  
    else
        error('plot_ecg_strip:FileNotFound', 'File %s not found\n', ECG )
    end
    
elseif( isa(ECG, 'ECGwrapper') )
    % parse ECGwrapper object
    ECG_w = ECG;
    bWrapper_provided = true;
    
    if( isempty(heasig) )
        heasig = ECG_w.ECG_header;
    end
    
    ECG = [];
    
end

QRS_locations_wrapper = [];

if( isnumeric(QRS_locations) && ~isempty(QRS_locations) )
    QRS_locations_names = {'user_provided'};
    QRS_locations = {QRS_locations};
elseif( iscell(QRS_locations) && ~all(cellfun(@(a)(isempty(a)), QRS_locations))  )
    cant_QRS_loc_provided = length(QRS_locations);
    QRS_locations_names = arrayfun(@(a)(sprintf('user_provided %d',a)), colvec(1:cant_QRS_loc_provided), 'UniformOutput', false);
elseif( isa(QRS_locations, 'ECGwrapper') )
    QRS_locations_wrapper = QRS_locations;
    QRS_locations_names = {};
    QRS_locations = {};
elseif( bWrapper_provided )
    QRS_locations_wrapper = ECG_w;
    QRS_locations_names = {};
    QRS_locations = {};
end

if( ~isempty(QRS_locations_wrapper) )
    % parse ECGwrapper object
    % get the annotations from the wrapper.

    task_names = {'QRS_corrector' 'QRS_detection'; 'PPG_ABP_corrector' 'PPG_ABP_detector'};

    for kk = 1:size(task_names,1)
        
        cached_filenames = QRS_locations_wrapper.GetCahchedFileName( task_names(kk,:) );
        
        if( ~isempty(cached_filenames) )
        
            aux_annotations = load(cached_filenames{1});
            fnames = fieldnames(aux_annotations);
            aux_idx = find(cell2mat( cellfun(@(a)(~isempty(strfind(a, 'corrected_'))), fnames, 'UniformOutput', false)));
            if( isempty(aux_idx) )
                % no corrected annotations
                if( isfield(aux_annotations, 'series_quality') )
                    QRS_locations_found = length(aux_annotations.series_quality.ratios);
                    [~, aux_idx] = sort(aux_annotations.series_quality.ratios, 'descend');
                    aux_idx2 = 1:min(QRS_locations_found, 3);
                    QRS_locations_names = [QRS_locations_names aux_annotations.series_quality.AnnNames{aux_idx(aux_idx2),1} ];
                    for kk = rowvec(aux_idx(aux_idx2))
                        QRS_locations = [QRS_locations {aux_annotations.(aux_annotations.series_quality.AnnNames{kk,1}).(aux_annotations.series_quality.AnnNames{kk,2})}];
                    end
                end
            else
                aux_val = length(aux_idx);
                for kk = 1:aux_val
                    QRS_locations_names = [QRS_locations_names fnames(aux_idx(kk)) ];
                    QRS_locations = [QRS_locations {aux_annotations.(fnames{aux_idx(kk)}).time}];
                end
            end

        end
        
    end
   
    % check for other annotations, like manual
    aux_annotations = QRS_locations_wrapper.ECG_annotations;
    
    if( ~isempty(aux_annotations) && isfield(aux_annotations, 'time') && ~isempty(aux_annotations.time) )
        QRS_locations_names = [QRS_locations_names {'included'} ];
        QRS_locations = [QRS_locations {aux_annotations.time}];
    end
    
end

if( isa(annotations, 'ECGwrapper') )
    annotations_wrapper = annotations;
elseif( bWrapper_provided )
    annotations_wrapper = ECG_w;
else
    annotations_wrapper = [];
end

if( ~isempty(annotations_wrapper) )
    
    % get the ECG delineation from the wrapper.    
    annotations_wrapper.ECGtaskHandle = 'ECG_delineation';
    
    cached_filenames = annotations_wrapper.GetCahchedFileName();
    if( isempty(cached_filenames) )
        if( ~isempty(annotations) && bCalledFromWorkspace )
            cprintf('[1,0.5,0]', 'Delineation not found for this wrapper object.\n');
        end
        annotations = [];
    else
        annotations = load(cached_filenames{1});
        
        if( isfield(annotations, 'wavedet') )

            if( isfield(annotations.wavedet, 'multilead') )
                global_annotations = annotations.wavedet.multilead;
            end
            
            aux_fields = fieldnames(annotations.wavedet);
            [~, aux_idx ] = intersect( aux_fields, cellstr(strtrim(heasig.desc)) );
            
            aux_struct = [];
            for ii = rowvec(aux_idx)
                aux_struct = [aux_struct annotations.wavedet.(aux_fields{ii}) ];
            end
            
            annotations = aux_struct;
            
        end
    end    
end

% todo: handle multiple delineations structs
if( ~isempty(ECG_w) && ~isempty(ECG_w.ECG_delineation) )
    annotations = ECG_w.ECG_delineation;
end

if( isa(global_annotations, 'ECGwrapper') )
    global_annotations_wrapper = global_annotations;
elseif( bWrapper_provided )
    global_annotations_wrapper = ECG_w;
else
    global_annotations_wrapper = [];
end

if( ~isempty(global_annotations_wrapper) )
    
    % get the ECG delineation from the wrapper.    
    global_annotations_wrapper.ECGtaskHandle = 'ECG_delineation';
    
    cached_filenames = global_annotations_wrapper.GetCahchedFileName();
    if( isempty(cached_filenames) )
        if( ~isempty(global_annotations) && bCalledFromWorkspace  )
            cprintf('[1,0.5,0]', 'Multilead/global Delineation not found for this wrapper object.\n');
        end
        global_annotations = [];
    else
        global_annotations = load(cached_filenames{1});
        
        if( isfield(global_annotations, 'wavedet') )

            if( isfield(global_annotations.wavedet, 'multilead') )
                global_annotations = global_annotations.wavedet.multilead;
            else
                if( bCalledFromWorkspace )
                    cprintf('[1,0.5,0]', 'Multilead/global Delineation not found for this wrapper object.\n');
                end
                global_annotations = [];
            end
        else
            if( bCalledFromWorkspace )
                cprintf('[1,0.5,0]', 'Multilead/global Delineation not found for this wrapper object.\n');
            end
            global_annotations = [];
        end
    end    
end

if( isa(hb_labels, 'ECGwrapper') )
    hb_labels_wrapper = hb_labels;
elseif( bWrapper_provided )
    hb_labels_wrapper = ECG_w;
else
    
    if( isempty(hb_labels) )
        hb_labels_idx = [];
    else
        if( all(isfield(hb_labels, cAnnotationsFieldNamesRequired )) )
            hb_labels_idx = renumlab(hb_labels.anntyp, char(cBeatLabels) );
        else
            hb_labels = [];
            hb_labels_idx = [];
            if( bCalledFromWorkspace )
                cprintf('[1,0.5,0]', disp_option_enumeration('Missing fields within structure:', cAnnotationsFieldNamesRequired ) );
            end
        end
    end
    
    hb_labels_wrapper = [];
end

if( ~isempty(hb_labels_wrapper) )
    
    % get the heartbeat classification from the wrapper.    
    hb_labels_wrapper.ECGtaskHandle = 'ECG_heartbeat_classifier';
    
    cached_filenames = hb_labels_wrapper.GetCahchedFileName();
    if( isempty(cached_filenames) )
        if( ~isempty(hb_labels) && bCalledFromWorkspace )
            cprintf('[1,0.5,0]', 'Heartbeat classification not found for this wrapper object.\n\n');
        end
        hb_labels = [];
        hb_labels_idx = [];
    else
        
        hb_aux = load(cached_filenames{1});

        QRS_locations_names = [QRS_locations_names {'hb_classifier'} ];
        QRS_locations = [QRS_locations {hb_aux.time}];
        hb_labels_idx = renumlab(hb_aux.anntyp, char(cBeatLabels) );
        
%         % match heartbeat classification with global_annotations.qrs time
%         % reference.
% 
%         aux_time = sort(union(hb_aux.time, global_annotations.qrs));
%         aux_labels.time = aux_time;
%         aux_val2 = nan(size(aux_time));
%         aux_labels.anntyp = repmat('U',size(aux_time));
%         [~,~, aux_idx2] = intersect(hb_aux.time, aux_time);
%         aux_labels.anntyp(aux_idx2) = hb_aux.anntyp;
%         hb_labels = aux_labels;
%         
%         [~,~, aux_idx2] = intersect(global_annotations.qrs, aux_time);
%         
%         for fn = rowvec(fieldnames(global_annotations))
%             aux_val3 = aux_val2;
%             aux_val3(aux_idx2) = global_annotations.(fn{1});
%             aux_anns.(fn{1}) = aux_val3;
%         end
%         aux_anns.qrs = aux_time;
% 
%         hb_labels_idx = renumlab(hb_labels.anntyp, char(cBeatLabels) );
        
    end    
    
end

if( ( (isempty(QRS_locations) || (iscell(QRS_locations) && all(cellfun(@(a)(isempty(a)), QRS_locations)))) && isempty(annotations) && isempty(global_annotations) ) || ...
    ( isnumeric(QRS_locations) && (QRS_locations(end) - QRS_locations(1)) < (Default_time_window * heasig.freq) || (iscell(QRS_locations) && all(cellfun(@(a)((a(end) - a(1)) < (Default_time_window * heasig.freq)), QRS_locations)))) ...
  )
    %no annotations at all
    QRS_start_idx = [];
    cant_qrs = [];
    
    if( isempty(start_time) )
        start_time = 0;
    end

    if( isempty(end_time) )
        end_time = realmax;
    end
    
else
    % there are annotations
    
    if( isempty(start_time) && isempty(end_time) )
        
        if( isempty(QRS_start_idx) )
            QRS_start_idx = 1;
        end

        if( isempty(cant_qrs) )
            cant_qrs = min([ 10; colvec(cellfun(@(a)(length(a)), QRS_locations)) ]);
        end

    else

        if( isempty(start_time) )
            if( isempty(end_time) )
                start_time = 0;
            else
                start_time = end_time - Default_time_window;
            end
        end

        if( isempty(end_time) )
            end_time = start_time + Default_time_window;
        end
        
    end
end

if( isempty(heasig) && bCalledFromWorkspace )
%     warning('plot_ecg_strip:UnknownSamplingRate', 'Assuming sampling rate of 1000 Hz.')
    if( bCalledFromWorkspace )
        cprintf('[1,0.5,0]', 'Assuming sampling rate of 1000 Hz.\n')
    end
    heasig.freq = 1000;
end

if( isempty(fig_hdl) )
    fig_hdl = figure;
    maximize(fig_hdl)
else
    % preserve visible property, specially for report generation, when it
    % is not important to show the figure.
    visible_prev = get(fig_hdl, 'Visible');
    clf('reset')
    set(fig_hdl, 'Visible', visible_prev);
    
end
set(fig_hdl, 'ToolBar','none');

axes_hdl = gca;


if( ~isempty(linked_hdl) )
    
    if( isprop(linked_hdl, 'UserData') )
        ud = get(linked_hdl, 'UserData');
        if( isfield(ud, 'linked_hdl') && fig_hdl ~= linked_hdl )
            % set the recyprocal link
            ud.linked_hdl = fig_hdl;
            set(linked_hdl, 'UserData', ud);
            clear ud
        else
            % not a proper link
            linked_hdl = [];
        end
    else
        % not a proper link
        linked_hdl = [];
    end
    
end

if( isfield(heasig, 'recname') )
    set(fig_hdl, 'Name', heasig.recname);
else
    set(fig_hdl, 'Name', 'Unknown ECG recording');
end

% maximize(fig_hdl);

if( isempty(ECG) )
    % in case not all signals are ECG.
    if( bOnlyECG )
        [ECG_signals_idx, heasig] = get_ECG_idx_from_header(heasig);
        if( isempty(ECG_signals_idx) && bCalledFromWorkspace  )
            cprintf('[1,0.5,0]', disp_option_enumeration('Could not find any ECG signal, check the lead description of the recording:', cellstr(heasig.desc) ) );
            fprintf(1, '\n');
        end
        
        cant_leads = length(ECG_signals_idx);
    else
        ECG_signals_idx = get_ECG_idx_from_header(heasig);
    end
    cant_samp = heasig.nsamp;
    if(cant_samp == 0)
        cant_samp = realmax;
    end
    
    % only for big recordings, can it be readed all together ?
    aux_MaxIOread = (MaxIOread * 1024^2) / heasig.nsig / 2;
    if( cant_samp > aux_MaxIOread && isempty(start_time) && isempty(end_time) )
        cant_samp = aux_MaxIOread;
        if( bCalledFromWorkspace )
            cprintf('[1,0.5,0]','Recording too big, reading only %d samples. Use "Start_time" and "End_time" arguments to select other parts.\n', cant_samp);
        end
    end
    
    cant_leads = heasig.nsig;
else

    [cant_samp cant_leads ] = size(ECG);
    % assume more data than channels
    if( cant_leads > cant_samp )
        ECG = ECG';
        [cant_samp cant_leads] = size(ECG);
    end
    
    if( bOnlyECG )
        % retain only ECG signals
        [ECG_signals_idx, heasig] = get_ECG_signals_idx_from_header(heasig);
        if( isempty(ECG_signals_idx)  && bCalledFromWorkspace )
            cprintf('[1,0.5,0]', disp_option_enumeration('Could not find any ECG signal, check the lead description of the recording:', cellstr(heasig.desc) ) );
            fprintf(1, '\n');
        end
        cant_leads = length(ECG_signals_idx);
        ECG = ECG(:,ECG_signals_idx);
    else
        cant_leads = heasig.nsig;
        ECG_signals_idx = 1:cant_leads;
    end
    
end

ECG_signals_idx = rowvec(ECG_signals_idx);

if( ~isfield( heasig, 'nsig' ) )
    heasig.nsig = cant_leads;
end

if( ~isfield( heasig, 'nsamp' ) )
    heasig.nsamp = cant_samp;
end

if( ~isempty(global_annotations) )
    QRS_locations_names = [QRS_locations_names {'global'} ];
    QRS_locations = [QRS_locations {global_annotations.qrs}];
end

if( isempty(QRS_locations) )
    QRS_locations_names = {[]};
    QRS_locations = {[]};
end

cant_QRS_locations = cellfun(@(a)(length(a)), QRS_locations);

if( ~isfield( heasig, 'adczero' ) )
%     warning('plot_ecg_strip:UnknownADCzero', 'Assuming 0 ADC zero.')
    if( bCalledFromWorkspace )
        cprintf('[1,0.5,0]','Assuming 0 ADC zero.\n')
    end
    heasig.adczero = zeros(cant_leads,1);
end

if( ~isfield( heasig, 'gain' ) )
%     warning('plot_ecg_strip:UnknownADCgain', 'Assuming gain 1 uV/samp.')
    if( bCalledFromWorkspace )
        cprintf('[1,0.5,0]','Assuming gain 1 uV/samp.\n')
    end
    heasig.gain = repmat(1, cant_leads,1);
end

% voltages to microvolts

if( isfield( heasig, 'units' ) )
    
    volt_idx = find(any(cell2mat(cellfun(@(b)(cell2mat(cellfun(@(a)(~isempty(strfind(a, b))), rowvec(upper(cellstr(heasig.units))), 'UniformOutput', false))), colvec(cUnitsVoltages), 'UniformOutput', false)),1));
    
%     volt_idx = intersect( ECG_signals_idx, volt_idx);
    
    for ii = volt_idx
        switch( upper(strtrim(heasig.units(ii,:))) )
            case {'NV', 'NANOVOLTS' , 'NANOVOLTIOS' }
                heasig.gain(ii) = heasig.gain(ii) * 1e3;

            case {'UV', 'MICROVOLTS' , 'MICROVOLTIOS' }

                heasig.gain(ii) = heasig.gain(ii) * 1;

            case {'MV', 'MILIVOLTS' , 'MILIVOLTIOS' }

                heasig.gain(ii) = heasig.gain(ii) / 1e3;

            case {'V', 'VOLTS' , 'VOLTIOS' }

                heasig.gain(ii) = heasig.gain(ii) / 1e6;

            otherwise
                if( bCalledFromWorkspace )
                    cprintf('[1,0.5,0]', 'Unknown units, assuming uV.\n')            
                end

        end
    end
    
    if( length(volt_idx) ~= cant_leads )
        
        
    end
    
else
%     warning('plot_ecg_strip:UnknownADCunits', 'Assuming uV.')
    if( bCalledFromWorkspace )
        cprintf('[1,0.5,0]', 'Unknown units, assuming uV.\n')   
    end
    volt_idx = 1:cant_leads;
    heasig.units = repmat('uV', cant_leads,1);
    
end

aux_val = cellstr(heasig.units);
aux_val(volt_idx) = repmat({'uV'}, length(volt_idx),1);
heasig.units = char(aux_val);

if( ~isfield( heasig, 'btime' ) )
    heasig.btime = '00:00:00';
end

% if( ~isfield( heasig, 'bdate' ) )
%     heasig.bdate = '01/01/0001';
% end

if(sum(heasig.btime == ':') == 2 )
    formatIn = 'HH:MM:SS';
elseif(sum(heasig.btime == ':') == 3 )
    formatIn = 'HH:MM:SS:FFF';
else
    formatIn = [];
    if( bCalledFromWorkspace )
        cprintf('[1,0.5,0]', 'Unknown base time format %s.\n', heasig.btime);
    end
end

if( isempty(formatIn) )
    base_time = 1;
else
    base_time = (datenum( heasig.btime , formatIn) - datenum( '00:00:00' , 'HH:MM:SS')) * 60 * 60 * 24 * heasig.freq; % in samples
end
% Bool to control if the time is measured relative from the begining of the
% recording
bRelative_begin = false;
prev_base_time = base_time;

if( isempty(QRS_start_idx) )
    aux_idx = max(1,round(start_time * heasig.freq)):min(cant_samp, round(end_time * heasig.freq));
    
    if( length(aux_idx) < (min_cant_samp_seconds*heasig.freq) )
        error('plot_ecg_strip:TimeFewSamples', 'Time range should be higher than %3.2f seconds', min_cant_samp_seconds )
    end
    
    if(isempty(aux_idx))
        error('plot_ecg_strip:TimeOutOfBounds', 'Time should be between 0 and %u (%s) seconds', round(cant_samp / heasig.freq), Seconds2HMS(cant_samp / heasig.freq) )
    end
    
else
    
    if( all(cellfun(@(a)(isempty(a)), QRS_locations)) )
%         error('plot_ecg_strip:NoQRSlocations', 'We could not found any QRS complex locations to plot' )
        QRS_start_idx = [];
        QRS_end_idx = [];
        aux_idx = 1:cant_samp;
    else
        QRS_start_idx = max(1, QRS_start_idx);
        QRS_end_idx = min(max(cant_QRS_locations), QRS_start_idx+cant_qrs);
        aux_idx = max(1, min(cellfun( @(a)( protected_index(a,QRS_start_idx)),QRS_locations))) : min(cant_samp, max(cellfun( @(a)(protected_index(a,QRS_end_idx)),QRS_locations) ) );
    end
    
    
end

this_Xmargin = min(1*heasig.freq, round((aux_idx(end) - aux_idx(1))*0.1));

aux_idx =            max(1,aux_idx(1)-this_Xmargin):min(cant_samp, aux_idx(end)+this_Xmargin);

aux_idx_downsample = [ max(1,aux_idx(1)-1*heasig.freq) min(cant_samp, aux_idx(end)+1*heasig.freq) ];

kLinesAnns      = 1;
kBackColourAnns = 2;

ann_graph_mode = kBackColourAnns;

% Visualization Constants:
% Time to start showing detail marks
closeDetailSampSize = 10 * heasig.freq;
mediumDetailSampSize = 30 * heasig.freq;
farDetailSampSize = 60 * heasig.freq;

kNoDetail        = 0;
kCloseDetailSL   = 1;
kMediumDetailSL  = 2;
kCloseDetailML   = 3;
kMediumDetailML  = 4;
kCloseDetailAll  = 5;
kMediumDetailAll = 6;
        
switch(strDetail)
    case 'none'
        eDetailLevel  = kNoDetail;
    case 'single-lead'
        eDetailLevel  = kCloseDetailSL;
    case 'multilead'
        eDetailLevel  = kCloseDetailML;
    otherwise
        eDetailLevel  = kCloseDetailAll;
end

% Also decimate the signal for report generation. Avoid vectoized reports
% too heavy
prev_units = get(fig_hdl, 'Units');
set(fig_hdl, 'Units', 'centimeters');

% Downsample version: For efficient marks visualization and printing only
paper_size = get(fig_hdl, 'Position');
nsamp_target = paper_size(3) * target_res;
set(fig_hdl, 'Units', prev_units);

down_factor = max(1, ceil( (aux_idx(end)-aux_idx(1)+1) / nsamp_target));

if( down_factor > 1 )
    fir_coeffs = design_downsample_filter(down_factor);

    aux_idx_downsample_start = max(1, round((aux_idx(1) - aux_idx_downsample(1) + 1)/down_factor));
    aux_idx_downsample_end = length( linspace(aux_idx(1),aux_idx(end),round((aux_idx(end)-aux_idx(1)+1)/down_factor)) ) + aux_idx_downsample_start - 1;

    if( isempty(ECG) )

        % envelope resampling
%         aux_val1 = double(ECG_w.read_signal(aux_idx_downsample(1), aux_idx_downsample(end)));
%         aux_val1 = hilbert(aux_val1);
%         aux_val2 = abs(aux_val1);
%         ECGd = resample(aux_val2, 1, down_factor, fir_coeffs);
        
        % just decimation
%         aux_val1 = double(ECG_w.read_signal(aux_idx_downsample(1), aux_idx_downsample(end)));
%         aux_val2 = aux_idx_downsample(end) - aux_idx_downsample(1) + 1;
%         ECGd = aux_val1(round(linspace(1,aux_val2,round(aux_val2/down_factor))),:);
        % anti alias - decimation
        ECGd = resample(double(ECG_w.read_signal(aux_idx_downsample(1), aux_idx_downsample(end))), 1, down_factor, fir_coeffs);

        ECGd = ECGd(aux_idx_downsample_start:aux_idx_downsample_end,:);
        ECG = ECG_w.read_signal(aux_idx(1), aux_idx(end));

        [cant_samp, cant_leads] = size(ECG);

        %transform to real units
        ECG = bsxfun( @rdivide, bsxfun( @minus, double(ECG), rowvec(double(heasig.adczero)) ), rowvec(double(heasig.gain)) ) ;
        ECGd = bsxfun( @rdivide, bsxfun( @minus, ECGd, rowvec(double(heasig.adczero)) ), rowvec(double(heasig.gain)) ) ;

    else
        ECGd = resample( ECG(aux_idx_downsample(1):aux_idx_downsample(end),:), 1, down_factor, fir_coeffs );
        ECGd = ECGd(aux_idx_downsample_start:aux_idx_downsample_end,:);
        ECGd = bsxfun( @rdivide, bsxfun( @minus, double(ECGd), rowvec(double(heasig.adczero)) ), rowvec(double(heasig.gain)) ) ;

        %transform to real units
        ECG = bsxfun( @rdivide, bsxfun( @minus, double(ECG(aux_idx,:)), rowvec(double(heasig.adczero)) ), rowvec(double(heasig.gain)) ) ;
        [cant_samp, cant_leads] = size(ECG);
        
    end
    
else
    % no resampling needed.
    
    if( isempty(ECG) )
        ECG = ECG_w.read_signal(aux_idx(1), aux_idx(end));
        ECGd = ECG;
        
        [cant_samp, cant_leads] = size(ECG);

        %transform to real units
        ECG = bsxfun( @rdivide, bsxfun( @minus, double(ECG), rowvec(double(heasig.adczero)) ), rowvec(double(heasig.gain)) ) ;
        ECGd = bsxfun( @rdivide, bsxfun( @minus, double(ECGd), rowvec(double(heasig.adczero)) ), rowvec(double(heasig.gain)) ) ;
    else
        %transform to real units
        ECG = bsxfun( @rdivide, bsxfun( @minus, double(ECG(aux_idx,:)), rowvec(double(heasig.adczero)) ), rowvec(double(heasig.gain)) ) ;
        ECGd = ECG;
        cant_samp = length(aux_idx);
    end
end
    
if(bFilterECG)
    if( isempty(ECG_signals_idx) )
        warning('plot_ecg_strip:NotFilter', 'No filter was applied since no ECG leads found.')
    else
        ECG(:,ECG_signals_idx) = BaselineWanderRemovalMedian( ECG(:,ECG_signals_idx), heasig.freq);
    end
    
%     filtro = bandpass_filter_design( heasig.freq );
%     ECG(:,ECG_signals_idx) = filter(filtro, flipud() );
%     ECG(:,ECG_signals_idx) = filter(filtro, flipud(ECG(:,ECG_signals_idx)) );
end

% ECG ranges
if( all(cellfun(@(a)(isempty(a)), QRS_locations)) )
%     if( isempty(global_annotations) )
        this_qrs_ploted = [];
%     else
%         this_qrs_ploted = find(global_annotations(1).qrs >= (aux_idx(1) + this_Xmargin) & global_annotations(1).qrs <= (aux_idx(end) - this_Xmargin) );
%         QRS_locations = {global_annotations(1).qrs};
%     end
    this_qrs_ploted = {this_qrs_ploted};
else
    this_qrs_ploted = cellfun( @(a)(find(a >= aux_idx(1) & a <= aux_idx(end) )), QRS_locations, 'UniformOutput', false );
    if( ~isempty(QRS_start_idx) && ~isempty(cant_qrs) )
        this_qrs_ploted = cellfun( @(a,b)(intersect(a, QRS_start_idx:min(b, QRS_start_idx+cant_qrs))), this_qrs_ploted, num2cell(cant_QRS_locations), 'UniformOutput', false ); 
    end
end

[ecg_range, ecg_min, ecg_max, ecg_median] = CalcECG_range(ECG);

if( isempty(gains) )
    gains = ones( cant_leads, 1);
    bAux = ecg_range ~= 0;
    gains(bAux) = max(ecg_range(bAux))./ecg_range(bAux);
else
    % by the moment one gain for all
    if(length(gains) == 1 )
        gains = repmat( gains, cant_leads, 1);
    elseif(length(gains) ~= cant_leads )
        error('plot_ecg_strip:BadLeadGain', ['Lead gain must be a single numeric value or a ' num2str(cant_leads) ' x 1 vector'] )
    end
end

% ensure good visibility
if( isempty(lead_offset) )
    
    if( cant_leads > 1 )
        lead_offset = [0; 1.1*( (ecg_min(1:end-1) - ecg_median(1:end-1))  .* gains(1:end-1) - ( ecg_max(2:end) - ecg_median(2:end) ) .* gains(2:end) ) ];
        offsets = ecg_median .* gains + abs(cumsum(lead_offset));
    else
        lead_offset = 0;
        offsets = 0;
    end
else
    if(length(lead_offset) == 1 )
        lead_offset = repmat( lead_offset, cant_leads, 1);
    elseif(length(lead_offset) ~= cant_leads )
        error('plot_ecg_strip:BadLeadOffset', ['Lead offset must be a single numeric value or a ' num2str(cant_leads) ' x 1 vector'] )
    end
end

% signal margins, relative to the total height and width
plot_left_margin_width = 0.03;
plot_rigth_margin_width = 0.035;
plot_top_margin_width = 0.1;
plot_bottom_margin_width = 0.1;

[plotYrange, plotYmin, plotYmax] = CalcPlotYlimits(ecg_min, ecg_max, gains, offsets);

plotXmin = aux_idx(1); 
plotXmax = aux_idx(end);
plotXrange = plotXmax - plotXmin;
plotXmin = plotXmin - plot_left_margin_width * plotXrange; 
plotXmax = plotXmax + plot_rigth_margin_width * plotXrange; 

% Set the colormap
ColorOrder = my_colormap( 12 );
ColorOrder = repmat(ColorOrder, ceil(heasig.nsig/size(ColorOrder,1)), 1 );

set(axes_hdl, 'ColorOrder', ColorOrder);

hold(axes_hdl, 'on')

%plot ECG scaled and translated
ECG_hdl = plot(axes_hdl, aux_idx, bsxfun( @minus, bsxfun( @times, ECG, rowvec(gains)), rowvec(offsets) ), 'LineWidth', 1.3 );

% ylabel(['ECG (' heasig.units(1,:) ')' ]);
% xlabel(['Time (s) relative to start of segment' ]);
set(axes_hdl, 'Box', 'off' );
set(axes_hdl, 'Xtick', [] );
set(axes_hdl, 'Ytick', [] );
set(axes_hdl, 'Xcolor', [1 1 1] );
set(axes_hdl, 'Ycolor', [1 1 1] );
% set(axes_hdl, 'Visible', 'off' );
% set(axes_hdl, 'Visible', 'off' );
% set(axes_hdl, 'Visible', 'off' );

set(axes_hdl, 'Position', [ 0.005 0.01 0.99 0.98 ] );

% legend(ECG_hdl, cellstr(heasig.desc) );

% user data initialized
user_data.linked_hdl = linked_hdl;

% global variables
report_format = 'png';
report_format_idx = find(strcmpi(report_format,cKnownReportFormats),1);

lead_selected_idx = 1:cant_leads;

bPaperModeOn = false;
major_tick_values_time = round([0.1 0.2 0.5 1 2 5 10 30 60]*heasig.freq); % seconds
major_tick_values_time = unique([major_tick_values_time major_tick_values_time*60 ]);
major_tick_values_voltage = [ [1 2 5 ] * 10^(-3) [1 2 5 ] * 10^(-2) [1 2 5 ] * 10^(-1)  [1 2 5 ] [1 2 5 ] * 10^1 [1 2 5 ] * 10^2 [1 2 5] * 10^3 [1 2 5] * 10^4 [1 2 5] * 10^5 [1 2 5] * 10^6  ]; % seconds
paperModeHdl = {};
topLevelHdl = {};
ECGd_hdl = [];

bPlotScales = false;
PrevStateWindowButtonMotionFcn = [];
UserChnageViewHdls = {};
%store user data.
start_sample = aux_idx(1);
end_sample = aux_idx(end);
% qrs_ploted = cellfun( @(a,b)(a(b)), QRS_locations, this_qrs_ploted, 'UniformOutput', false );
qrs_ploted = QRS_locations;
qrs_ploted_names = QRS_locations_names;
prev_Xrange = plotXmax - plotXmin;
prev_Yrange = plotYmax - plotYmin;
original_plot_lims = [ plotXmin plotXmax plotYmin plotYmax ];
plotYtimeAxis = [];
plotYindexesAxis = [];
PwaveHdls = {};
TwaveHdls = {};
QRScplxHdls = {};
PwaveGlblHdls = cell(heasig.nsig,1);
TwaveGlblHdls = cell(heasig.nsig,1);
QRScplxGlblHdls = cell(heasig.nsig,1);
QRSfpHdls = {};
QRSfpFarHdls = {};
titleExtent = nan;
timeAxisTop = nan;
startSignalX = nan;
endSignalX = nan;
xTextOffset = nan;
yTextOffset = nan;
YlabelLeftPosition = nan;
scroll_mode = 'zoom';
gain_offset_mode = 'gain';
% 1V maximum
min_gain = 1e-12 / plotYrange;
% 1 nV maximum
max_gain = 1e15 / plotYrange;

my_timer = timer('TimerFcn',@timer_fcn, 'StopFcn', @timer_stop_fcn , 'StartDelay', 10);

%% variables and flags from Dragzoom

fNoZoom = false;

hAxes = axes_hdl;
% get info about all axes and create axes info struct
mAxesInfo = GetAxesInfo();
fIsSelectedCurrentAxes = true;

% Drag Options
mDragKeysX = 'normal';      % 'normal', 'reverse'
mDragKeysY = 'normal';      % 'normal', 'reverse'
mDragShiftStep = 3;         % step dragging on keys
mDragShiftStepInc = 1;      % increase speed dragging on keys
mDragStartX = [];
mDragStartY = [];
mDragShiftStepInc = [];

% Zoom Options
mZoomScroll = 'normal';     % 'normal', 'reverse'
mZoomMinPow = 0;            % min zoom percent 10 ^ mZoomMinPow
mZoomMaxPow = 5;            % max zoom perzent 10 ^ mZoomMaxPow
mZoomNum = 51;              % count steps of log zoom grid
mZoomExtendNum = 301;       % count steps of log grid zoom extend for 2D
mZoomKeysNum = 181;         % count steps of log grid zoom for keys for 2D

% Rubber Band Options
mRbEdgeColor = 'k';         % rubber band edge color
mRbFaceColor = 'none';      % rubber band face color
mRbFaceAlpha = 1;           % rubber band face alpha (transparency)
mRubberBand = [];

% Magnifier Options
mMgSize = 100;              % default size of magnifier (pixels)
mMgMinSize = 50;            % min size of magnifier
mMgMaxSize = 200;           % max size of magnifier
mMgZoom = 2;                % default zoom on magnifier
mMgMinZoom = 1;             % min zoom on magnifier
mMgMaxZoom = 100;           % max zoom on magnifier
mMgLinesWidth = 1;          % lines width on magnifier
mMgShadow = 0.95;           % shadow area without magnifier
mMgSizeStep = 15;           % step change in the magnifier size
mMgZoomStep = 1.2;          % step change in the magnifier zoom
mMagnifier = [];
mMgDirection = [];


mDefaultZoomGrid = [];
mDefaultZoomSteps = [];
mZoomGrid = [];
mZoomSteps = [];
mZoomIndexX = [];
mZoomIndexY = [];
mZoom3DStartX = [];
mZoom3DStartY = [];
mZoom3DBindX = [];
mZoom3DBindY = [];
mDefaultXLim = original_plot_lims(1:2);
mDefaultYLim = original_plot_lims(3:4);
mPointerCross = [];


% flags
% flags
bgColor = [251 248 230]/255;

fIsSelectedCurrentAxes = true;
fIsDragAllowed = false;
fIsRubberBandOn = false;
fIsPointerCross = false;
fIsAxesGrid = false;
fIsEnableZoomX = false;
fIsEnableZoomY = false;
fIsMagnifierOn = false;
fIsEnableControl = true;
fIsMouseOnLegend = false;

SetDefaultZoomGrid();
mDragSaveShiftStep = mDragShiftStep;

%% single lead annotations

if( isempty(annotations) )
    annotations = [];
    this_annotation = [];
else
    for ii = 1:length(annotations)
        this_annotation = annotations(ii);
        aux_struct = [];
        
        if( isempty(this_qrs_ploted) )
            this_qrs_ploted = find(this_annotation.qrs >= (aux_idx(1) + this_Xmargin) & this_annotation.qrs <= (aux_idx(end) - this_Xmargin) );

            if( ~isempty(this_qrs_ploted) )

                for field_name = cAnnotationFieldNames
                    aux_val = this_annotation.(field_name{1});
                    if( isempty(aux_val) )
                        aux_struct.(field_name{1}) = nan(length(this_qrs_ploted),1);
                    else
                        aux_struct.(field_name{1}) = aux_val(this_qrs_ploted);
                    end
                end
                annotations(ii) = aux_struct;
            end
        end
    end
    
end

%% multilead annotations

if( isempty(global_annotations) )
    global_annotations = [];
else
    
    if( isempty(this_qrs_ploted) )
        this_qrs_ploted = find(global_annotations.qrs >= (aux_idx(1) + this_Xmargin) & global_annotations.qrs <= (aux_idx(end) - this_Xmargin) );
    
        if( isempty(this_qrs_ploted) )
            global_annotations = [];
            this_qrs_ploted = [];
        else

            for field_name = cAnnotationFieldNames
                if( isfield(global_annotations, field_name{1} ) )
                    aux_val = global_annotations.(field_name{1});
                    aux_struct.(field_name{1}) = aux_val(this_qrs_ploted);
                end
            end
            global_annotations = aux_struct;

            this_qrs_ploted = {global_annotations.qrs};

        end
        
    end
    
end

set(axes_hdl, 'Xlim', [plotXmin plotXmax]); 
if( plotYmin < plotYmax )
    set(axes_hdl, 'Ylim', [plotYmin plotYmax]); 
end

% this_qrs_ploted = cellfun( @(a)( find(a >= plotXmin & a <= plotXmax) ), qrs_ploted, 'UniformOutput', false );

if( isfield(heasig, 'recname') )
    recname = heasig.recname;
else
    recname = 'no_name';
end

if( isempty(strTitle))
    strTitle = [ 'Recording ' recname ' - ' Seconds2HMS( (aux_idx(1) + base_time ) / heasig.freq) ' : ' Seconds2HMS( (aux_idx(end) + base_time ) / heasig.freq) ];
else
    strTitle = strTitle;
end

update_axis_plot_ecg_strip( this_qrs_ploted );

uistack(ECG_hdl, 'top');

hold(axes_hdl, 'off');

if(bPrettyPrint)
    PaperModeOn()
end

% user data stored in figure.
set(fig_hdl, 'UserData', user_data);

% %Zoom and Pan handles
% zoom_hdl = zoom(fig_hdl);
% zoom reset;
% set(zoom_hdl, 'ActionPostCallback', @(a,b)( UserChangeView(a,b, 'zoom')) );
% 
% pan_hdl = pan(fig_hdl);
% set(pan_hdl, 'ActionPostCallback', @(a,b)( UserChangeView(a,b, 'pan')) );
% 
% set(fig_hdl,'WindowScrollWheelFcn',@ScrollWheel);            
% % set(fig_hdl,'WindowKeyPressFcn',@KeyPress);            
% addlistener(fig_hdl,'WindowKeyPressEvent', @KeyPress);

    set(fig_hdl,'CloseRequestFcn',@my_closefcn)

set(fig_hdl, ...
    'WindowButtonDownFcn',      {@WindowButtonDownCallback2D}, ...
    'WindowButtonUpFcn',        {@WindowButtonUpCallback2D}, ...
    'WindowScrollWheelFcn',     {@WindowScrollWheelFcn2D}, ...
    'WindowKeyPressFcn',        {@WindowKeyPressCallback2D}, ...
    'CloseRequestFcn',          {@my_closefcn}, ...
    'WindowKeyReleaseFcn',      {@WindowKeyReleaseCallback2D});

if(~bPrettyPrint)
    disp_help()
end

if( nargout > 1 )
    returned_handles = ECG_hdl;
else
    clear returned_handles
end

%==========================================================================

%--------------------------------------------------------------------------

    function update_axis_plot_ecg_strip( qrs2plot )

%   Some useful constans

    xTextOffset = 0.005*plotXrange;
    yTextOffset = 0.01*plotYrange;

    if( bPaperModeOn )
        
        startSignalX = (plotXmin + xTextOffset);
        endSignalX   = (plotXmax - xTextOffset);
        
        width_legend = 0.07*plotXrange;
        height_legend = 0.03*plotYrange;
        left_legend = plotXmax - width_legend - 4*xTextOffset;
        bottom_legend = plotYmax - height_legend - 3*yTextOffset;

        titleYposition = plotYmax - 4*yTextOffset;
        
    else
        startSignalX = (plotXmin + 5*xTextOffset);
        endSignalX   = (plotXmax - 3*xTextOffset);
        
        width_legend = 0.07*plotXrange;
        height_legend = 0.03*plotYrange;
        left_legend = plotXmax - width_legend - 2*xTextOffset;
        bottom_legend = plotYmax - height_legend - 1*yTextOffset;
        
        titleYposition = plotYmax - 2*yTextOffset;
        
    end
    
    timeAxisTop = plotYmin + (0.06 * plotYrange);
    plotYtimeAxis = plotYmin + (0.04 * plotYrange);

    

    %% initialization
    
    topLevelHdl = [];
    
    %% Lead reference:
%   Plot leads label or order number close to the signals.
%   Check overlapps among signal labels.

    % text label of each axis
    
    if( bPaperModeOn )
        YlabelLeftPosition = startSignalX + 6*xTextOffset;
    else
        YlabelLeftPosition = startSignalX - 3*xTextOffset;
    end
    
    if( isfield(heasig, 'desc') )
        if( bPaperModeOn )
            str_aux = cellstr(strtrim(heasig.desc));
        else
            str_aux = strcat(cellstr(strtrim(heasig.desc)), repmat({' ('},heasig.nsig,1), cellstr(heasig.units), repmat({')'},heasig.nsig,1) );
        end
    else
        str_aux =  cellstr(num2str(colvec(1:cant_leads)));
    end

    % signal labels
    if( bPaperModeOn )
        aux_hdl = cellfun( @(a,b,ii)( text( YlabelLeftPosition, b, a, 'FontSize', 8, 'HorizontalAlignment', 'left', 'Interpreter', 'none', 'BackgroundColor', [1 1 1], 'EdgeColor', ColorOrder(ii,:), 'Margin', 3 ) ), str_aux, num2cell(-offsets + ecg_max .* gains),num2cell((1:heasig.nsig)'), 'UniformOutput', false);
    else
        % white patch to cover the background        
        UserChnageViewHdls = [UserChnageViewHdls; {patch('Faces', [1 2 3 4], 'Vertices', [[0; 0; 3*xTextOffset; 3*xTextOffset ] + plotXmin [plotYmin; plotYmax - plotYmin; plotYmax - plotYmin;plotYmin ] ], 'FaceColor', [1 1 1], 'EdgeColor', [1 1 1] )} ];
        aux_hdl = cellfun( @(a,b,ii)( text( YlabelLeftPosition, -b + (ecg_min(ii) + ecg_range(ii)/2) * gains(ii), a, 'FontSize', 8, 'HorizontalAlignment', 'center', 'Rotation', 90, 'Interpreter', 'none', 'BackgroundColor', [1 1 1], 'EdgeColor', ColorOrder(ii,:), 'Margin', 3 ) ), str_aux, num2cell(offsets),num2cell((1:heasig.nsig)'), 'UniformOutput', false);
    end
    
    UserChnageViewHdls = [UserChnageViewHdls; colvec(aux_hdl)];
    
    % r2014 change handle to objects.
    aux_hdll = aux_hdl{1};
    for iii = 2:length(aux_hdl)
        aux_hdll(iii) = aux_hdl{iii};
    end
    aux_hdl = aux_hdll;
    
    if( length(aux_hdl) > 1 )
        aux_extent = cell2mat(get(aux_hdl, 'Extent'));
    else
        aux_extent = get(aux_hdl, 'Extent');
    end
    
    if( size(aux_extent,1) > 1 && sum(aux_extent(:,4)) > plotYrange )
        % not enough room for all descriptions -> overlap
        % fix first and last, and share the rest of space
        aux_val = cell2mat(get(aux_hdl, 'Position'));
        set(aux_hdl(1),   'Position', [aux_val(1,1) aux_val(1,2) - ( aux_extent(1,2) - (plotYmax - aux_extent(1,4))) ] )
        set(aux_hdl(end), 'Position', [aux_val(end,1) aux_val(end,2) - ( aux_extent(end,2) - (plotYmin + aux_extent(end,4))) ] )
        
        if( size(aux_extent,1) > 2 )
            aux_val = cell2mat(get(aux_hdl, 'Position'));
            if( length(aux_hdl) > 1 )
                aux_extent = cell2mat(get(aux_hdl, 'Extent'));
            else
                aux_extent = get(aux_hdl, 'Extent');
            end
            aux_space = aux_extent(1,2) - (aux_extent(end,2)+aux_extent(end,4));
            each_weight = aux_extent(2:end-1,4);
            each_weight = each_weight ./ sum(each_weight);
            each_space = aux_space .* each_weight;
            arrayfun( @(ii,a)( set(aux_hdl(ii), 'Position', [aux_val(ii,1) aux_val(ii,2) - ( aux_extent(ii,2) - (aux_extent(1,2) - a) ) ]) ), 1+(1:length(each_space))' , cumsum(each_space) )
        end
        
    end
    
    
    % save this handles to top them later
    topLevelHdl = [topLevelHdl; colvec(UserChnageViewHdls(end-cant_leads+1:end))];


    if( ~bPaperModeOn )

        % vertical axes for each signal

        %check overlapp
        if( length(offsets) > 1 )
            aux_solap = [ (-offsets(1:end-1) - colvec(ecg_range(1:end-1))/2) (-offsets(2:end) + colvec(ecg_range(2:end))/2) ];
            aux_solap = aux_solap(:,1) < aux_solap(:,2);
            aux_solap = [ aux_solap(1); aux_solap];
            aux_jitter = aux_solap .* colvec(linspace(0,2,cant_leads) * xTextOffset);
        else
            aux_jitter = 0;
        end
        
        %build vertical axis
        aux_Xaxis = [ -xTextOffset 0 0 -xTextOffset] + startSignalX;
        aux_Yaxis = [ ecg_max ecg_max ecg_min ecg_min ]';
        UserChnageViewHdls = [UserChnageViewHdls; colvec(arrayfun(@(a,b)(plot(axes_hdl, a, b, 'LineWidth', 0.25 )), bsxfun(@plus, repmat(colvec(aux_Xaxis),1,cant_leads), rowvec(aux_jitter)), bsxfun( @minus, bsxfun( @times, aux_Yaxis, rowvec(gains) ), rowvec(offsets)), 'UniformOutput', false )) ];       

    end
       
    % user_data.UserChnageViewHdls = [user_data.UserChnageViewHdls; text( sample_reference + xTextOffset, plotYmin + (0.02 * plotYrange), 'Index #' , 'FontSize', 8)];

    %% Time reference: 
%     Start/End of the excerpt.
%     QRS location/indexes if were provided
    
    if( plotXrange < mediumDetailSampSize )

        sample_reference = ceil( (plotXmin + min( heasig.freq, plot_left_margin_width*plotXrange )) * 1/heasig.freq) * heasig.freq;
        
        bAux = any( colvec(cellfun( @(a)(isempty(a)), qrs2plot ) ));
        if( bAux )
            if( sample_reference < plotXmin || sample_reference > plotXmin + plot_left_margin_width*plotXrange )
                sample_reference = ceil(plotXmin + min( heasig.freq, plot_left_margin_width*plotXrange));
            end
        else
            aux_val = min( cellfun( @(a,b)(a(b(1))), qrs_ploted, qrs2plot ) );
            if( sample_reference >= aux_val ) 
                sample_reference = ceil(plotXmin + min( heasig.freq, plot_left_margin_width*plotXrange));
            end
        end
        
        sample_last_reference = floor( (plotXmax - min( heasig.freq, plot_rigth_margin_width*plotXrange)) * 1/heasig.freq) * heasig.freq ;
        bAux = any( colvec(cellfun( @(a)(isempty(a)), qrs2plot ) ));
        if( bAux )
            if( sample_last_reference > plotXmax || sample_last_reference < plotXmax - min( heasig.freq,plot_rigth_margin_width*plotXrange)  )
                sample_last_reference = floor(plotXmax - min( heasig.freq, plot_rigth_margin_width*plotXrange));
            end
        else
            aux_val = max( cellfun( @(a,b)(a(b(end))), qrs_ploted, qrs2plot ) );
            if( sample_last_reference <= aux_val )
                sample_last_reference = floor(plotXmax - min( heasig.freq,plot_rigth_margin_width*plotXrange));
            end
        end
        
        if( ~bPaperModeOn )
            
            %obsolete: try to handle visibility with uistack
%             % plot a white patch as background
%             UserChnageViewHdls = [UserChnageViewHdls; patch('Faces', [1 2 3 4], 'Vertices', [[0; plotXrange; plotXrange; 0 ] + plotXmin [0; 0; timeAxisTop; timeAxisTop ] + plotYmin ], 'FaceColor', [1 1 1], 'EdgeColor', [1 1 1] ) ];

            aux_hdl = text( sample_reference - xTextOffset, plotYmin + (0.02 * plotYrange), ... 
                            'Index #' , ...
                            'FontSize', 8, ...
                            'HorizontalAlignment', 'right' ...
                            );
            aux_extent = get(aux_hdl, 'Extent');
            if(aux_extent(1) < plotXmin )
                aux_position = get(aux_hdl, 'Position');
                set(aux_hdl, 'Position', [ plotXmin+aux_extent(3) aux_position(2:3)]);
            end
            UserChnageViewHdls = [UserChnageViewHdls; {aux_hdl}];

            % time reference

            % start - end of ECG excerpt
            aux_hdl = text( sample_reference - xTextOffset, plotYtimeAxis, Seconds2HMS( (sample_reference + base_time ) / heasig.freq, 3), 'FontName', 'Arial', 'FontSize', 8, 'HorizontalAlignment', 'right');
            aux_extent = get(aux_hdl, 'Extent');
            if(aux_extent(1) < plotXmin )
                aux_position = get(aux_hdl, 'Position');
                set(aux_hdl, 'Position', [ plotXmin+aux_extent(3) aux_position(2:3)]);
            end
            UserChnageViewHdls = [UserChnageViewHdls; {aux_hdl}];

            aux_hdl = text( sample_last_reference + xTextOffset, plotYtimeAxis, Seconds2HMS((sample_last_reference + base_time ) / heasig.freq, 3), 'FontName', 'Arial', 'FontSize', 8);
            aux_extent = get(aux_hdl, 'Extent');
            if( aux_extent(1)+aux_extent(3) > plotXmax )
                aux_position = get(aux_hdl, 'Position');
                set(aux_hdl, 'Position', [ plotXmax - aux_extent(3) aux_position(2:3)]);
            end
            UserChnageViewHdls = [UserChnageViewHdls; {aux_hdl}];

            aux_seq = cellfun( @(a,b)(limit_sequence(a(b),[sample_reference sample_last_reference])), qrs_ploted, qrs2plot, 'UniformOutput', false );

            % plot indexes
            plotYindexesAxis = plotYmin + (0.02 * plotYrange);

            cHorAllignment = {'Left' 'Right' 'Center' };
            laux_seq = length(aux_seq);
            aux_hallign = repmat(cHorAllignment, 1, ceil(laux_seq/3) );
            aux_hallign = aux_hallign(1:laux_seq);

            linespec_color = rowvec(cellfun( @(a)(a{4}), cLinespecs(1:length(aux_seq)), 'UniformOutput', false ));

            aux_vall = cellfun( @(a,b,c,e,f,g)( colvec(arrayfun( @(d)(text( a(b(d)), plotYindexesAxis + e, num2str(b(d)) , 'FontSize', 8, 'HorizontalAlignment', f, 'Color', g)), c, 'UniformOutput', false)) ), qrs_ploted, qrs2plot, aux_seq, num2cell(0.02 * plotYrange * linspace(-0.1, 0.1, length(qrs_ploted))), aux_hallign, linespec_color, 'UniformOutput', false );
            % r2014 change handle to objects.
            laux_seq = aux_vall{1};
            for iii = 2:length(aux_vall)
                laux_seq = [laux_seq; colvec(aux_vall{iii})];
            end
            
            UserChnageViewHdls = [ ... 
                                UserChnageViewHdls; laux_seq
                                ];
            
            % plot times relative to the reference
            aux_vall = colvec(cellfun( @(a,b,c,e,f,g)( colvec(arrayfun( @(d)(text( a(b(d)), plotYtimeAxis + e, Seconds2HMS( (a(b(d)) - sample_reference + 1)/heasig.freq, 3 ) , 'FontSize', 8, 'HorizontalAlignment', f, 'Color', g )), c, 'UniformOutput', false)) ), qrs_ploted, qrs2plot, aux_seq, num2cell(0.02 * plotYrange * linspace(-0.1, 0.1, length(qrs_ploted))), aux_hallign, linespec_color, 'UniformOutput', false ));
            % r2014 change handle to objects.
            laux_seq = aux_vall{1};
            for iii = 2:length(aux_vall)
                laux_seq = [laux_seq; colvec(aux_vall{iii})];
            end
            
            UserChnageViewHdls = [ ... 
                                UserChnageViewHdls; laux_seq
                                ];
                            
            % black tips at the beginning/end
            UserChnageViewHdls = [UserChnageViewHdls; colvec(arrayfun(@(a,b)(plot(axes_hdl, a, b , 'k-', 'LineWidth', 0.25 )), [sample_reference sample_last_reference; sample_reference sample_last_reference], [ repmat(plotYmin + (0.01 * plotYrange), 1, 2); repmat(plotYmin + (0.06 * plotYrange), 1,2) ], 'UniformOutput', false)) ];       
        
        end
        
    else

        % for far zoom views

        sample_reference = ceil( (plotXmin + plot_left_margin_width*plotXrange ) * 1/heasig.freq) * heasig.freq;
        sample_last_reference = floor( (plotXmax - plot_rigth_margin_width*plotXrange) * 1/heasig.freq) * heasig.freq ;

        if( ~bPaperModeOn )

            aux_seq = linspace( sample_reference, sample_last_reference, 6);

            aux_seq = round( aux_seq/heasig.freq ) * heasig.freq ;

            plotYtimeAxis = plotYmin + (0.02 * plotYrange);

            UserChnageViewHdls = [ UserChnageViewHdls; {text( aux_seq(1) - xTextOffset, plotYtimeAxis, Seconds2HMS( (aux_seq(1) + base_time )/heasig.freq ) , 'FontSize', 8, 'HorizontalAlignment', 'right') } ];
            for jj = rowvec(aux_seq(2:end-1))
                UserChnageViewHdls = [ UserChnageViewHdls; {text( jj, plotYtimeAxis, Seconds2HMS( (jj-aux_seq(1))/heasig.freq ) , 'FontSize', 8, 'HorizontalAlignment', 'center')}];
            end
            UserChnageViewHdls = [ UserChnageViewHdls; {text( aux_seq(end) + xTextOffset, plotYtimeAxis, Seconds2HMS( (aux_seq(end) + base_time )/heasig.freq ) , 'FontSize', 8, 'HorizontalAlignment', 'left')}];

            % black tips at the beginning/end
            UserChnageViewHdls = [UserChnageViewHdls; colvec(arrayfun(@(a,b)(plot(axes_hdl, a, b, 'k-', 'LineWidth', 0.25 )), [sample_reference sample_last_reference; sample_reference sample_last_reference], [ repmat(plotYmin + (0.01 * plotYrange), 1, 2); repmat(plotYmin + (0.06 * plotYrange), 1,2) ], 'UniformOutput', false )) ];       
            
        end
    end

    bWaveLegendPlotted = false;
    bHBclassLegendPlotted = false;
    
    %% multilead or global annotations
    
    if( isempty(global_annotations) )

        aux_seq = cellfun( @(a)({1:length(a)}), qrs2plot );
        
        % for far zoom views
        if( ~isempty(QRSfpFarHdls) )
            cellfun(@(a)( CheckAndDeleteHdl(a) ), QRSfpFarHdls);
        end
        
        bAux = eDetailLevel ~= kNoDetail && ( ( eDetailLevel == kMediumDetailML || eDetailLevel == kMediumDetailAll  ) && ( plotXrange > mediumDetailSampSize && plotXrange < farDetailSampSize ) );
        if( bAux )
            QRSfpFarHdls = cellfun( @(a,b,c,d)( plot(axes_hdl, rowvec(a(b(c))), repmat(titleYposition - (d* 0.1 * plotYrange), 1,length(c)) )), qrs_ploted, qrs2plot, aux_seq, num2cell(linspace(0.9, 1.1, length(aux_seq))), 'UniformOutput', false );
            cellfun( @(a,b)(set_a_linespec(a, b)), QRSfpFarHdls, cLinespecsNone(1:length(QRSfpFarHdls)) );
        end

        % for closer zoom views
        if( ~isempty(QRSfpHdls) ) 
            cellfun(@(a)( CheckAndDeleteHdl(a) ), QRSfpHdls);
        end

        bAux = eDetailLevel ~= kNoDetail && ( (eDetailLevel == kCloseDetailML || eDetailLevel == kCloseDetailAll ) && plotXrange < mediumDetailSampSize);
        if( bAux )
            bAux2 = any( colvec(cellfun( @(a,b,c)(~isempty(a(b(c)))), qrs_ploted, qrs2plot, aux_seq ) ));
            if( bAux2 )
                QRSfpHdls = colvec(cellfun( @(a,b,c,d)( plot(axes_hdl, repmat(rowvec(a(b(c))), 2,1), [ repmat(plotYmin + d*(0.06 * plotYrange), 1,length(c)); repmat(titleYposition - d*(0.1 * plotYrange), 1,length(c)) ] ) ), qrs_ploted, qrs2plot, aux_seq, num2cell(linspace(0.9, 1.1, length(aux_seq))), 'UniformOutput', false));
                cellfun( @(a,b)(set_a_linespec(a, b)), QRSfpHdls, cLinespecs(1:length(QRSfpHdls)) );
                cellfun(@(a)(set(a, 'LineWidth', 0.25)), QRSfpHdls);
            end
        end
        
    else
        
        bAux = eDetailLevel ~= kNoDetail && ( (eDetailLevel == kCloseDetailML || eDetailLevel == kCloseDetailAll ) && plotXrange <= closeDetailSampSize );
        
        if( bAux )
            
            bWaveLegendPlotted = true;
            % Waves frame
            legend_hdl = {patch([left_legend left_legend [left_legend left_legend]+width_legend left_legend ], [bottom_legend [bottom_legend bottom_legend]+height_legend bottom_legend bottom_legend], [1 1 1], 'EdgeColor', [0 0 0])};
            legend_hdl = [legend_hdl; {text( left_legend + width_legend/4, bottom_legend + height_legend/2, 'P', 'FontSize', 8, 'HorizontalAlignment', 'center', 'BackgroundColor', color_Pwave_global)}];
            legend_hdl = [legend_hdl; {text( left_legend + width_legend/2, bottom_legend + height_legend/2, 'QRS', 'FontSize', 8, 'HorizontalAlignment', 'center', 'BackgroundColor', color_QRScomplex_global )}];
            legend_hdl = [legend_hdl; {text( left_legend + width_legend*3/4, bottom_legend + height_legend/2, 'T', 'FontSize', 8, 'HorizontalAlignment', 'center', 'BackgroundColor', color_Twave_global )}];
            UserChnageViewHdls = [UserChnageViewHdls; legend_hdl];
            
        end
        
        aux_seq = cellfun( @(a)({1:length(a)}), qrs2plot );
        
        % when annotations are provided
        % P wave
        bAux = eDetailLevel ~= kNoDetail && ( ( eDetailLevel == kCloseDetailML || eDetailLevel == kCloseDetailAll ) && plotXrange <= closeDetailSampSize );
        if( isempty(PwaveHdls) )
            if( bAux )
%                 PwaveHdls = [ PwaveHdls; PlotGlobalWaveMarks({'Pon' 'P' 'Poff'}, [ plotYmin + (0.1*plotYrange)  plotYmax - (0.05 * plotYrange) ], color_Pwave_global )];
                for jj = ECG_signals_idx
                    PwaveHdls = [ PwaveHdls; PlotWaveMarks(global_annotations, {'Pon' 'P' 'Poff'}, jj, 0.5*yTextOffset, PwaveColor*0.95 ) ];
                end
            end
        else
            if( bAux )
                str_aux = 'on';
            else
                str_aux = 'off';
            end
            cellfun( @(a)(set(a, 'Visible', str_aux )), PwaveHdls);
        end

        % QRS complex

        % for far zoom views
        bAux = eDetailLevel ~= kNoDetail && ( ( eDetailLevel == kMediumDetailML || eDetailLevel == kMediumDetailAll ) && plotXrange > mediumDetailSampSize && plotXrange < farDetailSampSize );
        if( isempty(QRSfpFarHdls) )
            if( bAux )
                QRSfpFarHdls = colvec(cellfun( @(a,b,c)( plot(axes_hdl, rowvec(a(b(c))), repmat(titleYposition - (0.1 * plotYrange), 1,length(c)) )), qrs_ploted, qrs2plot, aux_seq, 'UniformOutput', false ));
%                 set_rand_linespec(QRSfpFarHdls, 'v', 'none', [], 6 );
                cellfun( @(a,b)(set_a_linespec(a, b)), QRSfpFarHdls, cLinespecsNone(1:length(QRSfpHdls)) );

            end
        else
            if( bAux )
                str_aux = 'on';
            else
                str_aux = 'off';
            end
            
            cellfun(@(a)(set(a, 'Visible', str_aux )), QRSfpFarHdls);
            
        end

        bAux = eDetailLevel ~= kNoDetail && ( ( eDetailLevel == kCloseDetailML || eDetailLevel == kCloseDetailAll ) && plotXrange < mediumDetailSampSize);
        if( isempty(QRSfpHdls) ) 
            bAux2 = any( colvec(cellfun( @(a,b,c)(~isempty(a(b(c)))), qrs_ploted, qrs2plot, aux_seq )) );
            if( bAux && bAux2 )
                QRSfpHdls = colvec(cellfun( @(a,b,c)( plot(axes_hdl, repmat(rowvec(a(b(c))), 2,1), [ repmat(plotYmin + (0.06 * plotYrange), 1,length(c)); repmat(titleYposition - (0.1 * plotYrange), 1,length(c)) ] ) ), qrs_ploted, qrs2plot, aux_seq , 'UniformOutput', false));
%                 aux_val = cellfun( @(a)(set_rand_linespec(a(1), '^', ':', [], 5 )), QRSfpHdls, 'UniformOutput', false);
                cellfun( @(a,b)(set_a_linespec(a, b)), QRSfpHdls, cLinespecs(1:length(QRSfpHdls)) );
                cellfun( @(a)(set(a, 'LineWidth', 0.25)), QRSfpHdls);
                hb_classifier_idx = find(strcmpi(qrs_ploted_names, 'hb_classifier'));
                if( ~isempty(hb_classifier_idx) )
                    
                    aux_hb_classifier_val = colvec(qrs_ploted{hb_classifier_idx});
                    aux_seq_classifier = aux_seq{hb_classifier_idx};
                    qrs2plot_classifier = qrs2plot{hb_classifier_idx};
                    QRSfpHdls = [ QRSfpHdls; ...
                                arrayfun( @(b)( ...
                                text( aux_hb_classifier_val(b), titleYposition - (0.08 * plotYrange), cBeatLabels{hb_labels_idx(b)}, ...
                                    'FontSize', 8, ...
                                    'VerticalAlignment', 'middle', ...
                                    'HorizontalAlignment', 'left', ...
                                    'Interpreter', 'none', ...
                                    'Color', 1-cBeatLabelsColorCode{hb_labels_idx(b)}, ... 
                                    'BackgroundColor', cBeatLabelsColorCode{hb_labels_idx(b)}, ... 
                                    'EdgeColor', cBeatLabelsColorCode{hb_labels_idx(b)} ) ...
                                ), qrs2plot_classifier(aux_seq_classifier), 'UniformOutput', false ) ];
                end
                
            end
        else
            if( bAux )
                str_aux = 'on';
            else
                str_aux = 'off';
            end
            cellfun(@(a)(set(a, 'Visible', str_aux )), QRSfpHdls);
        end

        bAux = eDetailLevel ~= kNoDetail && ( ( eDetailLevel == kCloseDetailML || eDetailLevel == kCloseDetailAll ) && plotXrange <= closeDetailSampSize );
        if( isempty(QRScplxHdls) )
            if( bAux )
%                 QRScplxHdls = [ QRScplxHdls; PlotGlobalWaveMarks({'QRSon' 'qrs' 'QRSoff'}, [ plotYmin + (0.08*plotYrange)  titleYposition - (0.1*plotYrange) ], color_QRScomplex_global)];
                for jj = ECG_signals_idx
                    QRScplxHdls = [ QRScplxHdls; PlotWaveMarks(global_annotations, {'QRSon' 'qrs' 'QRSoff'}, jj, 2*yTextOffset, QRScplxColor*0.95 ) ];
                end
            end
        else
            if( bAux )
                str_aux = 'on';
            else
                str_aux = 'off';
            end
            cellfun(@(a)(set(a, 'Visible', str_aux )), QRScplxHdls);
            
        end

        % T wave
        bAux = eDetailLevel ~= kNoDetail && ( ( eDetailLevel == kCloseDetailML || eDetailLevel == kCloseDetailAll ) && plotXrange <= closeDetailSampSize );
        if( isempty(TwaveHdls) )
            if( bAux )
%                 TwaveHdls = [ TwaveHdls; PlotGlobalWaveMarks({'Ton' 'T' 'Toff'}, [ plotYmin + (0.09*plotYrange)  plotYmax - (0.07*plotYrange) ], color_Twave_global)];
                for jj = ECG_signals_idx
                    TwaveHdls = [ TwaveHdls; PlotWaveMarks(this_annotation, {'Ton' 'T' 'Toff'}, jj, yTextOffset, TwaveColor*0.95 ) ];
                end
            end
        else
            if( bAux )
                str_aux = 'on';
            else
                str_aux = 'off';
            end
            cellfun(@(a)(set(a, 'Visible', str_aux )), TwaveHdls);
            
        end    

    end

    bAux = eDetailLevel ~= kNoDetail && ( (eDetailLevel == kCloseDetailML || eDetailLevel == kCloseDetailAll ) && plotXrange <= mediumDetailSampSize );
    
    % QRS class labels
    if( bAux && ~isempty(hb_labels_idx) )

        if(bWaveLegendPlotted)
            left_hb_legend = left_legend - 1.3*width_legend;
        else
            left_hb_legend = left_legend - 0.05*width_legend;
        end

        width_hb_legend = 1.1*width_legend;
        height_hb_legend = 1.1*height_legend;
        bottom_hb_legend = bottom_legend - 0.05*height_legend;

        bHBclassLegendPlotted = true;

        % Heartbeats class frame
        legend_hdl = {patch([left_hb_legend left_hb_legend [left_hb_legend left_hb_legend]+width_hb_legend left_hb_legend ], [bottom_hb_legend [bottom_hb_legend bottom_hb_legend]+height_hb_legend bottom_hb_legend bottom_hb_legend], [1 1 1], 'EdgeColor', [0 0 0])};
        lcBeatLabels = length(cBeatLabels);
        for jj = 1:lcBeatLabels
            legend_hdl = [legend_hdl; {text( left_hb_legend + jj*width_hb_legend/(lcBeatLabels+1), bottom_hb_legend + height_hb_legend/2, cBeatLabels{jj}, 'FontSize', 8, 'HorizontalAlignment', 'center', 'BackgroundColor', QRScplxColor, 'EdgeColor', cBeatLabelsColorCode{jj} )} ];
        end
        UserChnageViewHdls = [UserChnageViewHdls; legend_hdl];
    end

%     bAux = eDetailLevel ~= kNoDetail && ( (eDetailLevel == kCloseDetailML || eDetailLevel == kCloseDetailAll ) && plotXrange <= closeDetailSampSize );
    bAux = eDetailLevel ~= kNoDetail && ( (eDetailLevel == kCloseDetailML || eDetailLevel == kCloseDetailAll ) );
    
    % QRS annotations labels
    bAux2 = any( colvec(cellfun( @(a,b,c)(~isempty(a(b(c)))), qrs_ploted, qrs2plot, aux_seq )) );
    
    if( bAux && bAux2 && ~isempty(qrs_ploted_names) )

        cant_qrs_names = length(qrs_ploted_names);

        if(bWaveLegendPlotted)
            left_qrs_names_legend = left_legend;
            width_qrs_names_legend = width_legend;
            height_qrs_names_legend = height_legend * cant_qrs_names;
            bottom_qrs_names_legend = bottom_legend - 0.5*height_legend - height_qrs_names_legend;
        elseif(bHBclassLegendPlotted)
            left_qrs_names_legend = left_hb_legend;
            width_qrs_names_legend = width_hb_legend;
            height_qrs_names_legend = height_hb_legend * cant_qrs_names;
            bottom_qrs_names_legend = bottom_hb_legend - 0.5*height_hb_legend - height_qrs_names_legend;
        else
            left_qrs_names_legend = left_legend ;
            width_qrs_names_legend = width_legend;
            height_qrs_names_legend = height_legend * cant_qrs_names;
            bottom_qrs_names_legend = bottom_legend + height_legend - height_qrs_names_legend;
        end

        % QRS names frame
        legend_hdl = {patch([left_qrs_names_legend left_qrs_names_legend [left_qrs_names_legend left_qrs_names_legend]+width_qrs_names_legend left_qrs_names_legend ], [bottom_qrs_names_legend [bottom_qrs_names_legend bottom_qrs_names_legend]+height_qrs_names_legend bottom_qrs_names_legend bottom_qrs_names_legend], [1 1 1], 'EdgeColor', [0 0 0])};
        for jj = 1:cant_qrs_names
            aux_ls = cLinespecs{jj};
            legend_hdl = [legend_hdl; {text( left_qrs_names_legend + 0.1*width_qrs_names_legend, bottom_qrs_names_legend + height_qrs_names_legend * (jj)/(cant_qrs_names+1), adjust_string( qrs_ploted_names{jj}, 15 ), ...
                                            'FontSize', 8, ...
                                            'VerticalAlignment', 'middle', ...
                                            'HorizontalAlignment', 'left', ...
                                            'Interpreter', 'none', ...
                                            'Color', 1-aux_ls{4}, ...
                                            'BackgroundColor', aux_ls{4}, ... 
                                            'EdgeColor', aux_ls{5} )} ... 
                            ];
        end
        UserChnageViewHdls = [UserChnageViewHdls; legend_hdl];
    end    

    %% single-lead annotations
    
    if( ~isempty(annotations) )

        if( eDetailLevel ~= kNoDetail && ( (eDetailLevel == kCloseDetailSL || eDetailLevel == kCloseDetailAll ) && plotXrange <= closeDetailSampSize ) )
            bWaveLegendPlotted = true;
            % frame
            legend_hdl = {patch([left_legend left_legend [left_legend left_legend]+width_legend left_legend ], [bottom_legend [bottom_legend bottom_legend]+height_legend bottom_legend bottom_legend], [1 1 1], 'EdgeColor', [0 0 0])};
            legend_hdl = [legend_hdl; {text( left_legend + width_legend/4, bottom_legend + height_legend/2, 'P', 'FontSize', 8, 'HorizontalAlignment', 'center', 'BackgroundColor', PwaveColor)}];
            legend_hdl = [legend_hdl; {text( left_legend + width_legend/2, bottom_legend + height_legend/2, 'QRS', 'FontSize', 8, 'HorizontalAlignment', 'center', 'BackgroundColor', QRScplxColor )}];
            legend_hdl = [legend_hdl; {text( left_legend + width_legend*3/4, bottom_legend + height_legend/2, 'T', 'FontSize', 8, 'HorizontalAlignment', 'center', 'BackgroundColor', TwaveColor )}];
            UserChnageViewHdls = [UserChnageViewHdls; legend_hdl];
        end
        
        for jj = 1:length(ECG_signals_idx)

            this_annotation = annotations(jj);
            
            % P wave
            bAux = ( eDetailLevel ~= kNoDetail && ( ( eDetailLevel == kCloseDetailSL || eDetailLevel == kCloseDetailAll ) && plotXrange <= closeDetailSampSize ) );
            if( isempty(PwaveGlblHdls{jj}) )
                if( bAux )
%                     PwaveGlblHdls{jj} = [ PwaveGlblHdls{jj}; PlotWaveMarks(this_annotation, {'Pon' 'P' 'Poff'}, jj, 0.5*yTextOffset, ColorOrder(jj,:) )];
                    PwaveGlblHdls{jj} = [ PwaveGlblHdls{jj}; PlotWaveMarks(this_annotation, {'Pon' 'P' 'Poff'}, ECG_signals_idx(jj), 0.5*yTextOffset, PwaveColor )];
                end
            else
                if( bAux )
                    str_aux = 'on';
                else
                    str_aux = 'off';
                end
                cellfun(@(a)(set(a, 'Visible', str_aux )), PwaveGlblHdls{jj});
                
            end
            
            % QRS complex

            bAux = eDetailLevel ~= kNoDetail && ( ( eDetailLevel == kCloseDetailSL || eDetailLevel == kCloseDetailAll ) && plotXrange <= closeDetailSampSize);
            if( isempty(QRScplxGlblHdls{jj}) )
                if( bAux )
                    QRScplxGlblHdls{jj} = [ QRScplxGlblHdls{jj}; PlotWaveMarks(this_annotation, {'QRSon' 'qrs' 'QRSoff'}, ECG_signals_idx(jj), 2*yTextOffset, QRScplxColor)];
                end
            else
                if( bAux )
                    str_aux = 'on';
                else
                    str_aux = 'off';
                end
                cellfun(@(a)(set(a, 'Visible', str_aux )), QRScplxGlblHdls{jj});
                
            end

            % T wave
            bAux = eDetailLevel ~= kNoDetail && ( ( eDetailLevel == kCloseDetailSL || eDetailLevel == kCloseDetailAll ) && plotXrange <= closeDetailSampSize );
            if( isempty(TwaveGlblHdls{jj}) )
                if( bAux )
%                     TwaveGlblHdls{jj} = [ TwaveGlblHdls{jj}; PlotWaveMarks(this_annotation, {'Ton' 'T' 'Toff'}, jj, yTextOffset, ColorOrder(jj,:) )];
                    TwaveGlblHdls{jj} = [ TwaveGlblHdls{jj}; PlotWaveMarks(this_annotation, {'Ton' 'T' 'Toff'}, ECG_signals_idx(jj), yTextOffset, TwaveColor )];
                end
            else
                if( bAux )
                    str_aux = 'on';
                else
                    str_aux = 'off';
                end
                cellfun(@(a)(set(a, 'Visible', str_aux )), TwaveGlblHdls{jj});
                
            end    
        end    
    end

    %% Time/Voltage scales

    % ECG voltage
    if( bPlotScales )
        
        PlotECGscale(limits);

        % time
        rect_scale_height = 4*yTextOffset;
        rect_scale_width = 0.06 * plotXrange;
        rect_scale_X = plotXmax - (0.07 * plotXrange);
        rect_scale_Y = plotYmin + (0.07 * plotYrange);
        UserChnageViewHdls = [UserChnageViewHdls; {rectangle('Position',[rect_scale_X, rect_scale_Y, rect_scale_width, rect_scale_height], 'FaceColor', [1 1 1])}];

        XscaleSize = floor( 0.03 * plotXrange * 1/heasig.freq) * heasig.freq ; 
        if( XscaleSize == 0 )
            % multiple of mult milliseconds
            mult = 800;
            XscaleSize_ms = 0;
            while( XscaleSize_ms == 0 && mult > 1)
                XscaleSize_ms =  floor(floor( 0.03 * plotXrange * 1e3 * 1/heasig.freq) / mult) * mult;
                mult = round(mult / 2);
            end
            XscaleSize = XscaleSize_ms * heasig.freq / 1e3;
            UserChnageViewHdls = [ UserChnageViewHdls; {text( rect_scale_X + round(rect_scale_width/2), rect_scale_Y + yTextOffset, [ num2str(XscaleSize_ms) ' ms' ], 'FontSize', 8, 'HorizontalAlignment', 'center')}];
        else
            % multiple of mult seconds
            mult = 60;
            XscaleSize_s = 0;
            while( XscaleSize_s == 0)
                XscaleSize_s =  floor(floor( 0.03 * plotXrange * 1/heasig.freq) / mult) * mult;
                mult = round(mult / 2);
            end
            XscaleSize = XscaleSize_s * heasig.freq;
            UserChnageViewHdls = [ UserChnageViewHdls; {text( rect_scale_X + round(rect_scale_width/2), rect_scale_Y + yTextOffset, [ num2str(XscaleSize_s) ' s' ], 'FontSize', 8, 'HorizontalAlignment', 'center')}];
        end
        UserChnageViewHdls = [UserChnageViewHdls; {plot(axes_hdl, [ -(XscaleSize/2) -(XscaleSize/2) XscaleSize/2 XscaleSize/2] + rect_scale_X + round(rect_scale_width/2) , [ -yTextOffset 0 0 -yTextOffset] + rect_scale_Y + 3*yTextOffset, 'k-',  'LineWidth', 0.25 )}];
    end

    %% Title
    
    if( bPaperModeOn )
        UserChnageViewHdls = [UserChnageViewHdls; ...
                { text( plotXmin + 0.5 * plotXrange, titleYposition, ... 
                        strTitle, ...
                        'BackGroundColor', [0.99 0.92 0.8], ...
                        'EdgeColor', [1 0 0], ...
                        'LineWidth', 1.2, ...
                        'Interpreter', 'none', ...
                        'FontSize', 9, ...
                        'HorizontalAlignment', 'center' )}; ];
    else                
        UserChnageViewHdls = [UserChnageViewHdls; ...
                { text( plotXmin + 0.5 * plotXrange, titleYposition, ... 
                        strTitle, ...
                        'BackGroundColor', [0.702 0.78 1], ...
                        'EdgeColor', [0.078 0.169 0.549], ...
                        'LineWidth', 1.2, ...
                        'Interpreter', 'none', ...
                        'FontSize', 9, ...
                        'HorizontalAlignment', 'center' )}; ];
    end
    
    titleExtent = get(UserChnageViewHdls{end}, 'Extent');

    topLevelHdl = [topLevelHdl; UserChnageViewHdls(end)];
    

    end

%--------------------------------------------------------------------------

%==========================================================================

    function aux_seq = limit_sequence(qrs2plot, limits)

        aux_QRS2plot = 15;
        
        aux_idx = find( qrs2plot > limits(1) & qrs2plot < limits(2) );
        
        if( isempty(aux_idx) )
            aux_seq = [];
        else
            if( aux_QRS2plot < length(aux_idx) )
                aux_seq = unique(round(linspace( 1, length(aux_idx), aux_QRS2plot)));
                aux_seq = aux_idx(aux_seq);
            else
                aux_seq = aux_idx;
            end
        end
        
    end

%==========================================================================

%--------------------------------------------------------------------------

%==========================================================================

    function PlotECGscale(limits)

        plotXmin = limits(1);
        plotXmax = limits(2);
        plotYmin = limits(3);
        plotYmax = limits(4);

        plotXrange = plotXmax - plotXmin;
        plotYrange = plotYmax - plotYmin;

        xTextOffset = 0.005*plotXrange;

        rect_scale_height = 0.13 * plotYrange;
        rect_scale_width = 4*xTextOffset;
        rect_scale_X = plotXmax - 6*xTextOffset;
        rect_scale_Y = plotYmin + (0.13 * plotYrange);

        UserChnageViewHdls = [UserChnageViewHdls; {rectangle('Position',[rect_scale_X, rect_scale_Y, rect_scale_width, rect_scale_height], 'FaceColor', [1 1 1])}];

        YscaleSize = rect_scale_height / gains(volt_idx(1)); 
        if( YscaleSize < 1e-3 )
            % picovolts
            k = 1e-6;
            str_aux = 'p';
        elseif( YscaleSize < 1 && YscaleSize >= 1e-3  )
            % nanovolts
            k = 1e-3;
            str_aux = 'n';
        elseif( YscaleSize >= 1 && YscaleSize <= 1e3)
            % microvolts
            k = 1;
            str_aux = '{\mu}';
        elseif( YscaleSize < 1e6 && YscaleSize >= 1e3  )
            % milivolts
            k = 1e3;
            str_aux = 'm';
        elseif( YscaleSize >= 1e6  )
            % volts
            k = 1e6;
            str_aux = '';
        end

        mult = 800;
        YscaleSize = 0;
        while( YscaleSize == 0 && mult >= 1)
            YscaleSize =  floor(floor( rect_scale_height / lead_gain / k ) / mult ) * mult;
            mult = round(mult / 2);
        end
        YscaleSize_uv = YscaleSize * k * lead_gain;
        UserChnageViewHdls = [ UserChnageViewHdls; {text( plotXmax - 5*xTextOffset, rect_scale_Y + rect_scale_height/2, [ num2str(YscaleSize) ' ' str_aux 'V' ], 'FontSize', 8, 'HorizontalAlignment', 'center', 'Rotation', 90)}];

        UserChnageViewHdls = [UserChnageViewHdls; {plot(axes_hdl, [ -xTextOffset 0 0 -xTextOffset] + (plotXmax - 3*xTextOffset) , [ -YscaleSize_uv/2 -YscaleSize_uv/2 YscaleSize_uv/2 YscaleSize_uv/2 ] + rect_scale_Y + rect_scale_height/2, 'k-', 'LineWidth', 0.25 )}];

    end

%--------------------------------------------------------------------------

%==========================================================================

    function UserChangeView(obj, event_obj, id_caller)

        if( nargin < 3 )
            id_caller = '';
        end

        hold(axes_hdl, 'on')

        xlimits = get(axes_hdl,'Xlim'); 
        ylimits = get(axes_hdl,'Ylim');

        plotYmax = ylimits(2); 
        plotYmin = ylimits(1);
        plotXmax = xlimits(2); 
        plotXmin = xlimits(1);

        if( strcmpi(id_caller, 'pan') )

            if( plotXmin < original_plot_lims(1) )
                plotXmin = original_plot_lims(1);
                plotXmax = original_plot_lims(1) + prev_Xrange;
            end
            if( plotXmax > original_plot_lims(2) )
                plotXmax = original_plot_lims(2);
                plotXmin = original_plot_lims(2) - prev_Xrange;
            end
            aux_lim = (original_plot_lims(3) * max(gains)) - max(offsets);
            if( plotYmin < aux_lim )
                plotYmin = aux_lim;
                plotYmax = aux_lim + prev_Yrange;
            end
            aux_lim = (original_plot_lims(4) * min(gains)) - min(offsets);
            if( plotYmax > aux_lim )
                plotYmax = aux_lim;
                plotYmin = aux_lim - prev_Yrange;
            end

        elseif( strcmpi(id_caller, 'zoom') )
%             plotYrange = diff(ylimits);

            % additional space for x axis
%             plotYmin = plotYmin - 0.04*plotYrange;

            plotYmax = min(plotYmax, original_plot_lims(4));
            plotYmin = max(plotYmin, original_plot_lims(3));
            % plotYmin = min(plotYmin, user_data.original_plot_lims(4) - 0.2*plotYrange);
            % plotYmax = max(plotYmax, user_data.original_plot_lims(3) + 0.2*plotYrange);

        end

        plotYrange = plotYmax - plotYmin;
        plotXrange = plotXmax - plotXmin;

        prev_Yrange = plotYrange;
        prev_Xrange = plotXrange;

        if( plotYmin < plotYmax )
            set(axes_hdl, 'Ylim', [plotYmin plotYmax]); 
        end

        aux_qrs_ploted = cellfun( @(a)( find(a >= plotXmin & a <= plotXmax) ), qrs_ploted, 'UniformOutput', false );

        % additional space for x axis
%         Xmargin = plotXrange * 0.1;
%         plotXmin = max(plotXmin - Xmargin, original_plot_lims(1));
%         plotXmax = min(plotXmax + Xmargin, original_plot_lims(2));
%         plotXrange = plotXmax - plotXmin;

        set(axes_hdl, 'Xlim', [ plotXmin plotXmax]);

        cellfun(@(a)( CheckAndDeleteHdl(a) ), UserChnageViewHdls);
        UserChnageViewHdls = [];
%         if( ~isempty(UserChnageViewHdls) )
%             if( ishandle(UserChnageViewHdls) )
%                 delete(UserChnageViewHdls);
%             end
%             UserChnageViewHdls = [];
%         end

        update_axis_plot_ecg_strip( aux_qrs_ploted );

        if( ishandle(ECG_hdl) )
            uistack(ECG_hdl, 'top');
        end
        if( ishandle(ECGd_hdl) )
            uistack(ECGd_hdl, 'top');
        end
        
        hold(axes_hdl, 'off')

%         set(obj, 'UserData', user_data);

        if( ~strcmpi(event_obj, 'JustUpdate') )
            %save focus
            aux_fig = gcf;
            % update linked plots if any
            for jj = rowvec(linked_hdl)
                figure(jj)
                axes_hdl = get(jj,'CurrentAxes');
                set(axes_hdl, 'Xlim', xlimits);
                UserChangeView(jj, 'JustUpdate', id_caller)
            end
            %restore focus
            figure(aux_fig)

        end

    end

%--------------------------------------------------------------------------

%==========================================================================

    function [ecg_range ecg_min ecg_max ecg_median] = CalcECG_range(ECG)
        
    ecg_prctiles = cell2mat(arrayfun( @(a)(prctile( randsample(ECG(:,a), max( 100, round(cant_samp * cant_leads / 100) ) ), [ 2.5 50 97.5 ] )), (1:cant_leads)', 'UniformOutput', false) )';
    ecg_range = ecg_prctiles(3,:) - ecg_prctiles(1,:);
    ecg_max =  colvec(ecg_prctiles(3,:) + 0.7 * ecg_range);
    ecg_min =  colvec(ecg_prctiles(1,:) - 0.4 * ecg_range);
    ecg_median =  colvec(ecg_prctiles(2,:));
    ecg_range = ecg_max - ecg_min;

end

%--------------------------------------------------------------------------

%==========================================================================

    function [plotYrange plotYmin plotYmax] = CalcPlotYlimits(ecg_min, ecg_max, lead_gain, lead_offset)

%         plotYmax = max(ecg_max*(lead_gain) - lead_offset); 
%         plotYmin = min(ecg_min*(lead_gain) - lead_offset);
        plotYmax = max(colvec(ecg_max).*lead_gain - colvec(lead_offset) ) ; 
        plotYmin = min(colvec(ecg_min).*lead_gain - colvec(lead_offset) );
        plotYrange = plotYmax - plotYmin;
        % additional space for x axis
        
        if( plotYrange == 0 )
            plotYmax = plotYmin + 0.1; 
            plotYmin = plotYmin - 0.1;
            plotYrange = plotYmax - plotYmin;
        end
        
        plotYmin = plotYmin - plot_bottom_margin_width*plotYrange;
        plotYmax = plotYmax + plot_top_margin_width*plotYrange;

    end

%--------------------------------------------------------------------------

%==========================================================================
    function WindowButtonDownCallback2D(src, evnt)    %#ok
        %WindowButtonDownCallback2D
        
        if fIsMouseOnLegend, return; end
        
        clickType = get(src, 'SelectionType');
        
        switch clickType
            case 'normal'
                DragMouseBegin();
                mMgDirection = 'plus';
            case 'open'
                if fIsMagnifierOn
                    MagnifierSizeChange(mMgDirection);
                else
                    ResetAxesToOrigView();
                end
            case 'alt'
                RubberBandBegin();
                mMgDirection = 'minus';
            case 'extend'
                if fIsMagnifierOn
                    MagnifierReset();
                else
                    
                position_mouse = get(axes_hdl, 'CurrentPoint');
                prev_lead_idx = lead_selected_idx;
                lead_selected_idx = max(1,heasig.nsig) * abs(ecg_max(1) * gains(1) - position_mouse(1,2)) / (ecg_max(1) * gains(1) - (-offsets(cant_leads) + ecg_min(cant_leads)*gains(cant_leads)) );
%                 update_title_efimero( num2str(lead_selected_idx), 5 );
                
                lead_selected_idx = min(cant_leads, max(1, 1+floor( lead_selected_idx ) ) );
                
                if( length(prev_lead_idx) == length(lead_selected_idx) )
                    % disable lead selection
                    lead_selected_idx = 1:cant_leads;
                    update_title_efimero( 'All leads', 5 );
                else
                    update_title_efimero( heasig.desc(lead_selected_idx,:), 5 );
                end
                
                    
                    %just measure on graph
                    fNoZoom = true;
                    RubberBandBegin();
%                     ZoomMouseExtendBegin();
                end
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function WindowButtonUpCallback2D(src, evnt)      %#ok
        %WindowButtonUpCallback2D
        
        DragMouseEnd();
%         ZoomMouseExtendEnd();
        RubberBandEnd();
        fNoZoom = false;
    end
%--------------------------------------------------------------------------

%==========================================================================
    function WindowButtonMotionCallback2D(src, evnt)  %#ok
        %WindowButtonMotionCallback2D
                
        if ~(fIsMagnifierOn || fIsDragAllowed || fIsRubberBandOn)
            % set current axes under cursor
            SelectAxesUnderCursor();
        end
        
        if fIsEnableControl
            DragMouse();
            RubberBandUpdate();
            MagnifierUpdate();
        end
        
%         ZoomMouseExtend();
        PointerCrossUpdate();
    end
%--------------------------------------------------------------------------

%==========================================================================
    function WindowScrollWheelFcn2D(src, evnt)        %#ok
        %WindowScrollWheelFcn2D
        
        if fIsMouseOnLegend, return; end
        
        % Update Zoom Info
        % because it can be changed function 'zoom'
        UpdateCurrentZoomAxes();
        
        switch mZoomScroll
            case 'normal'
                directions = {'minus', 'plus'};
            case 'reverse'
                directions = {'plus', 'minus'};
        end
        
        verScrollCount = evnt.VerticalScrollCount;
        
        if (verScrollCount > 0)
            direction = directions{1};
        elseif (verScrollCount < 0)
            direction = directions{2};
        else
            return;
        end
        
        if( strcmpi(scroll_mode, 'gain/offset') )
            
            changeGainOffset(verScrollCount);
            set(src, 'UserData', user_data);
            
        elseif( strcmpi(scroll_mode, 'zoom') )
            ZoomMouse(direction);
            PointerCrossUpdate();
            MagnifierZoomChange(direction);
            UserChangeView( fig_hdl, [], 'zoom');
        end
        
    end
%--------------------------------------------------------------------------

    function changeGainOffset(verScrollCount)
        
        prev_offset = lead_offset;
        prev_gain = gains;

        bLeadMask = false(heasig.nsig,1);
        bLeadMask(lead_selected_idx) = true;
        
        if( strcmp(gain_offset_mode, 'offset' ) )
            % offset
            k_offset = verScrollCount * ecg_range .* gains * 0.05;

            this_lead_offset = lead_offset + k_offset;
%             this_lead_offset = max( 0, lead_offset + k_offset);

            if( all(prev_offset == this_lead_offset) )
                return
            else
                lead_offset = this_lead_offset ;
            end
            
            % retain all not selected
            lead_offset(~bLeadMask) = prev_offset(~bLeadMask);
            
%             offsets = abs(cumsum(lead_offset));
            offsets = cumsum(lead_offset);
            
        else
            % gain
            this_lead_gain = 1.5^(-verScrollCount) * gains;

            bAux = this_lead_gain < min_gain | this_lead_gain > max_gain;
                
            gains = this_lead_gain;
            gains(bAux) = prev_gain(bAux);

            % retain all not selected
            gains(~bLeadMask) = prev_gain(~bLeadMask);
            
        end


        [plotYrange plotYmin plotYmax] = CalcPlotYlimits(ecg_min, ecg_max, gains, offsets);

        cla(axes_hdl);
        set(axes_hdl, 'ColorOrder', ColorOrder);

        hold(axes_hdl, 'on')
        %plot ECG
        ECG_hdl = plot(axes_hdl, start_sample:end_sample, bsxfun( @minus, bsxfun( @times, ECG, rowvec(gains)), rowvec(offsets) ), 'LineWidth', 1.3 );

        set(axes_hdl, 'Box', 'off' );
        set(axes_hdl, 'Xtick', [] );
        set(axes_hdl, 'Ytick', [] );
        set(axes_hdl, 'Xcolor', [1 1 1] );
        set(axes_hdl, 'Ycolor', [1 1 1] );

        ylim([plotYmin plotYmax]);

        hold(axes_hdl, 'off')

        PwaveHdls = [];
        TwaveHdls = [];
        QRScplxHdls = [];
        PwaveGlblHdls = cell(heasig.nsig,1);
        TwaveGlblHdls = cell(heasig.nsig,1);
        QRScplxGlblHdls = cell(heasig.nsig,1);
        QRSfpHdls = [];
        QRSfpFarHdls = [];
        
        if( strcmp(gain_offset_mode, 'offset' ) )
            update_title_efimero( num2str(rowvec(offsets)), 5 );
        else
            update_title_efimero( num2str(rowvec(gains)), 5 );
        end
    end


%==========================================================================
    function WindowKeyPressCallback2D(src, evnt)      %#ok
        %WindowKeyPressCallback2D
        
        modifier = evnt.Modifier;
        
        switch evnt.Key
            case '0'
                ResetAxesToOrigView();
            case {'equal', 'add'}
                ZoomKeys('plus');
            case {'hyphen', 'subtract'}
                ZoomKeys('minus');
            
            case 'a'
                
                switch( ann_graph_mode )

                    case kBackColourAnns
                        ann_graph_mode = kLinesAnns;
                        update_title_efimero( 'Annotations line mode', 5 );
                    
                    case kLinesAnns
                        ann_graph_mode = kBackColourAnns;
                        update_title_efimero( 'Annotations colour mode', 5 );
                    
                    otherwise
                        ann_graph_mode = kBackColourAnns;
                        update_title_efimero( 'Annotations colour mode', 5 );
                        
                end

                cellfun( @(a)(delete(a)), [PwaveGlblHdls; TwaveGlblHdls; QRScplxGlblHdls]);
                
                PwaveGlblHdls = cell(heasig.nsig,1);
                TwaveGlblHdls = cell(heasig.nsig,1);
                QRScplxGlblHdls = cell(heasig.nsig,1);
                
                UserChangeView( fig_hdl, [], 'pan');                
                
            case 'h'
                disp_help();
                
            case 'c'
                SetPointerCrossKeys();
                
            case 'd'

                switch( eDetailLevel )

                    case kNoDetail
                        eDetailLevel = kCloseDetailSL;
                        update_title_efimero( 'Close SL', 5 );
                    
                    case kCloseDetailSL
                        eDetailLevel = kMediumDetailSL;
                        update_title_efimero( 'Medium SL', 5 );
                    
                    case kMediumDetailSL
                        eDetailLevel = kCloseDetailML;
                        update_title_efimero( 'Close ML', 5 );

                    case kCloseDetailML
                        eDetailLevel = kMediumDetailML;
                        update_title_efimero( 'Medium ML', 5 );

                    case kMediumDetailML
                        eDetailLevel = kCloseDetailAll;
                        update_title_efimero( 'Close All', 5 );
                        
                    case kCloseDetailAll
                        eDetailLevel = kMediumDetailAll;
                        update_title_efimero( 'Medium All', 5 );

                    case kMediumDetailAll
                        eDetailLevel = kNoDetail;
                        update_title_efimero( 'No detail', 5 );
                        
                    otherwise
                        eDetailLevel = kNoDetail;
                        update_title_efimero( 'No detail', 5 );
                end
                
                UserChangeView( fig_hdl, [], 'pan');                
                
            case 'p'
                
                if( bPaperModeOn )
                    PaperModeOff();
                else
                    PaperModeOn();
                end
                
            case 'g'
                scroll_mode = 'gain/offset';
                gain_offset_mode = 'gain';
%                 update_title_efimero( 'Gain', inf );
                
            case 'o'
                scroll_mode = 'gain/offset';
                gain_offset_mode = 'offset';
%                 update_title_efimero( 'Offset', inf );
                
            case 'x'
                fIsEnableZoomX = ~fIsEnableZoomX;
                if( fIsEnableZoomX )
                    fIsEnableZoomY = false;
                    update_title_efimero( 'X mode', 10 );
                else
                    update_title_efimero( 'XY mode', 10 );
                end
                
                
            case 'y'
                fIsEnableZoomY = ~fIsEnableZoomY;
                if( fIsEnableZoomY )
                    fIsEnableZoomX = false;
                    update_title_efimero( 'Y mode', 10 );
                else
                    update_title_efimero( 'XY mode', 10 );
                end
                
            case 'm'
                if fIsEnableControl
                    MagnifierOn();
                    update_title_efimero( 'Magnifier', 5 );
                end
                
            case 'r'
                
                report_format_idx = report_format_idx + 1;
                if(report_format_idx > length(cKnownReportFormats) || report_format_idx < 1 )
                    report_format_idx = 1;
                end
                report_format = cKnownReportFormats{report_format_idx};
                
                update_title_efimero( report_format, 5 );
                
            case 's'
                
                SaveReport();

            case 'b'
                
                if( bRelative_begin )
                    bRelative_begin = false;
                    base_time = prev_base_time;
                    update_title_efimero( 'Absolute time', 5 );
                else
                    bRelative_begin = true;
                    prev_base_time = base_time;
                    base_time = 1;
                    update_title_efimero( 'Relative from start time', 5 );
                end
                
        end
    end
%--------------------------------------------------------------------------


    function SaveReport()

        if( isempty(report_filename) ) 

            if( isempty(ECG_w) )
                report_path = [pwd filesep];
            else
                report_path = fileparts(ECG_w.recording_name);
                report_path = [report_path filesep];
            end

            if( isfield(heasig, 'recname') )
                report_filename = [report_path heasig.recname '_captured_' datestr(now, 'dd_mm_yy-HH_MM_SS' )  '.' report_format];
            else
                report_filename = [report_path 'ECG_strip_captured_' datestr(now, 'dd_mm_yy-HH_MM_SS' ) '.' report_format];
            end

        else

            report_path = fileparts(report_filename);
            report_path = [report_path filesep];

        end

        if( exist(report_path, 'dir') )

            init_ghostscript();

            export_fig(report_filename, '-nocrop', ['-' report_format], fig_hdl);

            update_title_efimero( ['Exported to ' report_filename], 5 );

        else
            if( bCalledFromWorkspace )
                cprintf('[1,0.5,0]', 'Could not create report file: folder %s does not exist\n', report_path );
            end
        end        
        
    end

%==========================================================================
    function PaperModeOff()
        bPaperModeOn = false;
        bAux = cellfun(@(a)(ishandle(a)), paperModeHdl );
        if( any(bAux) )
            cellfun(@(a)(delete(a)), paperModeHdl(bAux) );
            paperModeHdl = {};
        end
        if( ishandle(ECGd_hdl) )
            delete(ECGd_hdl);
        end
        
        set(ECG_hdl, 'Visible', 'on')
        
        UserChangeView( fig_hdl, [], 'zoom');
    end

%--------------------------------------------------------------------------

%==========================================================================
    function PaperModeOn()

        % definitions
        top_frame = plotYmax - 2* yTextOffset ;
        bottom_frame = plotYmin + 2* yTextOffset ;
        left_frame = startSignalX;
        right_frame = endSignalX;
        
        % some updates in the view
        bPaperModeOn = true;
        
        set(ECG_hdl, 'Visible', 'off')
        
        hold(axes_hdl, 'on');
        ECGd_hdl = plot(axes_hdl, round(linspace(start_sample, end_sample, size(ECGd,1))), bsxfun( @minus, bsxfun( @times, ECGd, rowvec(gains)), rowvec(offsets) ), 'LineWidth', 1.3 );
        hold(axes_hdl, 'off');
        
        UserChangeView( fig_hdl, [], 'zoom');
        
        % time grid
        [~, major_tick_idx] = sort( abs((plotXrange./major_tick_values_time) - 10));
        major_tick = major_tick_values_time(major_tick_idx(1));
        
        this_start = ceil(max(0,plotXmin)/major_tick)*major_tick;
        this_end = floor(min(heasig.nsamp,plotXmax)/major_tick)*major_tick;
        major_tick_x = this_start:major_tick:this_end;
        
        minor_tick = major_tick / 5;
        this_start = ceil(max(1,plotXmin)/minor_tick)*minor_tick;
        this_end = floor(min(heasig.nsamp,plotXmax)/minor_tick)*minor_tick;
        minor_tick_x = this_start:minor_tick:this_end;
        
        minor_tick_x = setdiff(minor_tick_x, major_tick_x);

        % voltage grid
        [~, major_tick_idx] = sort( abs(( (plotYrange / gains(1)) ./major_tick_values_voltage) - 10));
        major_tick = major_tick_values_voltage(major_tick_idx(1));
        [voltage_major_tick, voltage_major_tick_preffix]= microVoltsTransformer(major_tick);
        
        major_tick_voltage = bottom_frame:major_tick * gains(1):top_frame;
        minor_tick = major_tick * gains(1) / 5;
        minor_tick_voltage = bottom_frame:minor_tick:top_frame;
        
        minor_tick_voltage = setdiff(minor_tick_voltage, major_tick_voltage);

        bAux = ishandle(paperModeHdl);
        if( any(bAux) )
            delete(paperModeHdl(bAux));
            paperModeHdl = {};
        end
        
        hold(axes_hdl, 'on');

        % background
        % left
        paperModeHdl = [paperModeHdl; ...
            {patch('Faces', [1 2 3 4], ... 
            'Vertices', [ [plotXmin; plotXmin; left_frame; left_frame ] [plotYmin; plotYmax; plotYmax; plotYmin ] ], ... 
            'FaceColor', [1 1 1], 'EdgeColor', [1 1 1] )} ];
        
        % right
        paperModeHdl = [paperModeHdl; ...
            {patch('Faces', [1 2 3 4], ... 
            'Vertices', [ [right_frame; right_frame; plotXmax; plotXmax ] [plotYmin; plotYmax; plotYmax; plotYmin ] ], ... 
            'FaceColor', [1 1 1], 'EdgeColor', [1 1 1] )} ];
        
        % top
        paperModeHdl = [paperModeHdl; ...
            {patch('Faces', [1 2 3 4], ... 
            'Vertices', [ [plotXmin; plotXmax; plotXmax; plotXmin ] [plotYmax; plotYmax; top_frame; top_frame ] ], ... 
            'FaceColor', [1 1 1], 'EdgeColor', [1 1 1] )} ];
        
        % bottom
        paperModeHdl = [paperModeHdl; ...
            {patch('Faces', [1 2 3 4], ... 
            'Vertices', [ [plotXmin; plotXmax; plotXmax; plotXmin ] [bottom_frame; bottom_frame; plotYmin; plotYmin ] ], ... 
            'FaceColor', [1 1 1], 'EdgeColor', [1 1 1] )} ];
        
        % grid
        % time
        aux_grid_idx = length(paperModeHdl) + 1;
        paperModeHdl = [paperModeHdl; colvec(arrayfun(@(a)(plot(axes_hdl, repmat(a,2,1), [bottom_frame; top_frame], '-', 'Color', [1 0.6 0.6], 'LineWidth', 0.5 )), major_tick_x, 'UniformOutput', false)) ];
        paperModeHdl = [paperModeHdl; colvec(arrayfun(@(a)(plot(axes_hdl, repmat(a,2,1), [bottom_frame; top_frame], ':', 'Color', [1 0.6 0.6], 'LineWidth', 0.3 )), minor_tick_x, 'UniformOutput', false))];
        
        %voltage
        paperModeHdl = [paperModeHdl; colvec(arrayfun(@(a)(plot(axes_hdl, [left_frame; right_frame], repmat(a,2,1), 'Color', [1 0.6 0.6], 'LineWidth', 0.5 )), major_tick_voltage, 'UniformOutput', false)) ];
        paperModeHdl = [paperModeHdl; colvec(arrayfun(@(a)(plot(axes_hdl, [left_frame; right_frame], repmat(a,2,1), ':', 'Color', [1 0.6 0.6], 'LineWidth', 0.3 )), minor_tick_voltage, 'UniformOutput', false)) ];
        
        aux_grid_idx = aux_grid_idx:length(paperModeHdl);

        % Arrow function caughts an error situation that is preferred to
        % avoid.
        db_status = dbstatus();
        bRestoreErrorStatus = false;
        if( length(db_status) > 1 && strcmpi(db_status(end).cond, 'caught error')  )
            dbclear if caught error
            bRestoreErrorStatus = true;
        end
        
        % voltage scale
        if( length(major_tick_voltage) > 2 )
            Vscale_y = mean(major_tick_voltage([2,3]));
            Vscale_x = right_frame - 2*xTextOffset;
            Vscale_arrow_x = right_frame - 1*xTextOffset;
            paperModeHdl = [paperModeHdl; {text( Vscale_x, Vscale_y , sprintf([ '%d ' voltage_major_tick_preffix 'V' ], voltage_major_tick), 'FontSize', 8, 'HorizontalAlignment', 'center', 'Rotation', 90, 'BackGroundColor', [1 1 1], 'EdgeColor', [1 1 1] )}];
            paperModeHdl = [paperModeHdl; {arrow( [Vscale_arrow_x; major_tick_voltage(2)], [Vscale_arrow_x; major_tick_voltage(3)], 2, 1, [0 0 0], axes_hdl )}];
        end
        
        % restore error status
        if(bRestoreErrorStatus)
            dbstop if caught error
        end
        
        % time reference
        if( (plotXrange/heasig.freq) > 20  )
            precision = 0;
        elseif( (plotXrange/heasig.freq) > 10  )
            precision = 1;
        elseif( (plotXrange/heasig.freq) > 5  )
            precision = 2;
        else
            precision = 3;
        end
        
        if( ~isempty(major_tick_x) )
            
            paperModeHdl = [ ... 
                                paperModeHdl; ...
                                text( major_tick_x(1), bottom_frame + yTextOffset, Seconds2HMS( (major_tick_x(1) + base_time )/heasig.freq, precision ) , 'FontSize', 8, 'HorizontalAlignment', 'Left', 'BackgroundColor', [1 1 1]) ...
                                ];

            paperModeHdl = [ ... 
                                paperModeHdl; ...
                                text( major_tick_x(end), bottom_frame + yTextOffset, Seconds2HMS( (major_tick_x(end) + base_time )/heasig.freq, precision ) , 'FontSize', 8, 'HorizontalAlignment', 'Right', 'BackgroundColor', [1 1 1]) ...
                                ];

            paperModeHdl = [ ... 
                                paperModeHdl; ...
                                colvec(arrayfun( @(a)(text( a, bottom_frame + yTextOffset, Seconds2HMS( (a - major_tick_x(1) )/heasig.freq, precision ) , 'FontSize', 8, 'HorizontalAlignment', 'center', 'BackgroundColor', [1 1 1])), major_tick_x(2:end-1), 'UniformOutput', false)) ...
                                ];
        end
        
        % frame
        aux_frame_idx = length(paperModeHdl) + 1;
        
        paperModeHdl = [paperModeHdl; {plot(axes_hdl, [left_frame left_frame right_frame right_frame left_frame ], [bottom_frame top_frame top_frame bottom_frame bottom_frame], 'Color', [1 0.6 0.6], 'LineWidth', 2.5 )}];
        
        aux_frame_idx = aux_frame_idx:length(paperModeHdl);
        
        hold(axes_hdl, 'off');
        
        cellfun(@(a)(uistack(a,'bottom')), paperModeHdl(aux_grid_idx));
        cellfun(@(a)(uistack(a,'top')), [paperModeHdl(aux_frame_idx); topLevelHdl]);
        
    end
%--------------------------------------------------------------------------

%==========================================================================
    function WindowKeyReleaseCallback2D(src, evnt)    %#ok
        %WindowKeyReleaseCallback2D
        
        switch evnt.Key
            case {'leftarrow', 'rightarrow', 'uparrow', 'downarrow'}
                mDragShiftStep = mDragSaveShiftStep;
            case 'o'
                scroll_mode = 'zoom';
%                 update_title_efimero( 'Zoom', 5 );
                
            case 'g'
                scroll_mode = 'zoom';
%                 update_title_efimero( 'Zoom', 5 );
                
            case 'm'
                MagnifierOff();
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function DragMouseBegin()
        %DragMouseBegin begin draging
        
        if (~fIsDragAllowed && ~fIsMagnifierOn)
            
            if( bPaperModeOn )
                PaperModeOff();
            end
            
            [cx, cy] = GetCursorCoordOnWindow();
            
            mDragStartX = cx;
            mDragStartY = cy;
            
            fIsDragAllowed = true;
            PrevStateWindowButtonMotionFcn = get(fig_hdl, 'WindowButtonMotionFcn');
            set(fig_hdl, 'WindowButtonMotionFcn', @WindowButtonMotionCallback2D);
            
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function DragMouseEnd()
        %DragMouseEnd end draging

        if fIsDragAllowed
            fIsDragAllowed = false;
            set(fig_hdl, 'WindowButtonMotionFcn', PrevStateWindowButtonMotionFcn);
            
            UserChangeView( fig_hdl, [], 'pan');
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function DragMouse()
        %DragMouse

        if fIsDragAllowed
            [cx, cy] = GetCursorCoordOnWindow();
            
            pdx = mDragStartX - cx;
            pdy = mDragStartY - cy;
            
            mDragStartX = cx;
            mDragStartY = cy;
            
            DragAxes(pdx, pdy);
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function DragKeys(direction)
        %DragKeys
        
        dx = mDragShiftStep;
        dy = mDragShiftStep;
        
        % Increment of speed when you hold the button
        mDragShiftStep = mDragShiftStep + mDragShiftStepInc;
        
        directionsX = {'right', 'left'};
        directionsY = {'down', 'up'};
        
        switch mDragKeysX
            case 'normal'
            case 'reverse'
                directionsX = fliplr(directionsX);
        end
        switch mDragKeysY
            case 'normal'
            case 'reverse'
                directionsY = fliplr(directionsY);
        end
        
        switch direction
            case directionsX{1}
                DragAxes(-dx, 0);
            case directionsX{2}
                DragAxes(dx, 0);
            case directionsY{1}
                DragAxes(0, dy);
            case directionsY{2}
                DragAxes(0, -dy);
        end
        
        PointerCrossUpdate();
        UserChangeView( fig_hdl, [], 'pan');
        
    end
%--------------------------------------------------------------------------

%==========================================================================
    function DragAxes(pdx, pdy)
        %DragAxes
        
        [xLim, yLim] = GetAxesLimits();
        
        pos = GetObjPos(axes_hdl, 'Pixels');
        pbar = get(axes_hdl, 'PlotBoxAspectRatio');
        
        %NOTE: MATLAB Bug?
        % Fixed problem with AspectRatio and Position of Axes
        % MATLAB Function PAN is not correct works with rectangular images!
        % Here it is correctly.
        
        imAspectRatioX = pbar(2) / pbar(1);
        if (imAspectRatioX ~= 1)
            posAspectRatioX = pos(3) / pos(4);
            arFactorX = imAspectRatioX * posAspectRatioX;
            if (arFactorX < 1)
                arFactorX = 1;
            end
        else
            arFactorX = 1;
        end
        
        imAspectRatioY = pbar(1) / pbar(2);
        if (imAspectRatioY ~= 1)
            posAspectRatioY = pos(4) / pos(3);
            arFactorY = imAspectRatioY * posAspectRatioY;
            if (arFactorY < 1)
                arFactorY = 1;
            end
        else
            arFactorY = 1;
        end
        
        if fIsEnableZoomX
            % For log plots, transform to linear scale
            if strcmp(get(axes_hdl, 'xscale'), 'log')
                xLim = log10(xLim);
                xLim = FixInfLogLimits('x', xLim);
                isXLog = true;
            else
                isXLog = false;
            end
            
            dx = pdx * range(xLim) / (pos(3) / arFactorX);
            xLim = xLim + dx;
            
            % For log plots, untransform limits
            if isXLog
                xLim = 10.^(xLim);
            end
        end
        if fIsEnableZoomY
            if strcmp(get(axes_hdl, 'yscale'), 'log')
                yLim = log10(yLim);
                yLim = FixInfLogLimits('y', yLim);
                isYLog = true;
            else
                isYLog = false;
            end
            
            dy = pdy * range(yLim) / (pos(4) / arFactorY);
            
            yLim = yLim + dy; 
            
            if isYLog
                yLim = 10.^(yLim);
            end
        end
        
        SetAxesLimits(xLim, yLim);
    end
%--------------------------------------------------------------------------


%==========================================================================
    function ZoomMouse(direction)
        %ZoomMouse zooming axes with mouse
        
        if (IsZoomMouseAllowed && ~fIsMagnifierOn)
            [acx, acy] = GetCursorCoordOnAxes();
            ZoomAxes(direction, acx, acy)
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function ZoomKeys(direction)
        %ZoomKeys zooming axes with keyboard
        
        if( bPaperModeOn )
            PaperModeOff();
        end
        
        UpdateCurrentZoomAxes();
        
        [mZoomGrid, mZoomSteps] = ZoomLogGrid(mZoomMinPow, mZoomMaxPow, mZoomKeysNum);
        UpdateCurrentZoomAxes();
        
        [acx, acy] = GetCursorCoordOnAxes();
        
        ZoomAxes(direction, acx, acy)
        PointerCrossUpdate();
        SetDefaultZoomGrid();
    end
%--------------------------------------------------------------------------

%==========================================================================
    function ZoomAxes(direction, cx, cy)
        %ZoomAxes Zoom axes in 2D and image modes
        
        [xLim, yLim] = GetAxesLimits();
        
        if (fIsEnableZoomX || (~fIsEnableZoomX && ~fIsEnableZoomY ) )
            mZoomIndexX = ChangeZoomIndex(direction, mZoomIndexX);
            zoomPct = GetZoomPercent(mZoomIndexX);

            xLim = RecalcZoomAxesLimits('x', xLim, mDefaultXLim, cx, zoomPct);
        end
        if (fIsEnableZoomY || (~fIsEnableZoomX && ~fIsEnableZoomY ) )
            mZoomIndexY = ChangeZoomIndex(direction, mZoomIndexY);
            zoomPct = GetZoomPercent(mZoomIndexY);

            yLim = RecalcZoomAxesLimits('y', yLim, mDefaultYLim, cy, zoomPct);
        end
        
        SetAxesLimits(xLim, yLim);
                
    end
%--------------------------------------------------------------------------

%==========================================================================
    function zoomPct = GetZoomPercent(zoomIndex, zoomGrid)
        %GetZoomPercent get zoom percent

        if (nargin < 2)
            zoomGrid = mZoomGrid;
        end
        
        zoomPct = zoomGrid(zoomIndex);       
    end
%--------------------------------------------------------------------------

%==========================================================================
    function zoomIndex = ChangeZoomIndex(direction, zoomIndex, zoomSteps)
        %ChangeZoomIndex
        
        if (nargin < 3)
            zoomSteps = mZoomSteps;
        end
        
        switch direction
            case 'plus'
                if (zoomIndex < zoomSteps)
                    zoomIndex = zoomIndex + 1;
                end
            case 'minus'
                if (zoomIndex > 1)
                    zoomIndex = zoomIndex - 1;
                end
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function axLim = RecalcZoomAxesLimits(ax, axLim, axLimDflt, zcCrd, zoomPct)
        %RecalcZoomAxesLimits recalc axes limits
        
        if strcmp(get(axes_hdl, [ax, 'scale']), 'log')
            axLim = log10(axLim);
            axLim = FixInfLogLimits(ax, axLim);
            axLimDflt = log10(axLimDflt);
            zcCrd = log10(zcCrd);
            isLog = true;
        else
            isLog = false;
        end
                
        if (zcCrd < axLim(1)), zcCrd = axLim(1); end
        if (zcCrd > axLim(2)), zcCrd = axLim(2); end
        
        rf = range(axLim);
        ra = range([axLim(1), zcCrd]);
        rb = range([zcCrd, axLim(2)]);
        
        cfa = ra / rf; 
        cfb = rb / rf;
        
        newRange = range(axLimDflt) * 100 / zoomPct;
        dRange = newRange - rf;
        
        axLim(1) = axLim(1) - dRange * cfa;
        axLim(2) = axLim(2) + dRange * cfb;
        
        if isLog
            axLim = 10.^axLim;
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function UpdateCurrentZoomAxes()
        %UpdateCurrentZoomAxes
        
        [xLim, yLim] = GetAxesLimits();
        [curentZoomX, curentZoomY] = GetCurrentZoomAxesPercent(xLim, yLim);
        
        if (curentZoomX ~= GetZoomPercent(mZoomIndexX))
            [nu, mZoomIndexX] = min(abs(mZoomGrid - curentZoomX));  %#ok ([~, ...])
        end
        if (curentZoomY ~= GetZoomPercent(mZoomIndexY))
            [nu, mZoomIndexY] = min(abs(mZoomGrid - curentZoomY));  %#ok ([~, ...])
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function [curentZoomX, curentZoomY] = GetCurrentZoomAxesPercent(xLim, yLim)
        %GetCurrentZoomAxesPercent
        
        if strcmp(get(axes_hdl, 'xscale'), 'log')
            xLim = log10(xLim);
            defaultXLim = log10(mDefaultXLim);
        else
            defaultXLim = mDefaultXLim;
        end
        if strcmp(get(axes_hdl, 'yscale'), 'log')
            yLim = log10(yLim);
            defaultYLim = log10(mDefaultYLim);
        else
            defaultYLim = mDefaultYLim;
        end
        
        curentZoomX = range(defaultXLim) * 100 / range(xLim);
        curentZoomY = range(defaultYLim) * 100 / range(yLim);
    end
%--------------------------------------------------------------------------


%==========================================================================
    function SetDefaultZoomGrid()
        %SetDefaultZoomGrid set default zoom grid
        
        [mDefaultZoomGrid, mDefaultZoomSteps] = ...
            ZoomLogGrid(mZoomMinPow, mZoomMaxPow, mZoomNum);
        
        mZoomGrid = mDefaultZoomGrid;
        mZoomSteps = mDefaultZoomSteps;
        
        mZoomIndexX = find(mZoomGrid == 100);
        mZoomIndexY = mZoomIndexX;
        mZoom3DIndex = mZoomIndexX;
    end
%--------------------------------------------------------------------------


%==========================================================================
    function PointerCrossOn()
        %PointerCrossOn
        
        if ~fIsPointerCross
            SetPointer('fullcrosshair');
            
            % text objects
            h = [];
            
            for jj = 1:cant_leads
                h = [ h text('Parent', axes_hdl, 'BackgroundColor', bgColor, 'Color', ColorOrder(jj,:), 'EdgeColor', ColorOrder(jj,:)) ];
            end
            h = [ h text('Parent', axes_hdl, 'BackgroundColor', bgColor, 'EdgeColor', [0 0 0] ) ];
            
            % create pointer cross struct
            mPointerCross = struct(...
                'htext',    h ...
                );
            
            PointerCrossSetup();
            fIsPointerCross = true;
            PointerCrossUpdate();
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function PointerCrossOff()
        %PointerCrossOff
        
        if fIsPointerCross
            delete(mPointerCross.htext);
            SetPointer('arrow');
            fIsPointerCross = false;
            set(fig_hdl, 'WindowButtonMotionFcn', []);
            mPointerCross = [];
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function PointerCrossSetup()
        %PointerCrossSetup
        set(fig_hdl, 'WindowButtonMotionFcn', @WindowButtonMotionCallback2D);
        
    end
%--------------------------------------------------------------------------

%==========================================================================
    function PointerCrossUpdate()
        %PointerCrossUpdate
        if fIsPointerCross
            [this_xlim, this_ylim] = GetAxesLimits();
            [acx, acy] = GetCursorCoordOnAxes();
            
            acx = min( cant_samp, max(1, round(acx) - start_sample + 1 ));
            
            % each lead
            extents = nan(cant_leads,4);
            
            for jj = 1:cant_leads
                
                if( any(jj == volt_idx) )
                    [aux_val, str_unit_prefix]= microVoltsTransformer(ECG( acx, jj));
                    set(mPointerCross.htext(jj), 'String', sprintf(['%3.0f ' str_unit_prefix 'V'], aux_val ) );
                else
                    set(mPointerCross.htext(jj), 'String', sprintf('%3.2f %s', ECG( acx, jj), heasig.units(jj,:) ) );
                end
                extents(jj,:) = get(mPointerCross.htext(jj), 'Extent');
            end
            % time
            if( (diff(this_xlim)/heasig.freq) > 20  )
                precision = 0;
            else
                precision = 3;
            end
            
            set(mPointerCross.htext(cant_leads+1), 'String', Seconds2HMS( (acx + start_sample - 1 + base_time )/heasig.freq, precision));
            
            % each lead
            for jj = 1:cant_leads
                set(mPointerCross.htext(jj), 'Position', [ this_xlim(2) - extents(jj,3), -offsets(jj) ] );
            end
            % time
            extents = get(mPointerCross.htext(cant_leads+1), 'Extent');
            set(mPointerCross.htext(cant_leads+1), 'Position', [acx + start_sample - 1 + 0.5*extents(3) this_ylim(1)] );
            
            uistack(mPointerCross.htext,'top');
            
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function RubberBandBegin()
        %RubberBandBegin
        
        if (~fIsRubberBandOn && ~fIsMagnifierOn)
            [acx, acy] = GetCursorCoordOnAxes();
            
            % create rubber band struct
            mRubberBand = struct(...
                'obj',	[patch('Parent', axes_hdl), patch('Parent', axes_hdl)], ...
                'txt_start_hdl', text('String', Seconds2HMS((acx+ base_time )/heasig.freq, 3), 'Parent', axes_hdl, 'BackgroundColor', bgColor, 'EdgeColor', [0 0 0] ), ...                
                'txt_duration_hdl', text('String', Seconds2HMS(0, 0), 'Parent', axes_hdl, 'BackgroundColor', bgColor, 'EdgeColor', [0 0 0] ), ...                
                'txt_amp_hdl', text('String', '0', 'Parent', axes_hdl, 'BackgroundColor', bgColor, 'EdgeColor', [0 0 0] ), ...                
                'txt_end_hdl', text('String', Seconds2HMS((acx+ base_time)/heasig.freq, 3), 'Parent', axes_hdl, 'BackgroundColor', bgColor, 'EdgeColor', [0 0 0] ), ...                
                'x1',  	acx, ...
                'y1',  	acy, ...
                'x2',  	acx, ...
                'y2',  	acy);

            extents = get(mRubberBand.txt_start_hdl, 'Extent');
            set(mRubberBand.txt_start_hdl, 'Position', [acx - extents(3) acy - extents(4)] );
            
            extents = get(mRubberBand.txt_end_hdl, 'Extent');
            set(mRubberBand.txt_end_hdl, 'Position', [acx acy - extents(4)] );
            
            extents = get(mRubberBand.txt_duration_hdl, 'Extent');
            set(mRubberBand.txt_duration_hdl, 'Position', [acx + 0.5 * extents(3) acy - extents(4)] );
            
            set(mRubberBand.txt_amp_hdl, 'Position', [acx  acy ] );
            
            hAxes2d = GetHandlesAxes2D();
            if ~isempty(hAxes2d)
                set(hAxes2d, ...
                    'XLimMode', 'manual', ...
                    'YLimMode', 'manual');
            end
            
            RubberBandSetPos();
            RubberBandSetup();
            fIsRubberBandOn = true;
            PrevStateWindowButtonMotionFcn = get(fig_hdl, 'WindowButtonMotionFcn');
            set(fig_hdl, 'WindowButtonMotionFcn', @WindowButtonMotionCallback2D);
            
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function RubberBandEnd()
        %RubberBandEnd
        
        if fIsRubberBandOn
            fIsRubberBandOn = false;
            set(fig_hdl, 'WindowButtonMotionFcn', PrevStateWindowButtonMotionFcn);
            
            delete(mRubberBand.obj);          
            delete(mRubberBand.txt_start_hdl);          
            delete(mRubberBand.txt_end_hdl);          
            delete(mRubberBand.txt_duration_hdl);  
            delete(mRubberBand.txt_amp_hdl);  
            
            if(~fNoZoom)
                RubberBandZoomAxes();
            end
            PointerCrossUpdate();
            mRubberBand = [];
            
            if(~fNoZoom)
                UserChangeView(fig_hdl, [], 'zoom');            
            end
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function RubberBandUpdate()
        %RubberBandUpdate
        
        if fIsRubberBandOn
            [acx, acy] = GetCursorCoordOnAxes();
            
            if( fIsEnableZoomX )
                mRubberBand.x2 = acx;
                mRubberBand.y2 = mRubberBand.y1;
            elseif( fIsEnableZoomY )
                mRubberBand.y2 = acy;
                mRubberBand.x2 = mRubberBand.x1;
            else
                mRubberBand.x2 = acx;
                mRubberBand.y2 = acy;
            end
            RubberBandSetPos();
            
            this_start = min(mRubberBand.x1, mRubberBand.x2);
            this_end = max(mRubberBand.x1, mRubberBand.x2);
            this_dur = (this_end - this_start);
            
            aux_dur = (this_dur/heasig.freq);
            if( aux_dur < 1  )
                time_precision = 3;
                this_dur_str = sprintf( '%3.0f ms', aux_dur*1e3);
            else
                if( aux_dur < 2  )
                    time_precision = 3;
                elseif( aux_dur < 5  )
                    time_precision = 2;
                elseif( aux_dur < 10  )
                    time_precision = 1;
                else
                    time_precision = 0;
                end
                this_dur_str = Seconds2HMS( aux_dur, time_precision);
            end
            
            upper_part = max( mRubberBand.y1, mRubberBand.y2);
            lower_part = min( mRubberBand.y1, mRubberBand.y2);
            this_amp = (upper_part - lower_part)/gains(1);
            if(this_amp > 999)
                %milli
                str_unit_prefix = 'm';
                this_amp = this_amp / 1e3;
            elseif(this_amp > 999999)
                str_unit_prefix = '';
                this_amp = this_amp / 1e6;
            else
                %micro volts per default
                str_unit_prefix = '\\mu';
            end

            % check decimal precision now
            if(this_amp < 9)
                amp_decs = '2';
            elseif(this_amp < 99)
                amp_decs = '1';
            else
                amp_decs = '0';
            end

            % time
            
            set(mRubberBand.txt_start_hdl, 'String', Seconds2HMS( (this_start+ base_time)/heasig.freq, time_precision) );
            extents = get(mRubberBand.txt_start_hdl, 'Extent');
            set(mRubberBand.txt_start_hdl, 'Position', [this_start - extents(3) lower_part - extents(4)] );
            
            set(mRubberBand.txt_end_hdl, 'String', Seconds2HMS( (this_end+ base_time)/heasig.freq, time_precision) );
            extents = get(mRubberBand.txt_end_hdl, 'Extent');
            set(mRubberBand.txt_end_hdl, 'Position', [this_end lower_part - extents(4)] );
            
            set(mRubberBand.txt_duration_hdl, 'String', this_dur_str );
            extents = get(mRubberBand.txt_duration_hdl, 'Extent');
            set(mRubberBand.txt_duration_hdl, 'Position', [ this_start + this_dur/2  upper_part + extents(4)] );
            
            % amplitude
            set(mRubberBand.txt_amp_hdl, 'String', sprintf(['%3.' amp_decs 'f ' str_unit_prefix 'V'], this_amp ) );
            extents = get(mRubberBand.txt_amp_hdl, 'Extent');
            set(mRubberBand.txt_amp_hdl, 'Position', [ this_end + 0.3*extents(3) lower_part + (upper_part - lower_part)/2 ] );
            
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function RubberBandSetPos()
        %RubberBandSetPos set position of rubber band
        
        x1 = mRubberBand.x1;
        y1 = mRubberBand.y1;
        x2 = mRubberBand.x2;
        y2 = mRubberBand.y2;
        
        set(mRubberBand.obj, ...
            'XData', [x1 x2 x2 x1], ...
            'YData', [y1 y1 y2 y2]);
    end
%--------------------------------------------------------------------------

%==========================================================================
    function RubberBandSetup()
        %RubberBandSetup
        
        set(mRubberBand.obj(1), ...
            'EdgeColor', 'w', ...
            'FaceColor', 'none', ...
            'LineWidth', 1.5, ...
            'LineStyle', '-');
        
        set(mRubberBand.obj(2), ...
            'EdgeColor', mRbEdgeColor, ...
            'FaceColor', mRbFaceColor, ...
            'FaceAlpha', mRbFaceAlpha, ...
            'LineWidth', 0.5, ...
            'LineStyle', '-');    
    end
%--------------------------------------------------------------------------

%==========================================================================
    function RubberBandZoomAxes()
        %RubberBandZoomAxes apply zoom from rubber band
        
        if( fIsEnableZoomY )
            xLim = get(axes_hdl, 'Xlim'); 
        else
            xLim = sort([mRubberBand.x1, mRubberBand.x2]);
        end
        
        if( fIsEnableZoomX )
            yLim = get(axes_hdl, 'Ylim'); 
        else
            yLim = sort([mRubberBand.y1, mRubberBand.y2]);
        end
        
        if (range(xLim) == 0 || range(yLim) == 0)
            return;
        end
        
        [zoomPctX, zoomPctY] = GetCurrentZoomAxesPercent(xLim, yLim);
        
        cx = mean(xLim);
        cy = mean(yLim);
        
        xLim = RecalcZoomAxesLimits('x', xLim, mDefaultXLim, cx, zoomPctX);
        yLim = RecalcZoomAxesLimits('y', yLim, mDefaultYLim, cy, zoomPctY);
        
        SetAxesLimits(xLim, yLim);
    end
%--------------------------------------------------------------------------

%==========================================================================
    function MagnifierOn()
        %MagnifierCreate
        
        if ~fIsMagnifierOn    
            
            if( bPaperModeOn )
                PaperModeOff();
            end
            
            if fIsPointerCross
                isPointerCross = true;
                PointerCrossOff();
            else
                isPointerCross = false;
            end
            
            mMgDirection = 'plus';
            
            % create magnifier struct
            mMagnifier = struct(...
                'obj',          copyobj(axes_hdl, fig_hdl), ...
                'frame_obj',    [], ...
                'size',         mMgSize, ...
                'zoom',         mMgZoom);
            
            fIsMagnifierOn = true;
            
            set(fig_hdl, 'WindowButtonMotionFcn', @WindowButtonMotionCallback2D);
            
            MagnifierSetup();
            MagnifierUpdate();
            
            if isPointerCross
                PointerCrossOn();
            end
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function MagnifierOff()
        %MagnifierOff
        
        if fIsMagnifierOn
            fIsMagnifierOn = false;
            set(fig_hdl, 'WindowButtonMotionFcn', PrevStateWindowButtonMotionFcn);
            
            set(axes_hdl, 'Color', get(mMagnifier.obj, 'Color'));
            
            delete(mMagnifier.obj);
            mMagnifier = [];
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function MagnifierUpdate()
        %MagnifierUpdate
        
        if fIsMagnifierOn            
            % see original idea of magnify by Rick Hindman -- 7/29/04
            % http://www.mathworks.com/matlabcentral/fileexchange/5961
            
            [acx, acy] = GetCursorCoordOnAxes();
            [wcx, wcy] = GetCursorCoordOnWindow('pixels');
            
            [xLim, yLim] = GetAxesLimits();
            
            if strcmp(get(axes_hdl, 'xscale'), 'log')
                xLim = log10(xLim);
                xLim = FixInfLogLimits('x', xLim);
                acx = log10(acx);
                isXLog = true;
            else
                isXLog = false;
            end
            if strcmp(get(axes_hdl, 'yscale'), 'log')
                yLim = log10(yLim);
                yLim = FixInfLogLimits('y', yLim);
                acy = log10(acy);
                isYLog = true;
            else
                isYLog = false;
            end
            
            figPos = GetObjPos(fig_hdl, 'pixels');
            axPos = GetObjPos(axes_hdl, 'normalized');
            
            % always square magnifier
            pbar = get(axes_hdl, 'PlotBoxAspectRatio');
            af = pbar(1) / pbar(2);
            if (af == 1 && (pbar(1) == 1 && pbar(2) == 1))
                af = figPos(3) / figPos(4);
            end
            
            mgSizePix = round(mMagnifier.size);
            mgZoom = mMagnifier.zoom;
            
            mgSize = mgSizePix / figPos(3); % normalized size
            
            mgPos(3) = mgSize * 2;
            mgPos(4) = mgPos(3) * af;
            
            mg3 = round(mgPos(3) * figPos(3));
            mg4 = round(mgPos(4) * figPos(4));
            
            if (mg4 < mg3)
                mgSize = (mgSizePix * (mg3 / mg4)) / figPos(3);
            end
            
            mgPos(3) = mgSize * 2;
            mgPos(4) = mgPos(3) * af;
            
            mgPos(1) = wcx / figPos(3) - mgSize;
            mgPos(2) = wcy / figPos(4) - mgSize * af;
            
            mgXLim = acx + (1 / mgZoom) * (mgPos(3) / axPos(3)) * diff(xLim) * [-0.5 0.5];
            mgYLim = acy + (1 / mgZoom) * (mgPos(4) / axPos(4)) * diff(yLim) * [-0.5 0.5];
            
            SetObjPos(mMagnifier.obj, mgPos, 'normalized');
            
            if isXLog
                mgXLim = 10.^mgXLim;
            end
            if isYLog
                mgYLim = 10.^mgYLim;
            end
            
            set(mMagnifier.obj, ...
                'XLim', mgXLim, ...
                'YLim', mgYLim);
            
            MagnifierBorderUpdate();
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function MagnifierSetup()
        %MagnifierSetup

        set(mMagnifier.obj, ...
            'Box', 'on', ...
            'XMinorTick', 'on', ...
            'YMinorTick', 'on');
        
        title(mMagnifier.obj, '');
        xlabel(mMagnifier.obj, ''); 
        ylabel(mMagnifier.obj, '');
        
        
        hLines = findobj(mMagnifier.obj, 'Type', 'line');
        if ~isempty(hLines)
            if (mMgLinesWidth ~= 1)
                set(hLines, 'LineWidth', mMgLinesWidth);
            end
        end
        
        set(axes_hdl, 'Color', get(axes_hdl, 'Color')*mMgShadow);
    end
%--------------------------------------------------------------------------

%==========================================================================
    function MagnifierBorderUpdate()
        %MagnifierBorderUpdate
        
    end
%--------------------------------------------------------------------------

%==========================================================================
    function MagnifierSizeChange(direction)
        %MagnifierSizeChange
        
        if fIsMagnifierOn
            switch direction
                case 'plus'
                    if (mMagnifier.size < mMgMaxSize)
                        mMagnifier.size = mMagnifier.size + mMgSizeStep;
                    end
                case 'minus'
                    if (mMagnifier.size > mMgMinSize)
                        mMagnifier.size = mMagnifier.size - mMgSizeStep;
                    end
            end
            
            MagnifierUpdate();
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function MagnifierZoomChange(direction)
        %MagnifierZoomChange
        
        if fIsMagnifierOn
            switch direction
                case 'plus'
                    if (mMagnifier.zoom < mMgMaxZoom)
                        mMagnifier.zoom = mMagnifier.zoom * mMgZoomStep;
                    end
                case 'minus'
                    if (mMagnifier.zoom > mMgMinZoom)
                        mMagnifier.zoom = mMagnifier.zoom / mMgZoomStep;
                    end
            end
            
            MagnifierUpdate();
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function MagnifierReset()
        %MagnifierReset
        
        if fIsMagnifierOn
            mMagnifier.size = mMgSize;
            mMagnifier.zoom = mMgZoom;
            MagnifierUpdate();
        end
    end
%--------------------------------------------------------------------------


%==========================================================================
    function ResetAxesToOrigView()
        %ResetAxesToOrigView reset axes to original limits
        
        if( bPaperModeOn )
            PaperModeOff();
        end
        
        SetAxesLimits(mDefaultXLim, mDefaultYLim);
        PointerCrossUpdate();
        
        mZoomIndexX = find(mZoomGrid == 100);
        mZoomIndexY = mZoomIndexX;
        UserChangeView( fig_hdl, [], 'zoom')        
    end
%--------------------------------------------------------------------------


%==========================================================================
    function [x, y, z] = GetCursorCoordOnAxes()
        %GetCursorCoordOnAxImg
        
        crd = get(axes_hdl, 'CurrentPoint');
        x = crd(2,1);
        y = crd(2,2);
        z = crd(2,3);
    end
%--------------------------------------------------------------------------

%==========================================================================
    function [x, y] = GetCursorCoordOnWindow(units)
        %GetCursorCoordOnWindow
        
        if (nargin < 1), units = 'pixels'; end
        
        dfltUnits = get(fig_hdl, 'Units');
        set(fig_hdl, 'Units', units);
        
        crd = get(fig_hdl, 'CurrentPoint');
        x = crd(1); 
        y = crd(2);
        
        set(fig_hdl, 'Units', dfltUnits);
    end
%--------------------------------------------------------------------------

%==========================================================================
    function pos = GetObjPos(h, units)
        %GetObjPos get object position
        
        if (nargin < 2), units = get(h, 'Units'); end
        
        dfltUnits = get(h, 'Units');
        set(h, 'Units', units);
        pos = get(h, 'Position');
        set(h, 'Units', dfltUnits);
    end
%--------------------------------------------------------------------------

%==========================================================================
    function SetObjPos(h, pos, units)
        %SetObjPos set object position
        
        if (nargin < 3), units = get(h, 'Units'); end
        
        dfltUnits = get(h, 'Units');
        set(h, 'Units', units);
        set(h, 'Position', pos);
        set(h, 'Units', dfltUnits);
    end
%--------------------------------------------------------------------------

%==========================================================================
    function [xLim, yLim] = GetAxesLimits()
        %GetAxesLimits
        
        xLim = get(axes_hdl, 'XLim');
        yLim = get(axes_hdl, 'YLim');
    end
%--------------------------------------------------------------------------

%==========================================================================
    function SetAxesLimits(xLim, yLim)
        %SetAxesLimits
        
        set(axes_hdl, 'XLim', xLim);
        set(axes_hdl, 'YLim', yLim);
    end
%--------------------------------------------------------------------------

%==========================================================================
    function SetPointerCrossKeys()
        %SetPointerCrossKeys set pointer fullcross
        
        if( bPaperModeOn )
            PaperModeOff();
        end
        
        if fIsPointerCross
            PointerCrossOff();
        else
            PointerCrossOn();
        end
        
        UserData = get(fig_hdl, 'UserData');
        UserData.tools.pointercross = mPointerCross;
        set(fig_hdl, 'UserData', UserData);
    end
%--------------------------------------------------------------------------

%==========================================================================
    function SetPointer(pointerType)
        %SetPointer set pointer symbol
        
        set(fig_hdl, 'Pointer', pointerType);
    end
%--------------------------------------------------------------------------

%==========================================================================
    function SetAxesGridKeys()
        %SetAxesGridKeys on/off axes grid
        
        if fIsAxesGrid
            action = 'off';
            fIsAxesGrid = false;
        else
            action = 'on';
            fIsAxesGrid = true;
        end
        
        set(axes_hdl, 'XGrid', action, 'YGrid', action, 'ZGrid', action);
        
        if fIsMagnifierOn
            set(mMagnifier.obj, 'XGrid', action, 'YGrid', action);
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function [zg, st] = ZoomLogGrid(a, b, n)
        %ZoomLogGrid log zoom grid
        
        zg = unique(round(logspace(a, b, n)));
        
        zg(zg<100) = [];	% begin zoom == 100%
        st = length(zg);
        
        if isempty(find(zg == 100, 1))
            error('dragzoom:badZoomGridOptions', 'Options for zoom grid is bad.')
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function tf = IsZoomMouseAllowed()
        %IsZoomMouseAllowed
        
        [wcx, wcy] = GetCursorCoordOnWindow();
        figPos = get(fig_hdl, 'Position');
        
        if (wcx >= 1 && wcx <= figPos(3) && wcy >= 1 && wcy <= figPos(4))
            tf = true;
        else
            tf = false;
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function SelectAxesUnderCursor()
        %SelectAxesUnderCursor select axes under cursor as current
        
        axi = GetAxesIndexUnderCursor();
        
        if (axi > 0)
            fIsEnableControl = true;
            
            if ~mAxesInfo(axi).iscurrent
                caxi = GetCurrentAxesIndex();
                
                if isempty(caxi)
                    DeleteInvalidAxesInfo();
                    
                    axi = GetAxesIndexUnderCursor();
                    isCax2d = mAxesInfo(axi).is2d;
                else
                    isCax2d = mAxesInfo(caxi).is2d;
                end
                
                SetCurrentAxes(axi);
                
                % for fix "legend" axes capture
                if mAxesInfo(axi).islegend;
                    fIsMouseOnLegend = true;
                else
                    fIsMouseOnLegend = false;
                end
                
                % check callbacks
                if (isCax2d ~= mAxesInfo(axi).is2d)
                    % if dimension of axes has changed
                    SetCallbacks();
                    
                    if fIsPointerCross
                        % disable pointer cross
                        PointerCrossOff()
                    end
                else
                    if fIsPointerCross
                        % reset pointer cross
                        PointerCrossOff()
                        SetPointerCrossKeys()
                    end
                end
            end
        else
            fIsEnableControl = false;
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function SetCurrentAxes(axi)
        %SetCurrentAxes set current axes and work mode
        
        axes_hdl = mAxesInfo(axi).handle;
        
        set(fig_hdl, 'CurrentAxes', axes_hdl);
        for i = 1:numel(mAxesInfo)
            mAxesInfo(i).iscurrent = false;
        end
        mAxesInfo(axi).iscurrent = true;
        
        
        mDefaultAxPos = mAxesInfo(axi).position;
        mDefaultXLim = mAxesInfo(axi).xlim;
        mDefaultYLim = mAxesInfo(axi).ylim;
        
        % save info to work correctly after saving figures
        UserData = get(fig_hdl, 'UserData');
        UserData.axesinfo = mAxesInfo;
        set(fig_hdl, 'UserData', UserData);
    end
%--------------------------------------------------------------------------

%==========================================================================
    function axi = GetCurrentAxesIndex()
        %GetCurrentAxesIndex
        
        axi = [];
        
        for i = 1:numel(mAxesInfo)
            if (ishandle(mAxesInfo(i).handle) && mAxesInfo(i).iscurrent)
                axi = i;
                return;
            end
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function axi = GetAxesIndexUnderCursor()
        %FindAxesUnderCursor find current axes under cursor
        
        axi = GetCurrentAxesIndex();
        
        if ~fIsSelectedCurrentAxes
            caxi = GetCurrentAxesIndex();
            if ~IsInBoundsAxes(mAxesInfo(caxi).handle)
                axi = 0;
            end
            return;
        end
        
        for i = 1:numel(mAxesInfo)
            if (ishandle(mAxesInfo(i).handle) && IsInBoundsAxes(mAxesInfo(i).handle))
                axi = i;
                return;
            else
                axi = 0; % without axes
            end
        end
    end
%--------------------------------------------------------------------------


%==========================================================================
    function hAxes2d = GetHandlesAxes2D()
        %GetHandlesAxes2D Get handles of 2-D axes
        
        isAxes2d = arrayfun(@(x) x.is2d && ~x.islegend, mAxesInfo);
        hAxes2d = hAxes(isAxes2d);
        
        if ~isempty(hAxes2d)
            % Set current axes on first position
            hAxes2d(eq(hAxes2d, axes_hdl)) = 0;
            hAxes2d = sort(hAxes2d);
            hAxes2d(eq(hAxes2d, 0)) = axes_hdl;
        end
    end
%--------------------------------------------------------------------------


%==========================================================================
    function AxesInfo = GetAxesInfo()
        %GetAxesInfo make and get axes info struct
        
        countAxes = length(hAxes);
        
        AxesInfo = struct(...
            'handle',       cell(1, countAxes), ...
            'iscurrent',    cell(1, countAxes), ...
            'is2d',         cell(1, countAxes), ...
            'isimage',      cell(1, countAxes), ...  
            'isvisible',    cell(1, countAxes), ...  
            'isvis3d',      cell(1, countAxes), ...
            'islegend',     cell(1, countAxes), ...
            'position',     cell(1, countAxes), ...
            'normposition', cell(1, countAxes), ...
            'xlim',         cell(1, countAxes), ...
            'ylim',         cell(1, countAxes), ...
            'camtarget',    cell(1, countAxes), ...
            'camposition',  cell(1, countAxes));
        
        for i = 1:countAxes
            h = hAxes(i);
            
            AxesInfo(i).handle = h;
            AxesInfo(i).iscurrent = IsCurrentAxes(h);
            AxesInfo(i).is2d = IsAxes2D(h);
            AxesInfo(i).isimage = IsImageOnAxes(h);
            AxesInfo(i).isvisible = strcmpi(get(h, 'Visible'), 'on');
            AxesInfo(i).isvis3d = IsAxesVis3D(h);
            AxesInfo(i).islegend = IsLegendAxes(h);
            AxesInfo(i).position = GetObjPos(h, 'pixels');
            AxesInfo(i).normposition = GetObjPos(h, 'normalized');
            AxesInfo(i).xlim = get(h, 'XLim');
            AxesInfo(i).ylim = get(h, 'YLim');
            AxesInfo(i).camtarget = get(h, 'CameraTarget');
            AxesInfo(i).camposition = get(h, 'CameraPosition');
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function tf = IsImageOnAxes(ax)
        %IsImageOnAxes
        
        if (nargin < 1), ax = axes_hdl; end
        
        h = findobj(ax, 'Type', 'Image');
        
        if isempty(h)
            tf = false;
        else
            tf = true;
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function tf = IsAxes2D(ax)
        %IsAxes2D
        
        if (nargin < 1), ax = axes_hdl; end
        
        tf = is2D(ax); % (!!!) internal undocumented function
    end
%--------------------------------------------------------------------------

%==========================================================================
    function tf = IsLegendAxes(ax)
        %IsLegendAxes
        
        tf = strcmp(get(ax, 'Tag'), 'legend');
    end
%--------------------------------------------------------------------------

%==========================================================================
    function targetInBounds = IsInBoundsAxes(ax)
        %InBoundsAxes Check if the user clicked within the bounds of the axes. If not, do nothing
        
        targetInBounds = true;
        tol = 3e-16;
        cp = get(ax, 'CurrentPoint');
        
        XLims = get(ax, 'XLim');
        if ((cp(1,1) - min(XLims)) < -tol || (cp(1,1) - max(XLims)) > tol) && ...
                ((cp(2,1) - min(XLims)) < -tol || (cp(2,1) - max(XLims)) > tol)
            targetInBounds = false;
        end
        
        YLims = get(ax, 'YLim');
        if ((cp(1,2) - min(YLims)) < -tol || (cp(1,2) - max(YLims)) > tol) && ...
                ((cp(2,2) - min(YLims)) < -tol || (cp(2,2) - max(YLims)) > tol)
            targetInBounds = false;
        end
        
        ZLims = get(ax, 'ZLim');
        if ((cp(1,3) - min(ZLims)) < -tol || (cp(1,3) - max(ZLims)) > tol) && ...
                ((cp(2,3) - min(ZLims)) < -tol || (cp(2,3) - max(ZLims)) > tol)
            targetInBounds = false;
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function tf = IsCurrentAxes(ax)
        %IsCurrentAxes
        
        hcAx = get(fig_hdl, 'CurrentAxes');
        tf = eq(ax, hcAx);
    end
%--------------------------------------------------------------------------

%==========================================================================
    function tf = IsAxesVis3D(ax)
        %IsAxesVis3D
        
        visProp = {
            get(ax, 'PlotBoxAspectRatioMode')
            get(ax, 'DataAspectRatioMode')
            get(ax, 'CameraViewAngleMode')
            };
        
        tf = all(strcmpi(visProp, 'manual'));
    end
%--------------------------------------------------------------------------

%==========================================================================
    function maximize(fig)
        units=get(fig,'units');
        set(fig,'units','normalized','outerposition',[0 0.03 1 0.97]);
        set(fig,'units',units);
    end
%--------------------------------------------------------------------------

%==========================================================================
    function b=protected_index(a,idx)
        if(isempty(a))
            b = nan;
        else
            b = a(idx);
        end
    end
%--------------------------------------------------------------------------

%==========================================================================

    function this_hdl = PlotWaveMarks( this_annotation, field_names, lead, vertTextOffset, this_color)

        this_hdl = {};
       
        if( isfield(this_annotation, field_names{1} ) )
            aux_on = colvec(this_annotation.(field_names{1}));
            aux_on( aux_on < 1 | aux_on > heasig.nsamp) = nan;
        else
            aux_on = [];
        end
        
        if( isfield(this_annotation, field_names{2} ) )
            aux_peak = colvec(this_annotation.(field_names{2}));
            aux_peak( aux_peak < 1 | aux_peak > heasig.nsamp) = nan;
        else
            aux_peak = [];
        end
        
        if( isfield(this_annotation, field_names{3} ) )
            aux_off = colvec(this_annotation.(field_names{3}));
            aux_off( aux_off < 1 | aux_off > heasig.nsamp) = nan;
        else
            aux_off = [];
        end
        
%         aux_mat = match_positions(aux_mat);
        
        if( ann_graph_mode == kLinesAnns ) 
            
            % wave start
            bOn = ~isnan(aux_on) & aux_on >= start_sample & aux_on <= end_sample;
            aux_on_idx = find(bOn);
            this_hdl = [ this_hdl; colvec(arrayfun(@(a,b)(plot(axes_hdl, a, b, 'Color' , this_color, 'LineStyle', ':', 'Marker', '<' , 'MarkerSize', 2, 'LineWidth', 0.25) ), repmat(rowvec(aux_on(aux_on_idx)), 2, 1 ), bsxfun( @plus, repmat([-vertTextOffset; vertTextOffset], 1, length(aux_on_idx)), rowvec( (ECG(aux_on(aux_on_idx) - start_sample + 1, lead) * gains(lead)) - offsets(lead) ) ), 'UniformOutput', false ) ),  ];

            % wave end
            bOff = ~isnan(aux_off) & aux_off >= start_sample & aux_off <= end_sample;
            aux_off_idx = find(bOff);
            this_hdl = [ this_hdl; colvec(arrayfun(@(a,b)(plot(axes_hdl, a, b, 'Color' , this_color, 'LineStyle', ':', 'Marker', '>', 'MarkerSize', 2, 'LineWidth', 0.25 )), repmat(rowvec(aux_off(aux_off_idx)), 2, 1 ), bsxfun( @plus, repmat([-vertTextOffset; vertTextOffset], 1, length(aux_off_idx)), rowvec((ECG(aux_off(aux_off_idx) - start_sample + 1, lead) * gains(lead)) - offsets(lead)  ) ), 'UniformOutput', false)) ];

            % wave peak
            bPeak = ~isnan(aux_peak) & aux_peak >= start_sample & aux_peak <= end_sample;
            aux_peak_idx = find(bPeak);
            this_hdl = [ this_hdl; colvec(arrayfun(@(a,b)(plot(axes_hdl, a, b, 'Color' , this_color, 'LineStyle', ':', 'Marker', '^', 'MarkerSize', 2, 'LineWidth', 0.25 )), repmat(rowvec(aux_peak(aux_peak_idx)), 2, 1 ), bsxfun( @plus, repmat([-vertTextOffset; vertTextOffset]*1.3, 1, length(aux_peak_idx)), rowvec( (ECG(aux_peak(aux_peak_idx) - start_sample + 1, lead) * gains(lead)) - offsets(lead) ) ), 'UniformOutput', false)) ];
            this_hdl = [ this_hdl; arrayfun( @(a)(text(a + 0.5*xTextOffset, (ECG(a - start_sample + 1, lead) * gains(lead)) - offsets(lead) + sign(ECG(a - start_sample + 1, lead)) * yTextOffset, field_names{2}, 'FontSize', 8, 'Color', this_color ) ), aux_peak(aux_peak_idx), 'UniformOutput', false ) ];

            % wave conection between start-end
            aux_complete_idx = find(bOn & bOff);
            this_hdl = [ this_hdl; plot(axes_hdl,  repmat([rowvec(aux_on(aux_complete_idx));rowvec(aux_off(aux_complete_idx))], 1, 2 ) , ...
                                    [ [-vertTextOffset + rowvec( (ECG(aux_on(aux_complete_idx) - start_sample + 1, lead) * gains(lead)) - offsets(lead) ) ; -vertTextOffset + rowvec((ECG(aux_off(aux_complete_idx) - start_sample + 1, lead) * gains(lead)) - offsets(lead))] [ vertTextOffset + rowvec( (ECG(aux_on(aux_complete_idx) - start_sample + 1, lead) * gains(lead)) - offsets(lead) ); vertTextOffset + rowvec((ECG(aux_off(aux_complete_idx) - start_sample + 1, lead) * gains(lead)) - offsets(lead)) ] ], 'Color' , this_color, 'LineStyle', ':', 'Marker', 'none', 'LineWidth', 0.25 )];
                            
        elseif( ann_graph_mode == kBackColourAnns )
            
            % wave start
            bOn = ~isnan(aux_on) & aux_on >= start_sample & aux_on <= end_sample;
            % wave end
            bOff = ~isnan(aux_off) & aux_off >= start_sample & aux_off <= end_sample;
            
            % wave peak
            bPeak = ~isnan(aux_peak) & aux_peak >= start_sample & aux_peak <= end_sample;
            aux_peak_idx = find(bPeak);
            this_hdl = [ this_hdl; colvec(arrayfun(@(a,b)(plot(axes_hdl, a, b, 'Color' , this_color, 'LineStyle', ':', 'Marker', '^', 'MarkerSize', 2, 'LineWidth', 0.25 )), repmat(rowvec(aux_peak(aux_peak_idx)), 2, 1 ), bsxfun( @plus, repmat([-vertTextOffset; vertTextOffset]*1.3, 1, length(aux_peak_idx)), rowvec( (ECG(aux_peak(aux_peak_idx) - start_sample + 1, lead) * gains(lead)) - offsets(lead) ) ), 'UniformOutput', false)) ];
%             this_hdl = [ this_hdl; arrayfun( @(a)(text(a + 0.5*xTextOffset, (ECG(a - start_sample + 1, lead) * gains(lead)) - offsets(lead) + sign(ECG(a - start_sample + 1, lead)) * yTextOffset, field_names{2}, 'FontSize', 8, 'Color', this_color ) ), aux_peak(aux_peak_idx) ) ];

            % wave conection between start-end
            bOnOff = bOn & bOff;
            aux_complete_idx = find(bOnOff);
            aux_complete_idx2 = aux_complete_idx;
            aux_complete_idxx = arrayfun( @(a)( max(1, aux_on(a)):min(heasig.nsamp,aux_off(a)) ),aux_complete_idx, 'UniformOutput', false);
            
            % on-peak
            aux_complete_idx = find( ~bOnOff & bOn & bPeak);
            aux_complete_idx2 = [aux_complete_idx2;colvec(aux_complete_idx)];
            aux_complete_idxx = [ aux_complete_idxx; arrayfun( @(a)( aux_on(a):aux_peak(a) ),aux_complete_idx, 'UniformOutput', false) ];
            
            % peak-off
            aux_complete_idx = find(~bOnOff & bPeak & bOff);
            aux_complete_idx2 = [aux_complete_idx2;colvec(aux_complete_idx)];
            aux_complete_idxx = [ aux_complete_idxx; arrayfun( @(a)( aux_peak(a):min(heasig.nsamp,aux_off(a)) ),aux_complete_idx, 'UniformOutput', false) ];

            if( ~isempty(aux_complete_idxx) ) 
                
                aux_offset = (start_sample - 1);
                %patch around the signal
    %             this_hdl = [ this_hdl; cellfun( @(a)( patch( [a fliplr(a) ], [ (ECG(a-aux_offset, lead)* gains(lead) )- offsets(lead) + 0.5*yTextOffset ; flipud((ECG(a-aux_offset, lead)* gains(lead) )- offsets(lead)) - 0.5*yTextOffset ]', this_color, 'EdgeColor', 'none')), aux_complete_idxx) ];
                %box around the wave

                max_vals = cellfun( @(a)( max(ECG(a-aux_offset, lead)) ), aux_complete_idxx, 'UniformOutput', false);
                min_vals = cellfun( @(a)( min(ECG(a-aux_offset, lead)) ), aux_complete_idxx, 'UniformOutput', false);

                this_edge_color = repmat({0.8*this_color}, length(aux_complete_idxx), 1 );

                this_hdl = [ this_hdl; cellfun( @(a,b,c,d)( patch( [a(1) a(1) a(end) a(end) ], ( [ c b b c ] * gains(lead) )- offsets(lead), this_color, 'EdgeColor', d)), aux_complete_idxx, max_vals, min_vals, this_edge_color, 'UniformOutput', false) ];
            
            end
            
        end
        
%         uistack(this_hdl, 'bottom');
        
    end
%--------------------------------------------------------------------------


 %==========================================================================
   function this_hdl = PlotGlobalWaveMarks( field_names, limits, this_color)

        this_hdl = [];
        if( isfield(global_annotations, field_names{1} ) )
            aux_on = global_annotations.(field_names{1});
        else
            aux_on = [];
        end
        
        if( isfield(global_annotations, field_names{2} ) )
            aux_peak = global_annotations.(field_names{2});
        else
            aux_peak = [];
        end
        
        if( isfield(global_annotations, field_names{3} ) )
            aux_off = global_annotations.(field_names{3});
        else
            aux_off = [];
        end
                            
        if( ann_graph_mode == kLinesAnns ) 
            
            % wave start
            bOn = ~isnan(aux_on) & aux_on >= start_sample & aux_on <= end_sample;
            aux_on_idx = find(bOn);
%             this_hdl = [ this_hdl; plot(axes_hdl, repmat(rowvec(aux_on(aux_on_idx)), 2, 1 ), bsxfun( @plus, repmat([-vertTextOffset; vertTextOffset], 1, length(aux_on_idx)), rowvec( (ECG(aux_on(aux_on_idx) - start_sample + 1, lead) * gains(lead)) - offsets(lead) ) ), 'Color' , this_color, 'LineStyle', ':', 'Marker', '<' , 'MarkerSize', 2, 'LineWidth', 0.25)];
            this_hdl = [ this_hdl; plot(axes_hdl, repmat(rowvec(aux_on(aux_on_idx)), 2, 1 ), repmat(colvec(limits), 1, length(aux_on_idx)), 'Color' , this_color, 'LineStyle', ':', 'Marker', '<' , 'MarkerSize', 4, 'LineWidth', 0.25)];

            % wave end
            bOff = ~isnan(aux_off) & aux_off >= start_sample & aux_off <= end_sample;
            aux_off_idx = find(bOff);
%             this_hdl = [ this_hdl; plot(axes_hdl, repmat(rowvec(aux_off(aux_off_idx)), 2, 1 ), bsxfun( @plus, repmat([-vertTextOffset; vertTextOffset], 1, length(aux_off_idx)), rowvec((ECG(aux_off(aux_off_idx) - start_sample + 1, lead) * gains(lead)) - offsets(lead)  ) ), 'Color' , this_color, 'LineStyle', ':', 'Marker', '>', 'MarkerSize', 2, 'LineWidth', 0.25 )];
            this_hdl = [ this_hdl; plot(axes_hdl, repmat(rowvec(aux_off(aux_off_idx)), 2, 1 ), repmat(colvec(limits), 1, length(aux_off_idx)), 'Color' , this_color, 'LineStyle', ':', 'Marker', '>', 'MarkerSize', 4, 'LineWidth', 0.25 )];

            % wave peak
            bPeak = ~isnan(aux_peak) & aux_peak >= start_sample & aux_peak <= end_sample;
            aux_peak_idx = find(bPeak);
%             this_hdl = [ this_hdl; plot(axes_hdl, repmat(rowvec(aux_peak(aux_peak_idx)), 2, 1 ), bsxfun( @plus, repmat([-vertTextOffset; vertTextOffset]*1.3, 1, length(aux_peak_idx)), rowvec( (ECG(aux_peak(aux_peak_idx) - start_sample + 1, lead) * gains(lead)) - offsets(lead) ) ), 'Color' , this_color, 'LineStyle', ':', 'Marker', '^', 'MarkerSize', 2, 'LineWidth', 0.25 )];
            this_hdl = [ this_hdl; plot(axes_hdl, repmat(rowvec(aux_peak(aux_peak_idx)), 2, 1 ), repmat(colvec(limits) + [-0.02; 0.02] * diff(limits), 1, length(aux_peak_idx)), 'Color' , this_color, 'LineStyle', ':', 'Marker', '^', 'MarkerSize', 4, 'LineWidth', 0.25 )];
%             this_hdl = [ this_hdl; arrayfun( @(a)(text(a + 0.5*xTextOffset, (ECG(a - start_sample + 1, lead) * gains(lead)) - offsets(lead) + sign(ECG(a - start_sample + 1, lead)) * yTextOffset, field_names{2}, 'FontSize', 8, 'Color', this_color ) ), aux_peak(aux_peak_idx) ) ];

            % wave conection between start-end
            aux_complete_idx = find(bOn & bOff);
%             this_hdl = [ this_hdl; plot(axes_hdl,  repmat([rowvec(aux_on(aux_complete_idx));rowvec(aux_off(aux_complete_idx))], 1, 2 ) , ...
%                                     [ [-vertTextOffset + rowvec( (ECG(aux_on(aux_complete_idx) - start_sample + 1, lead) * gains(lead)) - offsets(lead) ) ; -vertTextOffset + rowvec((ECG(aux_off(aux_complete_idx) - start_sample + 1, lead) * gains(lead)) - offsets(lead))] [ vertTextOffset + rowvec( (ECG(aux_on(aux_complete_idx) - start_sample + 1, lead) * gains(lead)) - offsets(lead) ); vertTextOffset + rowvec((ECG(aux_off(aux_complete_idx) - start_sample + 1, lead) * gains(lead)) - offsets(lead)) ] ], 'Color' , this_color, 'LineStyle', ':', 'Marker', 'none', 'LineWidth', 0.25 )];
            this_hdl = [ this_hdl; plot(axes_hdl,  repmat([rowvec(aux_on(aux_complete_idx));rowvec(aux_off(aux_complete_idx))], 1, 2 ) , ...
                                    [ repmat(limits(1), 2, length(aux_complete_idx)) repmat(limits(2),2,length(aux_complete_idx)) ], 'Color' , this_color, 'LineStyle', ':', 'Marker', 'none', 'MarkerSize', 4, 'LineWidth', 0.25 )];
                                
                                
                                
%             aux_on_idx = find(~isnan(aux_on));
%             aux_off_idx = find(~isnan(aux_off));
%             aux_peak_idx = find(~isnan(aux_peak));
%             aux_complete_idx = find(~isnan(aux_on) & ~isnan(aux_off));
                            
        elseif( ann_graph_mode == kBackColourAnns )
            
            % wave start
            bOn = ~isnan(aux_on) & aux_on >= start_sample & aux_on <= end_sample;
            % wave end
            bOff = ~isnan(aux_off) & aux_off >= start_sample & aux_off <= end_sample;
            
            % wave peak
            bPeak = ~isnan(aux_peak) & aux_peak >= start_sample & aux_peak <= end_sample;
            aux_peak_idx = find(bPeak);
            this_hdl = [ this_hdl; plot(axes_hdl, repmat(rowvec(aux_peak(aux_peak_idx)), 2, 1 ), repmat(colvec(limits) + [-0.02; 0.02] * diff(limits), 1, length(aux_peak_idx)), 'Color' , this_color, 'LineStyle', ':', 'Marker', '^', 'MarkerSize', 4, 'LineWidth', 0.25 )];

            % wave conection between start-end
            bOnOff = bOn & bOff;
            aux_complete_idx = find(bOnOff);

            % normal sampled version
            aux_complete_idxx = arrayfun( @(a)( max(1, aux_on(a)):min(heasig.nsamp,aux_off(a)) ),aux_complete_idx, 'UniformOutput', false);

            % on-peak
            aux_complete_idx = find( ~bOnOff & bOn & bPeak);
            aux_complete_idxx = [ aux_complete_idxx; arrayfun( @(a)( aux_on(a):aux_peak(a) ),aux_complete_idx, 'UniformOutput', false) ];
            
            % peak-off
            aux_complete_idx = find(~bOnOff & bPeak & bOff);
            aux_complete_idxx = [ aux_complete_idxx; arrayfun( @(a)( aux_peak(a):aux_off(a) ),aux_complete_idx, 'UniformOutput', false) ];

            aux_offset = (start_sample - 1);
            %patch around the signal
%             this_hdl = [ this_hdl; cellfun( @(a)( patch( [a fliplr(a) ], [ (ECG(a-aux_offset, lead)* gains(lead) )- offsets(lead) + 0.5*yTextOffset ; flipud((ECG(a-aux_offset, lead)* gains(lead) )- offsets(lead)) - 0.5*yTextOffset ]', this_color, 'EdgeColor', 'none')), aux_complete_idxx) ];
            %box around the wave
            
            this_hdl = [ this_hdl; cellfun( @(a,b)( patch( [a(1) a(1) a(end) a(end) ] - aux_offset, [ limits(2) limits(1) limits(1) limits(2) ], this_color, 'EdgeColor', b)), aux_complete_idxx, this_edge_color) ];
            
        end
        
        uistack(this_hdl, 'bottom');        
    end
%--------------------------------------------------------------------------


%==========================================================================

    function [aux_val, str_unit_prefix] = microVoltsTransformer(aux_val)

        orig_aux_val = aux_val;
        aux_val = abs(aux_val);
        
        if(aux_val > 999)
            %milli
            str_unit_prefix = 'm';
            aux_val = orig_aux_val / 1e3;
        elseif(aux_val > 999999)
            str_unit_prefix = '';
            aux_val = orig_aux_val / 1e6;
        else
            %micro volts per default
            str_unit_prefix = '\\mu';
        end

    end
%--------------------------------------------------------------------------

    function disp_help()
        
        title_color = 'blue*';
        sub_title_color = 'magenta';
        
        disp_string_framed('*Blue', 'plot_ecg_strip help' );
        
        cprintf( title_color, 'Mouse actions:\n\n');
        
        cprintf( sub_title_color, '  Normal mode:\n' );
        
        fprintf(1, [... 
                '      single-click and holding LB : Activation Drag mode\n' ... 
                '      single-click and holding RB : Activation Rubber Band for region zooming\n' ... 
                '      single-click MB             : Activation ''Extend'' Zoom mode\n' ... 
                '      scroll wheel MB             : Activation Zoom mode\n' ... 
                '      double-click LB, RB, MB     : Reset to Original View\n' ... 
                    ] );
                       
        
        cprintf( sub_title_color, '  Magnifier mode:\n');
        
        fprintf(1, [... 
                '      single-click LB             : Not Used\n' ... 
                '      single-click RB             : Not Used\n' ... 
                '      single-click MB             : Reset Magnifier to Original View\n' ... 
                '      scroll MB                   : Change Magnifier Zoom\n' ... 
                '      double-click LB             : Increase Magnifier Size\n' ... 
                '      double-click RB             : Decrease Magnifier Size\n' ... 
                    ] );
                       
        fprintf(1, '\n');
        
        cprintf( title_color, 'Hotkeys in 2D mode:');

        fprintf(1, '\n\n');
        
        fprintf(1, [... 
                '      ''h''                         : Show help\n' ... 
                '      ''+''                         : Zoom plus\n' ... 
                '      ''-''                         : Zoom minus\n' ... 
                '      ''d''                         : Toggle the detail level of the annotations\n' ... 
                '      ''a''                         : Toggle the annotations graph mode\n' ... 
                '      ''0''                         : Set default axes (reset to original view)\n' ... 
                '      ''c''                         : On/Off pointer in crosshair mode\n' ... 
                '      ''g''                         : Change lead gain with scroll\n' ... 
                '      ''o''                         : Change lead offset with scroll\n' ... 
                '      ''x''                         : Zoom and drag works only for X axis\n' ... 
                '      ''y''                         : Zoom and drag works only for Y axis\n' ... 
                '      ''m''                         : If pressed and holding, Magnifier mode on\n' ... 
                '      ''p''                         : On/Off paper mode\n' ... 
                '      ''r''                         : Format of the exported file (PDF/PNG)\n' ... 
                '      ''s''                         : Export current view\n' ...                     ] );
                    ] );
                
    end

    function timer_stop_fcn(obj,event_obj)
         
%         if(bPreserveFix)
%             % never stop when editing
%             start(my_timer)
%         end
        
    end

    function timer_fcn(obj,event_obj)
        
        delete(findobj('Tag', 'title_efimero' ));

%         if(~bPreserveFix)
%             % allow edition of the closer wave
%             bFixedWave = false;
%         end
        
    end

    function update_title_efimero( strTitle, delay )
       
        delete(findobj('Tag', 'title_efimero' ))
        
        left_legend = plotXmin + 4*xTextOffset;
        bottom_legend = plotYmax - 2*yTextOffset;

        aux_hdl = text( left_legend , bottom_legend , strTitle, 'FontSize', 8, 'HorizontalAlignment', 'left', 'BackgroundColor', 'r' );

        set(aux_hdl, 'Tag', 'title_efimero')
            
        if( ~isinf(delay) && strcmpi(my_timer.Running, 'off') )
            my_timer.StartDelay = delay;
            start(my_timer)
        end
        
    end

    function my_closefcn(obj,event_obj)

        stop(my_timer)
        delete(my_timer)
        delete(fig_hdl)
        
    end

    function CheckAndDeleteHdl( this_hdl )

        if( ishandle(this_hdl) )
            delete(this_hdl);
        end
        
    end

end

%--------------------------------------------------------------------------
