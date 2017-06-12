%% (Internal) GUI for correcting QRS detections
% 
% Description: 
% This function implements a graphical user interface (GUI) to correct annotations
% possibly previous performed by an automatic algorithm. The idea is to
% easily visualize, compare and correct annotations. For this purpose the
% GUI presents several representations of the time evolution of the events
% (possibly but not necessarily heartbeats).
% 
% Arguments:
%     
%     +ECG: [numeric or cell]
%           
%           [numeric]: signal matrix of dimension [sig_length sig_size
%           repetitions_size] where:
%             - sig_length: time length in samples
%             - sig_size: number of ECG leads or number of signals.
%             - repetitions_size: number of repetitions of the same
%             signals. Typically used when time-synchronized events, like
%             heartbeats.  
%           
%           [cell]: cell array of length repetitions_size, where each cell
%           is (probably a time alligned event) a signal of dimension
%           [sig_length sig_size]
% 
%     +QRS_locations: [numeric] OPTIONAL. Default values enclosed in ()
%                 Synchronization sample. In ECG context, this values are
%                 the QRS fiducial point. (empty)
% 
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
%     +QRS_annotations: [cell] OPTIONAL. Default values enclosed in ()
%               Annotations to be included in the mosaic. The funcion
%               accepts 2 type of annotations: points and lines. 
% 
% Example:
% 
%    https://www.youtube.com/watch?v=qgWjvsvafVg&list=PLlD2eDv5CIe9sA2atmnb-DX48FIRG46z7&index=3
% 
% See also ECGtask_QRS_corrector
% 
% Limits and Known bugs:
%   Probably a lot :( ... but dont panic! send me feedback if you need help.
% 
% Author: Mariano Llamedo Soria (llamedom at frba.utn.edu.ar)
% Version: 0.1 beta
% Birthdate  : 16/2/2012
% Last update: 7/2/2014
% Copyright 2008-2015
% 
function ann_output = QRScorrector(varargin)

    %% Constants
    cached_filename = 'tmp_QRScorrector_cache.mat';

    cHeaderFieldNamesRequired = {'freq' 'nsamp' 'nsig' 'gain' 'adczero' };
    cAnnotationsFieldNamesRequired = {'time', 'qrs' };

    ann_output = [];

    %% Argument parsing

    %argument definition
    p = inputParser;   % Create instance of inputParser class.
    p.addParamValue('FilterSignals', [], @(x)( islogical(x)));
    p.addParamValue('BuildArtificial', false, @(x)( islogical(x)));
    p.addParamValue('recording', [], @(x)( ischar(x)));
    p.addParamValue('recording_path', [], @(x)( ischar(x)));
    p.addParamValue('recording_indexes', [], @(x)( isnumeric(x) && all(x > 0) ) );
    p.addParamValue('ECG', [], @(x)(isnumeric(x) || isa(x, 'ECGwrapper') ) );
    p.addParamValue('ECG_header', [], @(x)(isstruct(x)) );
    p.addParamValue('QRS_annotations', [], @(x)( (isnumeric(x) && all(x > 0)) || isstruct(x) ) );
    p.addParamValue('OutputVarName', 'payload', @(x)( ischar(x)));
    p.addParamValue('tmp_path', [], @(x)(ischar(x)) );
    p.addParamValue('Figure', [], @(x)(ishandle(x)) );

    try
        p.parse( varargin{:} );
    catch MyError
        rethrow(MyError);
    end

    hb_detail_window = 3;

    bFilter = p.Results.FilterSignals;
    bBuildArtificial = p.Results.BuildArtificial;
    recording = p.Results.recording;
    recording_path = p.Results.recording_path;
    recording_indexes = p.Results.recording_indexes;
    ECG = p.Results.ECG;
    ECG_header = p.Results.ECG_header;
    ECG_annotations = p.Results.QRS_annotations;
    tmp_path = p.Results.tmp_path;
    OutputVarName = p.Results.OutputVarName;
    fig_hdl = p.Results.Figure;

    title_efimero_hdl = [];
    RRserie_zoombars_hdl = [];
    RRserie_zoom_zoombars_hdl = [];
    Scatter_axes_hdl = [];
    RRserie_axes_hdl = [];
    RRserie_global_axes_hdl = [];
    RRserie_zoom_axes_hdl = [];
    ECG_axes_hdl = [];
    RRserie_selection_hdl = [];
    RRscatter_selection_hdl = [];
    copy_paste_buffer = [];
    
    % Dont know why this variable uses a lot of bytes to store at disk.
    clear p

    %% Argument parsing

    bLoadECG = true;

    ECG_w = [];

    start_idx = 1;
    win_size = 30; % minutes
    min_win_size = 1; % minutes
    max_win_size = 60; % minutes
    
    win_size_zoom = 5; % seconds
    min_win_size_zoom = 0.5; % seconds
    max_win_size_zoom = 30; % seconds
    
    
    if( ~isempty(ECG) )
        %% ECG already read

        if( isempty(ECG_annotations))
            disp_string_framed('*[1,0.5,0]', 'No annotations provided' );
        else
            if( isstruct(ECG_annotations) )
                if( any(isfield(ECG_annotations, cAnnotationsFieldNamesRequired )) )
                    % only one ann struct provided
                    ECG_struct.provided_anns = ECG_annotations;
                else
                    bFieldFound = false;
                    for fn = rowvec(fieldnames(ECG_annotations))
                        if( ~isempty(intersect( fieldnames(ECG_annotations.(fn{1})), cAnnotationsFieldNamesRequired )) )
                            bFieldFound = true;
                            break
                        end
                    end

                    if( bFieldFound )
                        ECG_struct = ECG_annotations;
                    else
                        error( 'QRScorrector:ArgCheck:InvalidAnnotations', disp_option_enumeration( 'Please provide the following fields in the annotations struct:',  cAnnotationsFieldNamesRequired) );
                    end                    

                end
            else
                ECG_struct.provided_anns.time = ECG_annotations;
            end
        end

        if( isa(ECG, 'ECGwrapper') )
            % parse ECGwrapper object
            ECG_w = ECG;

            if( isempty(ECG_header) )
                ECG_struct.header = ECG_w.ECG_header;
            end

            ECG_struct.signal = ECG_w.read_signal(start_idx, start_idx + round(win_size * 60 *ECG_struct.header.freq) + 10 * ECG_struct.header.freq );

        else
            
            if( isempty(ECG_header))
               error( 'QRScorrector:ArgCheck:InvalidHeader', 'Please provide the ECG header.\n\n' );
            else
                if( ~isfield(ECG_header, cHeaderFieldNamesRequired ) )
                    strAux = [ repmat(' + ', length(cHeaderFieldNamesRequired), 1) char(cHeaderFieldNamesRequired) repmat('\n', length(cHeaderFieldNamesRequired), 1 ) ];
                    error( 'QRScorrector:ArgCheck:InvalidHeader', ['Please provide the following fields in the header struct:\n ' rowvec(strAux') ] );
                end
            end

            ECG_struct.signal = ECG;
            ECG_struct.header = ECG_header;

        end        

        clear ECG ECG_header ECG_annotations
        
        recording_indexes = 1;
        rec_names.name = ECG_struct.header.recname;

        if( isempty(recording_path) )
            recording_path = [fileparts(mfilename('fullpath')) filesep 'tmp' filesep ];
        end   

        ratios = [];
        annotations_ranking = [];
        
        bLoadECG = false;
        
    elseif( ~isempty(recording) )

        recording_indexes = 1;
        if( isempty(recording_path) )
            [recording_path, rec_names.name, aux_ext] = fileparts(recording);
        else
            [~, rec_names.name, aux_ext] = fileparts(recording);
        end

        rec_names.name = [rec_names.name aux_ext];
        recording_path = [recording_path filesep];

    elseif( ~isempty(recording_path) )

        if( ~exist(recording_path, 'dir') )
            error('QRScorrector:ArgCheck:InvalidFolder', 'Folder does not exist.')
        end

        if( recording_path(end) ~= filesep )
            recording_path = [ recording_path filesep ];
        end

        rec_names = dir([recording_path '*.mat' ]);

        lrec_names = length(rec_names);

        if( lrec_names == 0)
            error([ 'No files found in ' recording_path ] )
        else

            if( isempty( recording_indexes ) )
                recording_indexes = 1:lrec_names;
            end

            lrec_names = length(recording_indexes);

            filename = [ recording_path 'tmp' filesep cached_filename ];

            bAux = (exist(filename, 'file') > 0);

            if( bAux )
                update_title_efimero(sprintf('Loading cached ratios from %s.\n', filename), 5 );
                
                aux_val2 = load(filename);
                ratios = aux_val2.ratios;
            end

            if( ~bAux || length(ratios) ~= lrec_names)

                ratios = nan(lrec_names,1);

                pb = progress_bar( 'Processing ratios', sprintf( '%d recordings found', lrec_names ) );

                pb.Loops2Do = lrec_names;

                for ii2 = recording_indexes
                    pb.start_loop();

                    pb.checkpoint(sprintf( 'Loading recording %d', ii2 ));
                    pepe = load([recording_path rec_names(ii2).name ], 'series_quality');
                    ratios(ii2) = max(pepe.series_quality.ratios);

                    pb.end_loop();
                end
                clear pb

                save(filename, 'ratios')
            end

            [ratios, worst_detections_idx] = sort(ratios);

            recording_indexes = recording_indexes(worst_detections_idx);

        end

    else
    %     strAux = help('QRScorrector'); %#ok<MCHLP>
        error( 'QRScorrector:ArgCheck:InvalidECGarg', 'Please provide an ECG recording as described in the documentation, help(''QRScorrector'') maybe could help you.\n' );
    end

    end_idx = start_idx + round((win_size * 60)*ECG_struct.header.freq);
    
    if( isempty(tmp_path) )
        tmp_path = recording_path;
    end

    %check path integrity.
    if(~exist(tmp_path, 'dir'))
        %try to create it
        if( ~mkdir(tmp_path) )
            error('QRScorrector:ArgCheck:InvalidPath', 'Invalid tmp_path. Please provide a valid path.\n' );
        end
    end
    
    major_tick_values_time = round([0.5 1 2 5 10 30 60]*ECG_struct.header.freq); % seconds
    major_tick_values_time = unique([major_tick_values_time major_tick_values_time*60 major_tick_values_time*60*60 ]);

    bAnnsEdited = false;
    bRecEdited = false;
    bSeries = false;
    
    all_annotations_selected = [];
    all_annotations_selected_serie_location = [];
    serie_location_mask = [];
    
    bFirstLoad = [];
    bSeriesChange = true;
    
    AnnNames = [];
    AnnNames_idx = [];
    all_annotations = [];
    estimated_labs = [];
    
    hb_idx = 1;
    lead_idx = 1;
    selected_hb_idx = [];
    RRscatter_hb_idx_hdl = [];   
    RRserie_hb_idx_hdl = [];   
    undo_buffer = [];
    undo_buffer_idx = 1;
    Pattern_hdl = [];
    ECG_limits = [];
    bLockECG = false;
    
    anns_under_edition = [];
    RRserie_filt = {};
    RRserie = {};
    RR_idx = {};
    anns_under_edition_idx = [];
    
    if( isfield(ECG_struct.header, 'btime') )
        aux_val2 = datevec(datenum(ECG_struct.header.btime, 'HH:MM:SS'));
        base_start_time = round((aux_val2(4) * 60 * 60 + aux_val2(5) * 60 + aux_val2(6)) * ECG_struct.header.freq);
    else
        base_start_time = 1;
    end
    
    ColorOrder = my_colormap( 12 );
    all_markers = {'+','o','*','.','x','s','d','^','v','<','>','p','h'};

    side_plot_hdl = [];

    fIsDragTimeAllowed = false;
    x_timeScroll_units = [];
    drag_timeScroll_start_x = [];
    k_timeScroll_units_pixel = [];
    bChangeWin = false;
    
    PrevStateWindowButtonMotionFcn = [];
    fIsDragAllowed = false;
    drag_start_x = [];
    drag_start_y = [];
    x_units = [];
    k_units_pixel = [];
    max_x_drag = [];
    min_y_drag = [];
    max_y_drag = [];
    
    filtro = [];
    rec_names = rec_names(recording_indexes,:);
    rec_idx = 1;
    recording_ratios = ratios;

    Xoffset = 0.09;
    Yoffset = 0.08;
    Xscale = 1.05;
    Yscale = 1.05;

    if( ECG_struct.header.nsamp > (win_size * 60 * ECG_struct.header.freq) ) 
        win_size = min(win_size, ECG_struct.header.nsamp / 60 / ECG_struct.header.freq / 3);
        size_y_RR_global = 0.13;
    else
        size_y_RR_global = 0;
    end
    
    annotation_list_control = [];
    leads_control = [];
    recordings_control = [];
    annotation_under_edition_label = [];

    maximized_size = [];
    
    my_timer = timer('TimerFcn',@timer_fcn, 'StartDelay', 10);
    
    if( isempty(fig_hdl) )
        fig_hdl = figure();
    else
        clf(fig_hdl,'reset')
    end

    set(fig_hdl, 'WindowButtonDownFcn', @WindowButtonDownCallback2D);            
    set(fig_hdl, 'WindowButtonUpFcn',   @WindowButtonUpCallback2D);            
    set(fig_hdl, 'WindowKeyPressFcn',   @KeyPress);            
    
    maximize(fig_hdl);
    maximized_size = get(fig_hdl, 'Position');

    set(fig_hdl, 'Position', [ maximized_size(3:4) maximized_size(3:4) ] .* [ 0.05 0.13 0.95 0.9] );
    
    set(fig_hdl,'CloseRequestFcn',@my_closefcn)

    DoRecording();


    function my_closefcn(obj,event_obj)

        if( bRecEdited )

           selection = questdlg('Unsaved results, edition will be lost. Exit anyway ?',...
              'Close Request Function',...
              'Yes','No','Yes'); 
           switch selection, 
              case 'Yes',
                    stop(my_timer)
                    delete(my_timer)
                    delete(fig_hdl)
              case 'No'
              return 
           end

        else
            stop(my_timer)
            delete(my_timer)
            delete(fig_hdl)
        end    
        
    end

    function DoRecording()

        [~, rec_name] = fileparts(rec_names(rec_idx).name);

        rec_path = [tmp_path rec_name ];

        if(bLoadECG)
            ECG_struct = load(rec_path);
            ECG_struct.header.recname = rec_name;
        end

        bAnnsEdited = false;
        bRecEdited = false;

        RRserie_filt = {[]};
        
        hb_idx = 1;
        lead_idx = 1;
        selected_hb_idx = [];
        RRscatter_hb_idx_hdl = [];   
        RRserie_hb_idx_hdl = [];   
        undo_buffer = [];
        undo_buffer_idx = 1;
        Pattern_hdl = [];

%         [~, ~, inv_dower_idx] = intersect({'I','II','V1','V2','V3','V4','V5','V6'}, cellstr(ECG_struct.header.desc));
%         if( length(inv_dower_idx) == 8 )
%             bUseDower = true;
%         else
%             inv_dower_idx = 1:ECG_struct.header.nsig;
%             bUseDower = false;
%         end

        ratios = [];
        AnnNames = [];
        all_annotations = [];
        estimated_labs = [];
        
        bFirstLoad = true;
        
        if( isempty(bFilter) )
            
            ECG_idx = get_ECG_idx_from_header(ECG_struct.header);

            if( isempty(ECG_idx) )
                bFilter = false;
                cprintf('[1,0.5,0]', disp_option_enumeration( 'No filter was applied since no ECG leads found. Lead description found:', cellstr(ECG_struct.header.desc) ) )
                fprintf(1, '\n')
            else
                bFilter = true;
            end
       
        end

        if( bFilter && isempty(filtro) )
            filtro = bandpass_filter_design( ECG_struct.header.freq );
        end
        
%         if( isfield(ECG_struct, 'noise_power') )
%             noise_power = ECG_struct.noise_power;
%         else
%             noise_power = [];
%         end

        if( isfield(ECG_struct, 'series_quality' ) ) 

            AnnNames = ECG_struct.series_quality.AnnNames;
            ratios = ECG_struct.series_quality.ratios;
            estimated_labs = ECG_struct.series_quality.estimated_labs;

            cant_anns = size(AnnNames,1);

            all_annotations = cell(cant_anns,1);
            for jj = 1:cant_anns
                all_annotations{jj} = ECG_struct.(AnnNames{jj,1}).(AnnNames{jj,2});
            end
            
        else
            
            calc_ratios();

        end

        if( isempty(ratios) )
            annotations_ranking = [];
        else
            
            [ratios, best_detections_idx] = sort(ratios, 'descend');

            aux_idx = find(cell2mat( cellfun(@(a)(~isempty(strfind(a, 'artificial'))), AnnNames(:,1), 'UniformOutput', false)));

            if( bBuildArtificial && isempty(aux_idx) && ~isempty(estimated_labs) )

                % generate artificial annotations combining K best annotations
                aux_idx = best_detections_idx(1:min(10, length(best_detections_idx)));
                artificial_annotations = combine_anns(all_annotations(aux_idx), estimated_labs(aux_idx), ECG_struct.header);

                if( ~isempty(artificial_annotations) )
                    for ii = length(artificial_annotations):-1:1
                        aux_str = ['artificial_' num2str(ii)];
                        ECG_struct.(aux_str) = artificial_annotations(ii);
                    end

                    calc_ratios();

                    [ratios, best_detections_idx] = sort(ratios, 'descend');
                end
            end

            aux_val = 1:length(ratios);

            [~, annotations_ranking] = sort(aux_val(best_detections_idx));

        end
        
        %     if( isfield(ECG_struct, 'corrected_anns' ) ) 
        %         AnnNames_idx = find(best_detections_idx == find(strcmp(AnnNames(:,1), 'corrected_anns') ));
        %         if( isfield(ECG_struct.corrected_anns, 'time' ) ) 
        %             anns_under_edition = unique(colvec(ECG_struct.corrected_anns.time));
        %         else
        %             anns_under_edition = unique(colvec(ECG_struct.corrected_anns));
        %         end
        %     else
        %         AnnNames_idx = 1;
        %         anns_under_edition = unique(round(colvec( ECG_struct.(AnnNames{AnnNames_idx,1}).(AnnNames{AnnNames_idx,2}) )));
        %     end

        AnnNames_idx = 1;
        
        if( isempty(AnnNames) )
            anns_under_edition = [];
        else
            aux_val = ECG_struct.(AnnNames{AnnNames_idx,1}).(AnnNames{AnnNames_idx,2});
            if( size(aux_val, 1) > 1 && size(aux_val, 2) > 1 )
                bSeries = true; 
                [anns_under_edition, aux_idx] = unique(round(colvec( aux_val(:,1) )));
                all_annotations_selected_serie_location = {aux_val(aux_idx,1) + round(aux_val(aux_idx,2) * ECG_struct.header.freq) };
            else
                anns_under_edition = unique(round(colvec( aux_val )));
            end
        end

        all_annotations_selected = {anns_under_edition};

        Redraw();

        set(annotation_list_control, 'Value', 1);
        set(leads_control, 'Value', 1);

    end

    function update_q_ratios(obj,event_obj) 


        calc_ratios();

        cant_anns = size(AnnNames,1);
        aux_str = repmat( ' - ',cant_anns,1);

        annotation_list_control = uicontrol( ... 
                      'style','listbox', ...
                      'units','normalized', ...
                      'string', [ char(cellstr(num2str((1:cant_anns)'))) aux_str char(AnnNames(:,1)) repmat( ' (',cant_anns,1) num2str(round(colvec(ratios * 1000))) aux_str num2str(colvec(annotations_ranking))  repmat( ')',cant_anns,1)  ] , ...
                      'position', [0.865 0.11 0.13 0.18] , ...
                      'min', 2, ...
                      'max', 4, ...
                      'callback', @ChangeAnnotationsSelected);

    end

    function calc_ratios()

        this_AnnNames = [];

        for fname = rowvec(fieldnames(ECG_struct))
            if( isfield(ECG_struct.(fname{1}), 'time') )
                this_AnnNames = [this_AnnNames; cellstr(fname{1}) cellstr('time')];
            end
            if( isfield(ECG_struct.(fname{1}), 'qrs') )
                this_AnnNames = [this_AnnNames; cellstr(fname{1}) cellstr('qrs')];
            end
        end

        if( isempty(AnnNames) )
            new_AnnNames_idx = 1:size(this_AnnNames,1);
        else
            [~, new_AnnNames_idx] = setdiff(this_AnnNames(:,1), AnnNames(:,1) );
        end
        
        cant_anns = length(new_AnnNames_idx);

        aux_val = cell(cant_anns,1);
        jj = 1;
        for ii = rowvec(new_AnnNames_idx)
            aux_val{jj} = ECG_struct.(this_AnnNames{ii,1}).(this_AnnNames{ii,2});
            jj = jj + 1;
        end

        if( isempty(aux_val) )
            this_ratios = [];
            this_estimated_labs = [];
        else
        %     ratios = CalcRRserieQuality(ECG_struct.signal, ECG_struct.header, all_annotations);
            [ this_ratios, this_estimated_labs ] = CalcRRserieRatio(aux_val, ECG_struct.header);
        end
        
        all_annotations = [all_annotations; aux_val ];
        ratios = [ratios; this_ratios ];
        estimated_labs = [estimated_labs; this_estimated_labs];
        AnnNames = [AnnNames; this_AnnNames(new_AnnNames_idx,:)];
        
        [ratios, best_detections_idx] = sort(ratios, 'descend');

        aux_val = 1:length(ratios);

        [~, annotations_ranking] = sort(aux_val(best_detections_idx));

        AnnNames = AnnNames(best_detections_idx,:);
        all_annotations = all_annotations(best_detections_idx);
        estimated_labs = estimated_labs(best_detections_idx);

    end

    function Redraw()
        
        % update only the one that can be edited

        if( bSeries )
            aux_val = all_annotations_selected_serie_location{1};
            % update the series values
            aux_RR = (aux_val - anns_under_edition) / ECG_struct.header.freq;

            all_annotations_selected(1) = {anns_under_edition};
            
            this_all_anns = all_annotations_selected_serie_location;
            aux_val = this_all_anns{1};
            aux_val(serie_location_mask) = nan;
            this_all_anns{1} = aux_val;
            
            RR_idx(1) = { find( isnan(anns_under_edition) | anns_under_edition >= start_idx &  anns_under_edition <= end_idx ) } ;
            
        else
            
            if( length(anns_under_edition) < 2 )
                aux_RR = [];
                this_all_anns = [];
            else
%                 aux_RR = colvec(diff(anns_under_edition));

                [aux_RR, aux_RR_filt ] = RR_calculation(anns_under_edition, ECG_struct.header.freq, RRserie_filt{1});               
                aux_RR = aux_RR * 1/ECG_struct.header.freq;          
                
                all_annotations_selected(1) = {anns_under_edition};
                this_all_anns = all_annotations_selected;
            end
            
            RR_idx(1) = { find( anns_under_edition >= start_idx &  anns_under_edition <= end_idx ) } ;
            
        end
        
        RRserie(1) = {aux_RR};
        RRserie_filt(1) = {aux_RR_filt};
        
        anns_under_edition_idx = RR_idx{1};
        
        if( isempty(aux_RR) || all(isnan(aux_RR)) )
            limits = [0.6 0.7];
        else        
            limits = prctile(aux_RR(anns_under_edition_idx), [5 97]);
        end
        
        aux_val = diff(limits);
        aux_val = max(0.01, 0.5*aux_val);
        limits(1) = max(min(aux_RR(anns_under_edition_idx)), limits(1) - aux_val);
        limits(2) = min(max(aux_RR(anns_under_edition_idx)), limits(2) + aux_val);

%         if( bLockScatter )
%             x_lims_scatter = get(Scatter_axes_hdl, 'Xlim');
%             y_lims_scatter = get(Scatter_axes_hdl, 'Ylim');
%         end

%         if( bLockRRserie )
%             x_lims_RRserie = get(RRserie_axes_hdl, 'Xlim');
%             y_lims_RRserie = get(RRserie_axes_hdl, 'Ylim');
%         end

%         figure(fig_hdl);

%         maximize(fig_hdl);
%         maximized_size = get(fig_hdl, 'Position');

%         set(fig_hdl, 'Position', [ maximized_size(3:4) maximized_size(3:4) ] .* [ 0.05 0.13 0.95 0.9] );

    %     set(fig_hdl, 'Toolbar', 'none');

        %% Axis de Poincar�

        if( isempty(Scatter_axes_hdl) )
            if(size_y_RR_global > 0 )
                Scatter_axes_hdl = axes('Position', ([0.55 0.41 0.355 0.52 ] - [ Xoffset Yoffset 0 0 ]) .* [ 1 1 Xscale Yscale ] - [0 0 0 size_y_RR_global], 'ColorOrder', ColorOrder, 'ButtonDownFcn',@inspect_scatter  ) ;
            else
                Scatter_axes_hdl = axes('Position', ([0.55 0.46 0.355 0.52 ] - [ Xoffset Yoffset 0 0 ]) .* [ 1 1 Xscale Yscale ], 'ColorOrder', ColorOrder, 'ButtonDownFcn',@inspect_scatter  ) ;
            end
        end
        aux_val_sc_op = get(Scatter_axes_hdl, 'OuterPosition');
        aux_val_sc = get(Scatter_axes_hdl, 'Position');
        
        cla(Scatter_axes_hdl)
        bRRempty = all(cellfun( @(a)(isempty(a)), RRserie));
        if( bRRempty )
            RRscatter_hdl = {};
            hold(Scatter_axes_hdl, 'on')
        else
            hold(Scatter_axes_hdl, 'on')
            RRscatter_hdl = cellfun( @(this_rr_serie, this_rr_idx, ii)( plot(Scatter_axes_hdl, this_rr_serie(this_rr_idx) , [this_rr_serie(this_rr_idx(2:end)); this_rr_serie(this_rr_idx(end)) ], 'Marker', all_markers{ii}, 'LineStyle', 'none', 'MarkerEdgeColor', ColorOrder(ii,:), 'Color', ColorOrder(ii,:), 'ButtonDownFcn',@inspect_scatter )  ), RRserie, RR_idx, num2cell((1:length(RRserie))'), 'UniformOutput', false );
        end

        if( ~bRRempty )
            aux_RR = RRserie{1};
            if( ~isempty(aux_RR) )
                aux_RR = aux_RR(anns_under_edition_idx);
                aux_RRnext = [aux_RR(2:end); aux_RR(end)];
                hold(Scatter_axes_hdl, 'on')
                RRscatter_selection_hdl = plot(Scatter_axes_hdl, aux_RR(selected_hb_idx), aux_RRnext(selected_hb_idx), 'xg', 'ButtonDownFcn',@inspect_scatter );
                plot(Scatter_axes_hdl, limits, limits, ':r', 'LineWidth', 0.5, 'ButtonDownFcn',@inspect_scatter );
                hold(Scatter_axes_hdl, 'off')
            end
        end    

%         if( bLockScatter )
%             set(Scatter_axes_hdl, 'Xlim', x_lims_scatter);
%             set(Scatter_axes_hdl, 'Ylim', y_lims_scatter);
%         else
            set(Scatter_axes_hdl, 'Xlim', limits);
            set(Scatter_axes_hdl, 'Ylim', limits);
%     %         zoom reset
%         end

        title(Scatter_axes_hdl, ['Poincar� plot interval of ' Seconds2HMS((end_idx-start_idx+1)/ECG_struct.header.freq) ]);
        xlabel(Scatter_axes_hdl, 'Current value')
        ylabel(Scatter_axes_hdl, 'Next value')

        %% Axis de RR serie

        if( isempty(RRserie_axes_hdl) )
            RRserie_axes_hdl = axes('Position', [(aux_val_sc(1) - 0.375)/2 (aux_val_sc(2)+aux_val_sc(4)/2) 0.375 aux_val_sc(4)/2 ] , 'ColorOrder', ColorOrder, 'ButtonDownFcn',@inspect_RRserie  );
        end

        cla(RRserie_axes_hdl)
        if( isempty(aux_RR) )
            RRserie_hdl = {};
            hold(RRserie_axes_hdl, 'on')
        else
            hold(RRserie_axes_hdl, 'on')
            RRserie_hdl = cellfun( @(this_anns, this_rr_serie, this_rr_idx, ii)( plot(RRserie_axes_hdl, this_anns(this_rr_idx) , this_rr_serie(this_rr_idx), 'Marker', all_markers{ii}, 'LineStyle', ':', 'MarkerEdgeColor', ColorOrder(ii,:), 'Color', ColorOrder(ii,:), 'ButtonDownFcn',@inspect_RRserie )  ), all_annotations_selected, RRserie, RR_idx, num2cell((1:length(RRserie))'), 'UniformOutput', false );
        end

        if( ~isempty(aux_RR) )
            aux_RR = RRserie{1};
            RRserie_selection_hdl = plot(RRserie_axes_hdl, anns_under_edition(anns_under_edition_idx(selected_hb_idx)), aux_RR(anns_under_edition_idx(selected_hb_idx)), 'og' );
        end

        if( isempty(anns_under_edition) )
            k_units_pixel = [];
        else
            prev_units = get(RRserie_axes_hdl, 'Units');
            set(RRserie_axes_hdl, 'Units', 'pixels');
            this_pos = get(RRserie_axes_hdl, 'Position');
            k_units_pixel = (anns_under_edition(anns_under_edition_idx(end)) - anns_under_edition(anns_under_edition_idx(1))) / this_pos(3);
            set(RRserie_axes_hdl, 'Units', prev_units);
        end
        
%         if( bLockRRserie )
%             xlim(RRserie_axes_hdl, x_lims_RRserie);
%             ylim(RRserie_axes_hdl, y_lims_RRserie);
%             this_xlims = x_lims_RRserie;
%             this_ylims = y_lims_RRserie;
%         else
            if( isempty(aux_RR) )
                x_max = 1;
                x_min = 0;
                this_ylims = [0.3 2];
            else
                ylim(RRserie_axes_hdl, limits);
                this_ylims = limits;
                x_max = max(anns_under_edition(anns_under_edition_idx));
                x_min = min(anns_under_edition(anns_under_edition_idx));
                x_range = x_max - x_min;
                xlim(RRserie_axes_hdl, [ x_min - 0.05 * x_range x_max + 0.05 * x_range ]);
            end
            this_xlims = [ x_min x_max ];
%         end
        
        % red box around
        this_xlims = this_xlims + 0.01*[ -diff(this_xlims) diff(this_xlims) ];
        this_ylims = this_ylims + 0.01*[ diff(this_ylims) -diff(this_ylims) ];
        set(fig_hdl,'CurrentAxes', RRserie_axes_hdl);
        aux_hdl = patch([this_xlims(1) this_xlims(1) this_xlims(2) this_xlims(2) this_xlims(1) ], [this_ylims(1) this_ylims(2) this_ylims(2) this_ylims(1) this_ylims(1)], [1 1 1], 'EdgeColor', [1 0 0], 'LineWidth', 1.5, 'ButtonDownFcn', @inspect_RRserie);
        set(aux_hdl, 'FaceColor', 'none')
        uistack(aux_hdl, 'bottom');
        hold(RRserie_axes_hdl, 'off')
        
        if( length(aux_RR) > 5 )
            cant_ticks = 5;
            aux_idx = round(linspace(x_min, x_max, cant_ticks));
            set(RRserie_axes_hdl, 'XTick', rowvec(aux_idx) );
            aux_str = cellstr(Seconds2HMS((aux_idx(1)+base_start_time-1)*1/ECG_struct.header.freq));    
            aux_str = [aux_str; cellstr(Seconds2HMS((aux_idx(2:(end-1)))*1/ECG_struct.header.freq))];    
            aux_str = [aux_str; cellstr(Seconds2HMS((aux_idx(end)+base_start_time-1)*1/ECG_struct.header.freq))];    
            set(RRserie_axes_hdl, 'XTickLabel', char(aux_str) );
        end
%         xlabel(RRserie_axes_hdl, 'Time')
        ylabel(RRserie_axes_hdl, 'Serie value')

        title(RRserie_axes_hdl, [ 'Interval of ' Seconds2HMS((end_idx-start_idx+1)/ECG_struct.header.freq)  ]);
        
        %% Axis de RR serie ZOOM

        if( isempty(RRserie_zoom_axes_hdl) )
            RRserie_zoom_axes_hdl = axes('Position', [(aux_val_sc(1) - 0.375)/2  aux_val_sc(2) 0.375 aux_val_sc(4)/2.3 ], 'ColorOrder', ColorOrder );
        end

        aux_val_RR_zoom = get(RRserie_zoom_axes_hdl, 'Position');
        max_x_drag = aux_val_RR_zoom(1) + aux_val_RR_zoom(3);
        min_y_drag = aux_val_RR_zoom(2);
        aux_val_RR = get(RRserie_axes_hdl, 'Position');
        max_y_drag = aux_val_RR(2)+aux_val_RR(4);
        
        if( isempty(RRserie) )
            cla(RRserie_zoom_axes_hdl)
        else
            UpdateRRserieZoom();
        end

        %% Axis de RR serie global

        if( size_y_RR_global > 0 )
            
            if( isempty(RRserie_global_axes_hdl) )
                RRserie_global_axes_hdl = axes('Position', [aux_val_RR(1) aux_val_sc(2)+1.13*aux_val_sc(4) aux_val_sc(1)+aux_val_sc(3)-aux_val_RR(1) size_y_RR_global ], 'ColorOrder', ColorOrder, 'ButtonDownFcn', @TimeOffsetClick );
            end
            
            set(fig_hdl, 'CurrentAxes', RRserie_global_axes_hdl);
            cla(RRserie_global_axes_hdl)
            
            if( isempty(aux_RR) )
                RRserie_hdl = [];
                hold(RRserie_global_axes_hdl, 'on')
            else
                hold(RRserie_global_axes_hdl, 'on')
                
                prev_units = get(RRserie_global_axes_hdl, 'Units');
                set(RRserie_global_axes_hdl, 'Units', 'pixels');

                % Downsample version: For efficient marks visualization and printing only
                target_res = 5; % samples/pixel

                axes_size = get(RRserie_global_axes_hdl, 'Position');
                nsamp_target = axes_size(3) * target_res;
                set(RRserie_global_axes_hdl, 'Units', prev_units);

                max_length = max(cellfun( @(this_rr_serie)( length(this_rr_serie)  ), RRserie ));
                
                down_factor = max(1, ceil( max_length / nsamp_target));
                
%                 RRserie_aux = cellfun( @(this_rr_serie)( resample(this_rr_serie, 1, down_factor, 60 )  ), RRserie, 'UniformOutput', false );
%                 cellfun( @(this_anns, this_rr_serie, ii)( plot(RRserie_global_axes_hdl, this_anns( round(linspace(1,length(this_anns),length(this_rr_serie))) )  , this_rr_serie, 'Marker', all_markers{ii}, 'LineStyle', ':', 'MarkerEdgeColor', ColorOrder(ii,:), 'Color', ColorOrder(ii,:), 'ButtonDownFcn', @TimeOffsetClick )  ), all_annotations_selected, RRserie_aux, num2cell((1:length(RRserie))') );

                % downsample this way to deal with NaN values.
                not_nan_idx = cellfun( @(this_rr_serie)( ~isnan(this_rr_serie) ), RRserie, 'UniformOutput', false );
                RRserie_aux = cellfun( @(this_anns, this_rr_serie, this_non_nan_idx)( interp1(this_anns(this_non_nan_idx), this_rr_serie(this_non_nan_idx), this_anns( round(linspace(1,length(this_anns), ceil(sum(this_non_nan_idx)/down_factor) ))), 'pchip')  ), all_annotations_selected, RRserie, not_nan_idx, 'UniformOutput', false );
                cellfun( @(this_anns, this_rr_serie, ii)( plot(RRserie_global_axes_hdl, this_anns( round(linspace(1,length(this_anns),length(this_rr_serie))) )  , this_rr_serie, 'Marker', all_markers{ii}, 'LineStyle', ':', 'MarkerEdgeColor', ColorOrder(ii,:), 'Color', ColorOrder(ii,:), 'ButtonDownFcn', @TimeOffsetClick )  ), all_annotations_selected, RRserie_aux, num2cell((1:length(RRserie))') );
                
                
                limits = prctile(aux_RR, [1 99]);
                ylim(RRserie_global_axes_hdl, limits);
                
%                 RRserie_zoombars_hdl = plot(RRserie_global_axes_hdl, repmat([ start_idx end_idx], 2, 1), repmat(limits,2,1)', 'LineWidth', 3, 'Color', 'r', 'ButtonDownFcn', @TimeOffsetClick );
                RRserie_zoombars_hdl = patch([start_idx start_idx end_idx end_idx start_idx ], [limits(1) limits(2) limits(2) limits(1) limits(1)], [241 183 171]/255, 'EdgeColor', [1 0 0], 'ButtonDownFcn', @TimeOffsetClick, 'LineWidth', 0.5);
                uistack(RRserie_zoombars_hdl, 'bottom');
                
                hold(RRserie_global_axes_hdl, 'off')

                max_end = min(cellfun( @(this_ann)( this_ann(end)  ), all_annotations_selected ));
                min_start = min(cellfun( @(this_ann)( this_ann(1)  ), all_annotations_selected ));
                max_range = max_end - min_start;
                
                xlim(RRserie_global_axes_hdl, [ min_start - 0.02*max_range max_end + 0.02*max_range ]);
                
                if( length(aux_RR) > 5 )
                    cant_ticks = 8;
                    
                    [~, major_tick_idx] = sort( abs(((anns_under_edition(end)-anns_under_edition(1))./major_tick_values_time - cant_ticks)));
                    major_tick = major_tick_values_time(major_tick_idx(1));
                    
                    aux_idx = round(1:major_tick:anns_under_edition(end));
                    set(RRserie_global_axes_hdl, 'XTick', rowvec(aux_idx) );
                    aux_str = cellstr(Seconds2HMS((aux_idx(1)+base_start_time-1)*1/ECG_struct.header.freq));    
                    aux_str = [aux_str; cellstr(Seconds2HMS((aux_idx(2:(end-1)))*1/ECG_struct.header.freq))];    
                    aux_str = [aux_str; cellstr(Seconds2HMS((aux_idx(end)+base_start_time-1)*1/ECG_struct.header.freq))];    
                    set(RRserie_global_axes_hdl, 'XTickLabel', char(aux_str) );
                end
%                 xlabel(RRserie_global_axes_hdl, 'Time')
                ylabel(RRserie_global_axes_hdl, 'Serie value')
                
                k_timeScroll_units_pixel = (anns_under_edition(end) - anns_under_edition(1)) / this_pos(3);
                
                title(RRserie_global_axes_hdl, [ 'Interval of ' Seconds2HMS(ECG_struct.header.nsamp/ECG_struct.header.freq) ]);
                
            end
            
        end

        
        %% Axis de ECG

        if( isempty(ECG_axes_hdl) )
            if(size_y_RR_global > 0 )
                ECG_axes_hdl = axes('Position', [aux_val_RR_zoom(1) 0.07*aux_val_RR_zoom(2) aux_val_sc(1)+aux_val_sc(3)-aux_val_RR_zoom(1) 0.75*aux_val_RR_zoom(2)], 'ColorOrder', ColorOrder  );
            else
                ECG_axes_hdl = axes('Position', [aux_val_RR_zoom(1) 0.07*aux_val_RR_zoom(2) aux_val_sc(1)+aux_val_sc(3)-aux_val_RR_zoom(1) 0.75*aux_val_RR_zoom(2)], 'ColorOrder', ColorOrder  );
            end
        end

        if( ~bLockECG )
            ECG_limits = [];
        end            
        
        if( isempty(anns_under_edition_idx) )
            [ECG_hdl, ECG_limits ] = plot_ecg_heartbeat(ECG_struct.signal, lead_idx, this_all_anns, start_idx, [] , hb_detail_window, ECG_struct.header, filtro, ECG_axes_hdl, ECG_limits);    
        else
            [ECG_hdl, ECG_limits ] = plot_ecg_heartbeat(ECG_struct.signal, lead_idx, this_all_anns, start_idx, anns_under_edition_idx(hb_idx) , hb_detail_window, ECG_struct.header, filtro, ECG_axes_hdl, ECG_limits);    
        end
        
        if( length(lead_idx) > 1 )
            aux_str = rowvec(colvec([repmat(',', length(lead_idx), 1) ECG_struct.header.desc(lead_idx,:) ]'));
            title(ECG_axes_hdl, ['Heartbeat ' num2str(hb_idx) ' : Leads ' aux_str(2:end) ] )
        else
            title(ECG_axes_hdl, ['Heartbeat ' num2str(hb_idx) ' : Lead ' ECG_struct.header.desc(lead_idx,:)] )
        end

%         xlabel(ECG_axes_hdl, 'Sample #')
        ylabel(ECG_axes_hdl, 'ECG')
    %     zoom reset

        %% Controls

        if( isempty(annotation_list_control) || bFirstLoad )

            bFirstLoad = false;
                      
            %% Recording list

            aux_str = repmat( ' - ', length(recording_indexes),1);

            recordings_control = uicontrol( ... 
                          'style','listbox', ...
                          'units','normalized', ...
                          'string', [ char(cellstr(num2str(colvec(recording_indexes)))) aux_str num2str(round(recording_ratios*1000)) aux_str char(rec_names.name) ] , ...
                          'position', [0.865 0.75 0.13 0.2] , ...
                          'min', 1, ...
                          'max', 1, ...
                          'Value', rec_idx, ...
                          'callback', @ChangeRecordingSelected);

            uicontrol( ...
                          'style','text', ...
                          'string', 'Available Recordings', ...
                          'units','normalized', ...
                          'position', [0.865 0.96 0.13 0.025] );  

            %% Signal list
            leads_control = uicontrol( ... 
                          'style','listbox', ...
                          'units','normalized', ...
                          'string', [ char(cellstr(num2str((1:ECG_struct.header.nsig)'))) repmat( ' - ',ECG_struct.header.nsig,1) ECG_struct.header.desc ] , ...
                          'position', [0.865 0.4 0.13 0.3] , ...
                          'min', 2, ...
                          'max', 4, ...
                          'Value', lead_idx, ...
                          'callback', @ChangeLeadsSelected);

            uicontrol( ...
                          'style','text', ...
                          'string', 'Available signals', ...
                          'units','normalized', ...
                          'position', [0.865 0.71 0.13 0.025] );        

            %% Annotation list

            if( isempty(AnnNames) )
                aux_str = 'manually_created';
                aux_label_str = 'Annotation under edition: manually-created (NaN)';
                AnnNames = [ {aux_str} {'time'} ];
                ECG_struct.(aux_str).time = [];
                
            else
                cant_anns = size(AnnNames,1);

                aux_str = repmat( ' - ',cant_anns,1);

                aux_str = [ char(cellstr(num2str((1:cant_anns)'))) aux_str char(AnnNames(:,1)) repmat( ' (',cant_anns,1) num2str(round(colvec(ratios * 1000))) aux_str num2str(colvec(annotations_ranking))  repmat( ')',cant_anns,1)  ];
                
                aux_label_str = [ 'Annotation under edition: ' char(AnnNames( AnnNames_idx ,1)) ' (' num2str(ratios(AnnNames_idx)) ')' ];
                
            end
            
            annotation_list_control = uicontrol( ... 
                          'style','listbox', ...
                          'units','normalized', ...
                          'string', aux_str, ...
                          'position', [0.865 0.11 0.13 0.18] , ...
                          'min', 2, ...
                          'max', 4, ...
                          'callback', @ChangeAnnotationsSelected);

            uicontrol( ...
                          'style','text', ...
                          'string', 'Annotations available', ...
                          'units','normalized', ...
                          'position', [0.865 0.30 0.13 0.025] );

            annotation_under_edition_label = uicontrol( ...
                          'style','text', ...
                          'string', aux_label_str, ...
                          'units','normalized', ...
                          'position', [0.865 0.35 0.13 0.05] );

            uicontrol( ... 
                            'style','pushbutton', ...
                            'string', 'Delete Annotations', ...
                            'units','normalized', ...
                            'position', [0.865 0.07 0.063 0.03], ...
                            'callback', @DeleteAnnotations);

            uicontrol( ... 
                            'style','pushbutton', ...
                            'string', 'Upd Q ratios', ...
                            'units','normalized', ...
                            'position', [0.93 0.07 0.063 0.03], ...
                            'callback', @update_q_ratios);

            uicontrol( ... 
                            'style','pushbutton', ...
                            'string', 'Merge Annotations', ...
                            'units','normalized', ...
                            'position', [0.865 0.04 0.063 0.03], ...
                            'callback', @MergeAnnotations);

        end

%         if( ~isempty(RRscatter_hdl) )
%             set(RRscatter_hdl,'ButtonDownFcn',@inspect_scatter);            
%         end
%         set(Scatter_axes_hdl,'ButtonDownFcn',@inspect_scatter);            
        
%         if( ~isempty(RRserie_hdl) )
%             set(RRserie_hdl, 'ButtonDownFcn',@inspect_RRserie);            
%         end

        
%         set(RRserie_axes_hdl,'ButtonDownFcn',@inspect_RRserie);
        
        cellfun(@(a)( set(a,'ButtonDownFcn',@inspect_ECG)), ECG_hdl);            
        set(ECG_axes_hdl,'ButtonDownFcn',@inspect_ECG);            

    end

    function DragMouseBegin()
        %DragMouseBegin begin draging
        
        if ( ~fIsDragAllowed )
            
            [drag_start_x, drag_start_y] = GetCursorCoordOnWindow();
           
            fIsDragAllowed = true;
            PrevStateWindowButtonMotionFcn = get(fig_hdl, 'WindowButtonMotionFcn');
            set(fig_hdl, 'WindowButtonMotionFcn', @WindowButtonMotionCallback2D);
            
%             fprintf(1, 'on\n');
            
        end
    end

    function DragMouseEnd()
        %DragMouseEnd end draging

        if fIsDragAllowed
            fIsDragAllowed = false;
            
            set(fig_hdl, 'WindowButtonMotionFcn', PrevStateWindowButtonMotionFcn);
%             fprintf(1, 'off\n');

            if ( ~fIsDragTimeAllowed )

% esto lo comenté porque generaba una doble entrada a      plot_ecg_heartbeat que no entendí para que estaba puesto           
% 
%                 if( bSeries )
%                     this_all_anns = all_annotations_selected_serie_location;
%                     aux_val = this_all_anns{1};
%                     aux_val(serie_location_mask) = nan;
%                     this_all_anns{1} = aux_val;
%                 else
%                     this_all_anns = all_annotations_selected;
%                 end
%                 
%                 if( ~bLockECG )
%                     ECG_limits = [];
%                 end            
%                 
%                 [ECG_hdl, ECG_limits ] = plot_ecg_heartbeat(ECG_struct.signal, lead_idx, this_all_anns, start_idx, anns_under_edition_idx(hb_idx) , hb_detail_window, ECG_struct.header, filtro, ECG_axes_hdl, ECG_limits);
% 
%                 cellfun(@(a)( set(a,'ButtonDownFcn',@inspect_ECG)), ECG_hdl);            

            end
        end
    end

    function [xp, yp ] = GetCursorCoordOnWindow()
        %GetCursorCoordOnWindow
        
        dfltUnits = get(fig_hdl, 'Units');
        
        set(fig_hdl, 'Units', 'pixels');
        
        crd = get(fig_hdl, 'CurrentPoint');
        xp = crd(1); 
        yp = crd(2);
        
        set(fig_hdl, 'Units', dfltUnits);
    end

    function WindowButtonMotionCallback2D(src, evnt)  %#ok
        %WindowButtonMotionCallback2D
                
        % RR serie part
        if( fIsDragAllowed && ~isempty(x_units) )
            
            [drag_x, drag_y] = GetCursorCoordOnWindow();
            
            deltax = drag_x - drag_start_x;
            
            if( bChangeWin )
                
                win_size_zoom = min(max_win_size_zoom, max( min_win_size_zoom, win_size_zoom + (( deltax * 0.1 ) * k_units_pixel / ECG_struct.header.freq ) ));
                
%                 set(RRserie_zoombars_hdl, 'Xdata', [start_idx start_idx end_idx end_idx start_idx ]);
            
                update_title_efimero( sprintf('%s', Seconds2HMS(win_size_zoom)), 5 );

                UpdateRRserieZoom();
                
            else
            
                this_x_units = x_units + ( deltax * 0.1 ) * k_units_pixel;

                [~, aux_val] = sort( abs( this_x_units - anns_under_edition(anns_under_edition_idx)) );
                aux_val = aux_val(1);

                hb_idx = max(1, min(length(anns_under_edition_idx), aux_val ));        

    %             disp(hb_idx)

                UpdateRRserieZoom();

%                 if( ishandle(RRserie_hb_idx_hdl) )
%                     delete(RRserie_hb_idx_hdl);
%                 end
                aux_RR = RRserie{1};

%                 hold(RRserie_axes_hdl, 'on')
%                 RRserie_hb_idx_hdl = plot(RRserie_axes_hdl, anns_under_edition(anns_under_edition_idx(hb_idx)), aux_RR(anns_under_edition_idx(hb_idx)), 'or' );
%                 hold( RRserie_axes_hdl, 'off')        

                set(RRserie_hb_idx_hdl, 'Xdata', anns_under_edition(anns_under_edition_idx(hb_idx)) );
                set(RRserie_hb_idx_hdl, 'Ydata', aux_RR(anns_under_edition_idx(hb_idx)) );
                
%                 [ECG_hdl, ECG_limits ] = plot_ecg_heartbeat(ECG_struct.signal, lead_idx, all_annotations_selected, start_idx, anns_under_edition_idx(hb_idx) , hb_detail_window, ECG_struct.header, filtro, ECG_axes_hdl);    
            
            end
            
        else
            
%             fprintf(1, 'wait\n');
            
        end
    
        % RR global part
        if( fIsDragTimeAllowed )
            
            [drag_x, drag_y] = GetCursorCoordOnWindow();
            
            deltax = drag_x - drag_timeScroll_start_x;
            
            UpdateStartX(deltax);
            
        end
        
        
    
    end

    function UpdateStartX(deltax)

        if( bChangeWin )
            win_size = min(max_win_size, max( min_win_size, win_size + (( deltax * 0.1 ) * k_timeScroll_units_pixel * 1/ECG_struct.header.freq/60 ) ));
            update_title_efimero( sprintf('%s', Seconds2HMS(win_size*60)), 5);
        else
            start_idx = max(1, round(min( ECG_struct.header.nsamp - (win_size * 60 * ECG_struct.header.freq), x_timeScroll_units + ( deltax ) * k_timeScroll_units_pixel)) );
        end

        end_idx = max((min_win_size * 60 * ECG_struct.header.freq), min( ECG_struct.header.nsamp, start_idx + round((win_size * 60)*ECG_struct.header.freq)));

        set(RRserie_zoombars_hdl, 'Xdata', [start_idx start_idx end_idx end_idx start_idx ]);

        aux_idx = get(RRserie_axes_hdl, 'XTick' );
        aux_str = repmat({''},length(aux_idx),1);
        aux_str(1) = {Seconds2HMS((start_idx+base_start_time-1)*1/ECG_struct.header.freq)};  
        aux_str(end) = {Seconds2HMS((end_idx+base_start_time-1)*1/ECG_struct.header.freq)};  
        set(RRserie_axes_hdl, 'XTickLabel', char(aux_str) );       
        
    end


    function inspect_scatter(obj,event_obj)


        if( strcmp(get(fig_hdl,'SelectionType'),'alt'))
            % Delete annotation

            aux_RR = RRserie{1};
            aux_RR = aux_RR(anns_under_edition_idx);
            
            point = get(gca,'CurrentPoint');

            [~, hb_idx] = min(sum( bsxfun(@minus, point(1,1:2), [aux_RR, [aux_RR(2:end); aux_RR(end)]]).^2,2));

            PushUndoAction();

            if( bSeries )
                anns_under_edition(hb_idx) = nan;
            else
                anns_under_edition(hb_idx) = [];
                aux_rr_filt = RRserie_filt{1};
                aux_rr_filt(hb_idx) = [];
                RRserie_filt{1} = aux_rr_filt;
                
            end

            selected_hb_idx = [];

            Redraw();

            bAnnsEdited = true;
            bRecEdited = true;

        else

            if (strcmp(get(fig_hdl,'SelectionType'),'extend'))
                point1 = get(gca,'CurrentPoint');    % button down detected
                rbbox;
                point2 = get(gca,'CurrentPoint');    % button up detected            

                if( ~isempty(selected_hb_idx) )
                    delete(RRserie_selection_hdl)
                    delete(RRscatter_selection_hdl)
                end

                xlims = sort([ point1(1,1) point2(1,1) ]);
                ylims = sort([ point1(1,2) point2(1,2) ]);
                aux_RRcurr = RRserie{1};
                aux_RRcurr = aux_RRcurr(anns_under_edition_idx);
                aux_RRnext = [aux_RRcurr(2:end); aux_RRcurr(end)];
                selected_hb_idx = find( aux_RRcurr > xlims(1) & aux_RRcurr < xlims(2) & aux_RRnext > ylims(1) & aux_RRnext < ylims(2) );

                update_title_efimero([num2str(length(selected_hb_idx)) ' heartbeats selected.'], 5 );
                
                hold(Scatter_axes_hdl, 'on')
                RRnext = [RRserie(2:end); RRserie(end)];
                RRscatter_selection_hdl = plot(Scatter_axes_hdl, aux_RRcurr(anns_under_edition_idx(selected_hb_idx)), aux_RRnext(anns_under_edition_idx(selected_hb_idx)), 'xg' );
                hold(Scatter_axes_hdl, 'off')

                hold(RRserie_axes_hdl, 'on')
                RRserie_selection_hdl = plot(RRserie_axes_hdl, anns_under_edition(anns_under_edition_idx(selected_hb_idx)) , aux_RRcurr(selected_hb_idx), 'og' );
                hold(RRserie_axes_hdl, 'off')

                min_hb_idx = min(selected_hb_idx);
                max_hb_idx = max(selected_hb_idx);
                cant_hb_idx = max(selected_hb_idx) - min_hb_idx + 1;

                if( (anns_under_edition(anns_under_edition_idx(max_hb_idx)) - anns_under_edition(anns_under_edition_idx(min_hb_idx))) <= (10*ECG_struct.header.freq) )

                    if( bSeries )
                        this_all_anns = all_annotations_selected_serie_location;
                        aux_val = this_all_anns{1};
                        aux_val(serie_location_mask) = nan;
                        this_all_anns{1} = aux_val;
                    else
                        this_all_anns = all_annotations_selected;
                    end

                    if( ~bLockECG )
                        ECG_limits = [];
                    end            
                    
                    [ECG_hdl, ECG_limits ] = plot_ecg_heartbeat(ECG_struct.signal, lead_idx, this_all_anns, start_idx, anns_under_edition_idx(hb_idx) , hb_detail_window, ECG_struct.header, filtro, ECG_axes_hdl, ECG_limits);
                    cellfun(@(a)( set(a,'ButtonDownFcn',@inspect_ECG)), ECG_hdl);            
                end
                
                if( length(lead_idx) > 1 )
                    aux_str = rowvec(colvec([repmat(',', length(lead_idx), 1) ECG_struct.header.desc(lead_idx,:) ]'));
                    title(ECG_axes_hdl, ['Heartbeat ' num2str(hb_idx) ' : Leads ' aux_str(2:end) ] )
                else
                    title(ECG_axes_hdl, ['Heartbeat ' num2str(min_hb_idx) ' : Lead ' ECG_struct.header.desc(lead_idx,:)] )
                end

    %             zoom reset

                set(ECG_axes_hdl,'ButtonDownFcn',@inspect_ECG);            

                % Zoom RR serie
                UpdateRRserieZoom();

            else
                % Show ECG
                point = get(gca,'CurrentPoint');

                if( ishandle(RRscatter_hb_idx_hdl) )
                    delete(RRscatter_hb_idx_hdl)
                end
                if( ishandle(RRserie_hb_idx_hdl) )
                    delete(RRserie_hb_idx_hdl)
                end

                aux_RRcurr = RRserie{1};
                aux_RRcurr = aux_RRcurr(anns_under_edition_idx);
                
                [~, hb_idx] = min(sum( bsxfun(@minus, point(1,1:2), [aux_RRcurr(1:end-1), aux_RRcurr(2:end)]).^2,2));

                hold(Scatter_axes_hdl, 'on')
                RRscatter_hb_idx_hdl = plot(Scatter_axes_hdl, aux_RRcurr(hb_idx), aux_RRcurr(hb_idx+1), 'or' );
                hold(Scatter_axes_hdl, 'off')

                hold(RRserie_axes_hdl, 'on')
                RRserie_hb_idx_hdl = plot(RRserie_axes_hdl, anns_under_edition(anns_under_edition_idx(hb_idx)) , aux_RRcurr(hb_idx), 'or' );
                hold(RRserie_axes_hdl, 'off')

                if( bSeries )
                    this_all_anns = all_annotations_selected_serie_location;
                    aux_val = this_all_anns{1};
                    aux_val(serie_location_mask) = nan;
                    this_all_anns{1} = aux_val;
                else
                    this_all_anns = all_annotations_selected;
                end
                
                if( ~bLockECG )
                    ECG_limits = [];
                end            
                
                [ECG_hdl, ECG_limits ] = plot_ecg_heartbeat(ECG_struct.signal, lead_idx, this_all_anns, start_idx, anns_under_edition_idx(hb_idx) , hb_detail_window, ECG_struct.header, filtro, ECG_axes_hdl, ECG_limits);

                if( length(lead_idx) > 1 )
                    aux_str = rowvec(colvec([repmat(',', length(lead_idx), 1) ECG_struct.header.desc(lead_idx,:) ]'));
                    title(ECG_axes_hdl, ['Heartbeat ' num2str(hb_idx) ' : Leads ' aux_str(2:end) ] )
                else
                    title(ECG_axes_hdl, ['Heartbeat ' num2str(hb_idx) ' : Lead ' ECG_struct.header.desc(lead_idx,:)] )
                end

    %             zoom reset

                cellfun(@(a)( set(a,'ButtonDownFcn',@inspect_ECG)), ECG_hdl);            

                set(ECG_axes_hdl,'ButtonDownFcn',@inspect_ECG);            

                % Zoom RR serie
                UpdateRRserieZoom();

            end

        end


    end
    
    function inspect_ECG(obj,event_obj)


    %     disp(get(fig_hdl,'SelectionType'))

        if (strcmp(get(fig_hdl,'SelectionType'),'alt'))
            % Delete annotation

            point = get(gca,'CurrentPoint');

            point(1) = point(1) + start_idx - 1;
            
            [~, hb_idx] = min(abs( point(1) - anns_under_edition(anns_under_edition_idx)));

            PushUndoAction();

            if( bSeries )
                serie_location_mask(anns_under_edition_idx(hb_idx)) = true;
            else
                anns_under_edition(anns_under_edition_idx(hb_idx)) = [];
                aux_rr_filt = RRserie_filt{1};
                aux_rr_filt(anns_under_edition_idx(hb_idx)) = [];
                RRserie_filt{1} = aux_rr_filt;
            end
            
            selected_hb_idx( selected_hb_idx > length(anns_under_edition_idx) | selected_hb_idx == hb_idx ) = [];

            if( isempty(anns_under_edition) )
                hb_idx = [];
            else
                hb_idx = max(1, hb_idx - 1); 
            end

            Redraw();

            bAnnsEdited = true;
            bRecEdited = true;


        elseif (strcmp(get(fig_hdl,'SelectionType'),'normal'))

            % Show ECG

            point = get(gca,'CurrentPoint');

            point(1) = point(1) + start_idx - 1;
            
            PushUndoAction();

            if( bSeries )
                % only modify the location, not adding allowed.
                aux_val = all_annotations_selected_serie_location{1};
                
                [~, hb_idx] = min(abs( point(1) - aux_val(anns_under_edition_idx)));

                % refine the point, and enable it.
                aux_val(anns_under_edition_idx(hb_idx)) = round(point(1));
                serie_location_mask(anns_under_edition_idx(hb_idx)) = false;

                all_annotations_selected_serie_location{1} = aux_val;
                
            else
                % add a regular event when the anns are not series.
                aux_val = round(point(1));
                if( aux_val < start_idx || aux_val > end_idx )
                    update_title_efimero('Annotation out of red box bounds', 5 );
                    return
                end
                anns_under_edition = sort([anns_under_edition; aux_val ]);
                anns_under_edition_idx = find( anns_under_edition >= start_idx &  anns_under_edition <= end_idx );
                selected_hb_idx = find(aux_val == anns_under_edition(anns_under_edition_idx), 1);
                hb_idx = selected_hb_idx;
                aux_rr_filt = RRserie_filt{1};
                aux_rr_filt = aux_rr_filt([1:hb_idx hb_idx:end]);
                RRserie_filt{1} = aux_rr_filt;
                
            end
    
            Redraw();

            bAnnsEdited = true;
            bRecEdited = true;


        elseif (strcmp(get(fig_hdl,'SelectionType'),'extend'))

            % Show ECG

            point = get(gca,'CurrentPoint');
            aux_val = round(point(1));

            aux_val = aux_val + start_idx - 1;
            
            if( isempty(side_plot_hdl) || ~ishandle(side_plot_hdl) )
                side_plot_hdl = figure();
                set(side_plot_hdl, 'Position', [ maximized_size(3:4) maximized_size(3:4) ] .* [ 0 0 0.95 0.9] );
            else
                figure(side_plot_hdl);
            end

            plot_ecg_strip(ECG_struct.signal, ...
                            'ECG_header', ECG_struct.header, ...
                            'QRS_locations', anns_under_edition, ...
                            'Start_time', aux_val/ECG_struct.header.freq - win_size_zoom/2 -  1 , ...
                            'Figure_handle', side_plot_hdl , ...
                            'End_time', aux_val/ECG_struct.header.freq + win_size_zoom/2 + 1 );

            side_plot_hdl = gcf();
            
            figure(fig_hdl);

        end

    end
    
    function KeyPress(obj,event_obj)

        if (strcmp(event_obj.Key,'c') && strcmp(event_obj.Modifier,'control'))
            % copy selection
            if( ~isempty(selected_hb_idx) )
                copy_paste_buffer = anns_under_edition(selected_hb_idx);
                update_title_efimero('Selection copied', 5 );
                
            end

        elseif (strcmp(event_obj.Key,'v') && strcmp(event_obj.Modifier,'control'))
            % paste selection
            if( ~isempty(copy_paste_buffer) )
                bAnnsEdited = true;
                bRecEdited = true;
                PushUndoAction();
                anns_under_edition = sort( unique([anns_under_edition; colvec(copy_paste_buffer) ]) );
                selected_hb_idx = [];
                RRserie_filt = {[]};

                Redraw();
            end

        elseif (strcmp(event_obj.Key,'delete'))
            % delete selection
            if( ~isempty(selected_hb_idx) )
                bAnnsEdited = true;
                bRecEdited = true;
                PushUndoAction();

                aux_idx = [find(anns_under_edition < anns_under_edition(anns_under_edition_idx(selected_hb_idx(1))),1, 'last') find( anns_under_edition > anns_under_edition(anns_under_edition_idx(selected_hb_idx(end))), 1, 'first') ];

                if( isempty(aux_idx) ) 
                    aux_idx = 1;
                else
                    aux_idx = aux_idx(1);
                end
                
                aux_val = anns_under_edition(aux_idx);
                
                if( bSeries )
                    anns_under_edition(anns_under_edition_idx(selected_hb_idx)) = nan;
                else
                    anns_under_edition(anns_under_edition_idx(selected_hb_idx)) = [];
                end

                anns_under_edition_idx = find( anns_under_edition >= start_idx &  anns_under_edition <= end_idx );
                hb_idx = find(anns_under_edition(anns_under_edition_idx) == aux_val);
                selected_hb_idx = [];
                
                aux_rr_filt = RRserie_filt{1};
                aux_rr_filt(anns_under_edition_idx(selected_hb_idx)) = [];
                RRserie_filt{1} = aux_rr_filt;
                
                Redraw();
            end

        elseif (strcmp(event_obj.Key,'y') && ~isempty(event_obj.Modifier) && strcmp(event_obj.Modifier,'control'))
            % redo last change
            if( undo_buffer_idx < length(undo_buffer) )
                undo_buffer_idx = undo_buffer_idx + 1;
                anns_under_edition = undo_buffer{undo_buffer_idx};
                selected_hb_idx = [];
                RRserie_filt = {[]};
                Redraw();
            end

        elseif (strcmp(event_obj.Key,'z') && strcmp(event_obj.Modifier,'control'))
            % undo last change
            if( undo_buffer_idx > 1 )
                undo_buffer_idx = undo_buffer_idx - 1;
                anns_under_edition = undo_buffer{undo_buffer_idx};
                selected_hb_idx = [];
                RRserie_filt = {[]};
                
                Redraw();
            end

        elseif (strcmp(event_obj.Key,'rightarrow') && ~isempty(event_obj.Modifier) && strcmp(event_obj.Modifier,'control'))


            if(bRecEdited)

                update_annotations();

                if( bLoadECG )
                    update_title_efimero('Saving data ...', 5 );
                    save(rec_path, '-struct', 'ECG_struct');        
                    update_title_efimero(['Saved ' rec_path], 5 );

                    if( rec_idx >= 1 && rec_idx < length(recording_indexes) )
                        rec_idx = rec_idx + 1;
                    else
                        rec_idx = 1;
                    end

                    set(recordings_control, 'Value', rec_idx );

                    DoRecording();                

                else

                    if( isfield(ECG_struct, 'series_quality' ) ) 
                        anns_struct.series_quality = ECG_struct.series_quality;
                    end
                    for ii = 1:size(AnnNames,1)
                        anns_struct.(AnnNames{ii,1}) = ECG_struct.(AnnNames{ii,1});
                    end
                    assignin( 'caller', OutputVarName, anns_struct );
                    update_title_efimero(sprintf('Saving ''%s'' variable in caller workspace', OutputVarName), 5 );

                    bAnnsEdited = false;
                    bRecEdited = false;

                end
            end


        elseif (strcmp(event_obj.Key,'leftarrow') && ~isempty(event_obj.Modifier) && strcmp(event_obj.Modifier,'control'))

            if(bRecEdited)

                update_annotations();

                if( bLoadECG )
                    update_title_efimero('Saving data ...', 5 );
                    save(rec_path, '-struct', 'ECG_struct');        
                    update_title_efimero(['Saved ' rec_path], 5 );

                    if( rec_idx > 1 && rec_idx <= length(recording_indexes) )
                        rec_idx = rec_idx - 1;
                    else
                        rec_idx = length(recording_indexes);
                    end

                    set(recordings_control, 'Value', rec_idx );

                    DoRecording();                

                else
                    if( isfield(ECG_struct, 'series_quality' ) ) 
                        anns_struct.series_quality = ECG_struct.series_quality;
                    end
                    for ii = 1:size(AnnNames,1)
                        anns_struct.(AnnNames{ii,1}) = ECG_struct.(AnnNames{ii,1});
                    end
                    assignin( 'caller', OutputVarName, anns_struct );
                    update_title_efimero(sprintf('Saving ''%s'' variable in caller workspace', OutputVarName), 5 );

                    bAnnsEdited = false;
                    bRecEdited = false;

                end
            end

        elseif (strcmp(event_obj.Key, 'l') )
            % lock ECG Y axis
            
            bLockECG = ~bLockECG;
            
            if( bLockECG )
                update_title_efimero('ECG Y axis locked', 5 );
            else
                update_title_efimero('ECG Y axis unlocked', 5 );
            end            

        elseif (strcmp(event_obj.Key,'p') )

            update_title_efimero('Click and drag the time interval where the pattern is ...', 10 );
            
            set(fig_hdl, 'WindowButtonDownFcn', []);    
            set(ECG_axes_hdl,'ButtonDownFcn', []);            
            
            set(obj,'CurrentAxes', ECG_axes_hdl);
            waitforbuttonpress;
            
            try
                SearchPattern();
                set(fig_hdl, 'WindowButtonDownFcn', @WindowButtonDownCallback2D);            
                set(ECG_axes_hdl,'ButtonDownFcn',@inspect_ECG);            
                
            catch ME
                
                set(fig_hdl, 'WindowButtonDownFcn', @WindowButtonDownCallback2D);            
                set(ECG_axes_hdl,'ButtonDownFcn',@inspect_ECG);            
                rethrow(ME);
            end
            

        elseif ( ~isempty(event_obj.Modifier) && strcmp(event_obj.Key,'g') && strcmp(event_obj.Modifier,'control'))

            if(bRecEdited)
            
                update_annotations();
                
                if( bLoadECG )
                    update_title_efimero('Saving data ...', 5 );
                    save(rec_path, '-struct', 'ECG_struct');        
                    update_title_efimero(['Saved ' rec_path], 5 );
                else

                    if( isfield(ECG_struct, 'series_quality' ) ) 
                        anns_struct.series_quality = ECG_struct.series_quality;
                    end
                    for ii = 1:size(AnnNames,1)
                        anns_struct.(AnnNames{ii,1}) = ECG_struct.(AnnNames{ii,1});
                    end
                    assignin( 'caller', OutputVarName, anns_struct );
                    update_title_efimero(sprintf('Saving ''%s'' variable in caller workspace', OutputVarName), 5 );

                    bAnnsEdited = false;
                    bRecEdited = false;
                end
                
            end

        elseif (strcmp(event_obj.Key,'t'))
            %% Toggle ECG lead
            if( length(lead_idx) > 1 )
                lead_idx = 1;
            else
                if( lead_idx < ECG_struct.header.nsig )
                    lead_idx = lead_idx + 1;
                else
                    lead_idx = 1:ECG_struct.header.nsig;
                end
            end
            cla(ECG_axes_hdl);

            if( bSeries )
                this_all_anns = all_annotations_selected_serie_location;
                aux_val = this_all_anns{1};
                aux_val(serie_location_mask) = nan;
                this_all_anns{1} = aux_val;
            else
                this_all_anns = all_annotations_selected;
            end

            if( ~bLockECG )
                ECG_limits = [];
            end            
            
            [ECG_hdl, ECG_limits ] = plot_ecg_heartbeat(ECG_struct.signal, lead_idx, this_all_anns, start_idx, anns_under_edition_idx(hb_idx) , hb_detail_window, ECG_struct.header, filtro, ECG_axes_hdl, ECG_limits);
            
            if( length(lead_idx) > 1 )
                aux_str = rowvec(colvec([repmat(',', length(lead_idx), 1) ECG_struct.header.desc(lead_idx,:) ]'));
                title(ECG_axes_hdl, ['Heartbeat ' num2str(hb_idx) ' : Leads ' aux_str(2:end) ] )
            else
                title(ECG_axes_hdl, ['Heartbeat ' num2str(hb_idx) ' : Lead ' ECG_struct.header.desc(lead_idx,:)] )
            end

    %         zoom reset
            cellfun(@(a)( set(a,'ButtonDownFcn',@inspect_ECG)), ECG_hdl);            
            set(ECG_axes_hdl,'ButtonDownFcn',@inspect_ECG);            

        elseif (strcmp(event_obj.Key,'rightarrow'))
            %% right arrow

            if(bAnnsEdited)
                update_annotations();
            end

            AnnNames_idx = AnnNames_idx + 1;
            if( AnnNames_idx > length(AnnNames) )
                AnnNames_idx = 1;
            end

            update_title_efimero(['Using ' AnnNames{AnnNames_idx,1} ' annotations (' num2str(ratios(AnnNames_idx)) ')'], 5 );
            
            undo_buffer_idx = 1;

            anns_under_edition = unique(round(colvec( ECG_struct.(AnnNames{AnnNames_idx,1}).(AnnNames{AnnNames_idx,2}) )));

            bAnnsEdited = false;

            selected_hb_idx = [];

            Redraw();


        elseif (strcmp(event_obj.Key,'leftarrow'))
            %% left arrow

            if(bAnnsEdited)
                update_annotations();
            end

            AnnNames_idx = AnnNames_idx - 1;
            if( AnnNames_idx < 1  )
                AnnNames_idx = length(AnnNames);
            end

            update_title_efimero(['Using ' AnnNames{AnnNames_idx,1} ' annotations (' num2str(ratios(AnnNames_idx)) ')'], 5 );
            
            undo_buffer_idx = 1;

            anns_under_edition = unique(round(colvec( ECG_struct.(AnnNames{AnnNames_idx,1}).(AnnNames{AnnNames_idx,2}) )));

            bAnnsEdited = false;

            selected_hb_idx = [];

            Redraw();

        else
%             bLockRRserie = false;
%             bLockScatter = false;
%             lead_idx = 1:ECG_struct.header.nsig;
        end

    end
    
    function PushUndoAction()

        if( undo_buffer_idx < length(undo_buffer) )
            %me cargo todo lo que est� por rehacerse
            undo_buffer = undo_buffer(1:undo_buffer_idx);
        end

        undo_buffer{undo_buffer_idx} = anns_under_edition;
        undo_buffer_idx = undo_buffer_idx + 1;
    
    end

    function this_point = get_current_point()
        
        point = get(gca,'CurrentPoint');
        
        lanns_under_ed = length(anns_under_edition_idx);
        [~, aux_val] = sort( abs(point(1) - anns_under_edition(anns_under_edition_idx)) );
        aux_val = aux_val(1);

%         fprintf(1, '%d %d\n', point(1), aux_val);
        
        this_point = max(1, min(lanns_under_ed, aux_val ));        
    end

    function WindowButtonDownCallback2D(obj,event_obj)
        
        if (strcmp(get(fig_hdl,'SelectionType'),'normal'))

            prev_u = get(fig_hdl, 'units');
            set(fig_hdl, 'units','normalized');
            crd = get(fig_hdl, 'CurrentPoint');
            xp = crd(1); 
            yp = crd(2);
            set(fig_hdl, 'units', prev_u);

            if( ~isempty([min_y_drag max_x_drag min_y_drag]) && xp <= max_x_drag && yp <= max_y_drag && yp >= min_y_drag)
                DragMouseBegin()
            end
        end
        
    end

    function WindowButtonUpCallback2D(obj,event_obj)

        DragMouseEnd()

        if ( fIsDragTimeAllowed )

            fIsDragTimeAllowed = false;
        
            ECG_struct.signal = ECG_w.read_signal(start_idx, end_idx + 10 * ECG_struct.header.freq );

            hb_idx = 1;
            
            Redraw();
            
            set(fig_hdl, 'WindowButtonMotionFcn', PrevStateWindowButtonMotionFcn);
            
        end
            
    end

    function inspect_RRserie(obj,event_obj)

        if (strcmp(get(fig_hdl,'SelectionType'),'extend'))
            % Show ECG

            aux_RR = RRserie{1};
            aux_RR = aux_RR(anns_under_edition_idx);
            
            point1 = get(gca,'CurrentPoint');    % button down detected
            rbbox;
            point2 = get(gca,'CurrentPoint');    % button up detected            

            if( ~isempty(selected_hb_idx) )
                delete(RRserie_selection_hdl)
                delete(RRscatter_selection_hdl)
            end

            [~, aux_val] = sort( abs(point1(1,1) - anns_under_edition(anns_under_edition_idx)) );
            xlims = aux_val(1);
            [~, aux_val] = sort( abs(point2(1,1) - anns_under_edition(anns_under_edition_idx)) );
            xlims = sort([xlims aux_val(1)]);

            ylims = sort([ point1(1,2) point2(1,2) ]);
            lanns_under_ed = length(anns_under_edition_idx);
            bAux = false(lanns_under_ed,1);
            bAux(max(1,xlims(1)):min(lanns_under_ed,xlims(2))) = true;
            selected_hb_idx = find( bAux & aux_RR > ylims(1) & aux_RR < ylims(2) );

            if( isempty(selected_hb_idx) )
                % probar en el rango de x solamente
                selected_hb_idx = find( bAux );
            end

            if( ~isempty(selected_hb_idx) )

                update_title_efimero([num2str(length(selected_hb_idx)) ' heartbeats selected.'], 5 );

                hold(Scatter_axes_hdl, 'on')
                aux_RRnext = [aux_RR(2:end); aux_RR(end)];
                RRscatter_selection_hdl = plot(Scatter_axes_hdl, aux_RR(selected_hb_idx), aux_RRnext(selected_hb_idx), 'xg' );
                hold(Scatter_axes_hdl, 'off')

                hold(RRserie_axes_hdl, 'on')
                RRserie_selection_hdl = plot(RRserie_axes_hdl, anns_under_edition(anns_under_edition_idx(selected_hb_idx)) , aux_RR(selected_hb_idx), 'og' );
                hold(RRserie_axes_hdl, 'off')

                min_hb_idx = min(selected_hb_idx);
                max_hb_idx = max(selected_hb_idx);
                cant_hb_idx = max_hb_idx - min_hb_idx + 1;
                
                if( (anns_under_edition(anns_under_edition_idx(max_hb_idx)) - anns_under_edition(anns_under_edition_idx(min_hb_idx))) <= (10*ECG_struct.header.freq) )
                    
                    if( bSeries )
                        this_all_anns = all_annotations_selected_serie_location;
                        aux_val = this_all_anns{1};
                        aux_val(serie_location_mask) = nan;
                        this_all_anns{1} = aux_val;
                    else
                        this_all_anns = all_annotations_selected;
                    end

                    if( ~bLockECG )
                        ECG_limits = [];
                    end            
                    
                    [ECG_hdl, ECG_limits ] = plot_ecg_heartbeat(ECG_struct.signal, lead_idx, this_all_anns, start_idx, anns_under_edition_idx(hb_idx) , hb_detail_window, ECG_struct.header, filtro, ECG_axes_hdl, ECG_limits);
                    cellfun(@(a)( set(a,'ButtonDownFcn',@inspect_ECG)), ECG_hdl);            
                    
                end

                if( length(lead_idx) > 1 )
                    aux_str = rowvec(colvec([repmat(',', length(lead_idx), 1) ECG_struct.header.desc(lead_idx,:) ]'));
                    title(ECG_axes_hdl, ['Heartbeat ' num2str(hb_idx) ' : Leads ' aux_str(2:end) ] )
                else
                    title(ECG_axes_hdl, ['Heartbeat ' num2str(min_hb_idx) ' : Lead ' ECG_struct.header.desc(lead_idx,:)] )
                end

    %             zoom reset

                set(ECG_axes_hdl,'ButtonDownFcn',@inspect_ECG);            

                % Zoom RR serie
                UpdateRRserieZoom();

            end

        elseif (strcmp(get(fig_hdl,'SelectionType'),'normal'))

            % Show ECG

            bChangeWin = false;
            
            if( ishandle(RRscatter_hb_idx_hdl) )
                delete(RRscatter_hb_idx_hdl)
            end
            if( ishandle(RRserie_hb_idx_hdl) )
                delete(RRserie_hb_idx_hdl)
            end

            point = get(gca,'CurrentPoint');

            lanns_under_ed = length(anns_under_edition);
            [~, aux_val] = sort( abs(point(1) - anns_under_edition) );
            aux_val = aux_val(1);

            x_units = point(1);
%             fprintf(1, 'xunits\n');

            hb_idx = get_current_point();
%             x_units = hb_idx;
            
            hold(Scatter_axes_hdl, 'on')
            aux_RR = RRserie{1};
            aux_RR = aux_RR(anns_under_edition_idx);
            RRnext = [aux_RR(2:end); aux_RR(end)];
%             RRscatter_selection_hdl = plot(Scatter_axes_hdl, RRserie(selected_hb_idx), RRnext(selected_hb_idx), 'xg' );
            RRscatter_hb_idx_hdl = plot(Scatter_axes_hdl, aux_RR(hb_idx), RRnext(hb_idx), 'or' );
            hold(Scatter_axes_hdl, 'off')

            hold(RRserie_axes_hdl, 'on')
            RRserie_hb_idx_hdl = plot(RRserie_axes_hdl, anns_under_edition(anns_under_edition_idx(hb_idx)), aux_RR(hb_idx), 'or' );
            hold(RRserie_axes_hdl, 'off')        

            %Zoom RR serie
            UpdateRRserieZoom();

            if( bSeries )
                this_all_anns = all_annotations_selected_serie_location;
                aux_val = this_all_anns{1};
                aux_val(serie_location_mask) = nan;
                this_all_anns{1} = aux_val;
            else
                this_all_anns = all_annotations_selected;
            end

            if( ~bLockECG )
                ECG_limits = [];
            end            
            
            if( isempty(anns_under_edition_idx) )
                [ECG_hdl, ECG_limits ] = plot_ecg_heartbeat(ECG_struct.signal, lead_idx, this_all_anns, start_idx, [] , hb_detail_window, ECG_struct.header, filtro, ECG_axes_hdl, ECG_limits);    
            else
                [ECG_hdl, ECG_limits ] = plot_ecg_heartbeat(ECG_struct.signal, lead_idx, this_all_anns, start_idx, anns_under_edition_idx(hb_idx) , hb_detail_window, ECG_struct.header, filtro, ECG_axes_hdl, ECG_limits);    
            end

            if( length(lead_idx) > 1 )
                aux_str = rowvec(colvec([repmat(',', length(lead_idx), 1) ECG_struct.header.desc(lead_idx,:) ]'));
                title(ECG_axes_hdl, ['Heartbeat ' num2str(hb_idx) ' : Leads ' aux_str(2:end) ] )
            else
                title(ECG_axes_hdl, ['Heartbeat ' num2str(hb_idx) ' : Lead ' ECG_struct.header.desc(lead_idx,:)] )
            end

    %         zoom reset

            cellfun(@(a)( set(a,'ButtonDownFcn',@inspect_ECG)), ECG_hdl);            

            set(ECG_axes_hdl,'ButtonDownFcn',@inspect_ECG);            

            
        elseif (strcmp(get(fig_hdl,'SelectionType'), 'alt'))
            
            point = get(gca,'CurrentPoint');
            
            x_units = point(1);
            
            bChangeWin = true;
            
            DragMouseBegin()
            

        end
        
    end

    function SearchPattern()

        point1 = get(ECG_axes_hdl,'CurrentPoint');    % button down detected
        rbbox;
        point2 = get(ECG_axes_hdl,'CurrentPoint');    % button up detected            

        xlims = round(sort([ point1(1,1) point2(1,1) ]));
        aux_seq = xlims(1):xlims(2);

        if( ishandle(Pattern_hdl) )
            delete(Pattern_hdl)
        end
        
        update_title_efimero('Filtering ECG ...', 5 );
        % restrict only to the first signal
        lead_idx = lead_idx(1);
        llead_idx = length(lead_idx);
        
        if( isempty(filtro) )
            ECG = ECG_struct.signal(:,lead_idx);
        else
            ECG = filter(filtro, flipud(ECG_struct.signal(:,lead_idx)) );
            ECG = filter(filtro, flipud(ECG) );
        end

        hold(ECG_axes_hdl, 'on')
        Pattern_hdl = plot(ECG_axes_hdl, aux_seq, ECG(aux_seq,:), '.g' );
        hold(ECG_axes_hdl, 'off')        
        drawnow;

        % decimation factor
        if( isempty(ECG_w) )
            ndown = 1;
            pattern2detect = ECG_struct.signal(aux_seq,lead_idx);
            pattern2detect = bsxfun( @minus, pattern2detect, mean(pattern2detect));
            
        else            
            % resampling of the signal
            Nqrs = round( .09 * 250); % samples of a normal QRS sampled at 250 Hz
            fs_target = round(Nqrs / length(aux_seq) * ECG_struct.header.freq);
            ndown = round( ECG_struct.header.freq / fs_target);
            aux_pre = round(0.2*ECG_struct.header.freq);
            pattern2detect_orig = ECG_struct.signal(aux_seq,lead_idx);
            aux_seq = (xlims(1) - aux_pre ):(xlims(2) + aux_pre);
            pattern2detect = ECG_struct.signal(aux_seq,lead_idx);
            pattern2detect = resample(pattern2detect, 1, ndown, 4);
            pattern2detect = pattern2detect(round(aux_pre/ndown):(end-round(aux_pre/ndown)),:);
            pattern2detect = bsxfun( @minus, pattern2detect, mean(pattern2detect));
            pattern2detect = pattern2detect .* repmat(hamming(length(pattern2detect)), 1, llead_idx);
            
        end
        

%         [nsamp_pattern nsig_pattern] = size(pattern2detect);
%         max_idx = max_index(pattern2detect);

    %     cross_corr = 1/sqrt(var(ECG)*var(pattern2detect) )*conv(ECG, flipud(pattern2detect) , 'valid');
    %     pattern2detect_var = var(pattern2detect);
        update_title_efimero('Looking for the pattern ...', 5 );        
    %     similarity = arrayfun( @(a)( 1/sqrt(var(ECG(a-nsamp_pattern+1:a, :))* pattern2detect_var)*(rowvec(ECG(a-nsamp_pattern+1:a, :)-mean(ECG(a-nsamp_pattern+1:a, :))) * colvec(pattern2detect) ) ), colvec(nsamp_pattern:ECG_struct.header.nsamp));    

    %     similarity = arrayfun( @(a)( 1/sqrt(var(ECG(a-nsamp_pattern+1:a))* pattern2detect_var)*(rowvec(ECG(a-nsamp_pattern+1:a)-mean(ECG(a-nsamp_pattern+1:a))) * colvec(pattern2detect) ) ), colvec(nsamp_pattern:ECG_struct.header.nsamp));    
    %     similarity = [repmat(similarity(1,:),nsamp_pattern-1,1); similarity];

        
        if( isempty(ECG_w) )
            % short or easy memory handlable signals
            similarity = cellfun( @(a,b)( diff(conv( a, b, 'same' )) ), mat2cell(ECG_struct.signal(:,lead_idx), ECG_struct.header.nsamp, ones(1,llead_idx)), mat2cell(pattern2detect, diff(xlims)+1, ones(1,llead_idx)), 'UniformOutput', false);
            similarity = cell2mat(cellfun( @(a,b)( diff(conv( a, b, 'same' )) ), similarity, mat2cell(flipud(pattern2detect), diff(xlims)+1, ones(1,llead_idx)), 'UniformOutput', false));
        %     similarity = [repmat(similarity(1,:),round((nsamp_pattern-1)/2),1); similarity];
            similarity = abs(mean(similarity,2));
        else

            if( ndown > 1 )
                aux_w = ECGwrapper('recording_name', ECG_w.recording_name);
                aux_w.output_path = tempdir;
                aux_w.ECGtaskHandle = 'arbitrary_function';
                aux_w.cacheResults = false;
                aux_w.ECGtaskHandle.lead_idx = lead_idx;
                aux_w.ECGtaskHandle.signal_payload = true;
                aux_w.ECGtaskHandle.function_pointer = @(a)(resample(a, 1, ndown, 4));
                aux_w.ECGtaskHandle.sampling_rate_out = ECG_struct.header.freq / ndown;
                aux_w.user_string = ['resample_to_' num2str(round(aux_w.ECGtaskHandle.sampling_rate_out)) '_for_leads_' num2str(sort(lead_idx)) ];
                aux_w.Run
                aux_ecg_filename = char(aux_w.Result_files);
            else
                aux_ecg_filename = ECG_w.recording_name;
            end
            
            % similarity calculation
            aux_w = ECGwrapper('recording_name', aux_ecg_filename);
            aux_w.output_path = tempdir;
            aux_w.ECGtaskHandle = 'arbitrary_function';
            aux_w.cacheResults = false;
            % always one leaded signal
            aux_w.ECGtaskHandle.lead_idx = 1;
            aux_w.ECGtaskHandle.signal_payload = true;
            
            aux_w.user_string = ['similarity_calc_for_lead_' num2str(sort(lead_idx)) ];
            aux_w.ECGtaskHandle.sampling_rate_out = ECG_struct.header.freq / ndown;
            aux_w.ECGtaskHandle.function_pointer = @similarity_calculation;
            aux_w.ECGtaskHandle.payload = pattern2detect;
            aux_w.Run
            
            aux_similarity_filename = char(aux_w.Result_files);
        end
        
        prev_fig = gcf();
        
        fig2_hdl = figure(2);
        clf();
        set(fig2_hdl, 'Position', [ maximized_size(3:4) maximized_size(3:4) ] .* [ 0.05 0.13 0.95 0.9] );
        
        win_sample = round(3*ECG_struct.header.freq / ndown);
        break_sample = round(1*ECG_struct.header.freq / ndown);
        n_excerpts = 5;
        aux_idx = 1:win_sample;
        aux_windows = randsample( 1:(ECG_struct.header.nsamp/ndown)-win_sample, n_excerpts);
        sig_breaks = nan(break_sample, llead_idx + 1 );
        
        if( isempty(ECG_w) )
            
            aux_val = sig_breaks;
            for aux_start = (aux_windows+1)
%                 aux_val1 = [bsxfun( @minus, ECG(aux_start:(aux_start+win_sample),:), mean(ECG(aux_idx))) similarity(aux_start:(aux_start+win_sample))-mean(similarity(aux_start:(aux_start+win_sample))) ];
                aux_val1 = [bsxfun( @minus, ECG(aux_start:(aux_start+win_sample),:), mean(ECG(aux_idx))) similarity(aux_start:(aux_start+win_sample)) ];
                aux_val = [aux_val; aux_val1; sig_breaks ];
            end
            
        else
            
            aux_w = ECGwrapper('recording_name', aux_similarity_filename);
            aux_w_ecg = ECGwrapper('recording_name', aux_ecg_filename);
            
            aux_val = sig_breaks;
            for aux_start = (aux_windows+1)
%                 aux_val1 = [aux_w_ecg.read_signal( aux_start, aux_start + win_sample ) aux_w.read_signal( aux_start, aux_start + win_sample ) ];
%                 aux_val1 = [aux_val1(:, lead_idx) aux_w.read_signal( aux_start, aux_start + win_sample ) ];
                
                aux_val1 = aux_w_ecg.read_signal( aux_start, aux_start + win_sample );
                
                aux_val1 = aux_val1(:,lead_idx);
                
                aux_val1 = [aux_val1 aux_w.read_signal( aux_start, aux_start + win_sample ) ];
                
                aux_val = [aux_val; [bsxfun(@minus, aux_val1(:,1:end-1), mean(aux_val1(:,1:end-1)) ) aux_val1(:,end)]; sig_breaks ];
            end
        end
        
        aux_thr_scale = max(abs(aux_val));
        aux_val = bsxfun(@times, aux_val, 1./aux_thr_scale);
        aux_val(:,1:llead_idx) = (aux_val(:,1:llead_idx)*0.5) - 0.5;
        aux_hdls = plot(aux_val);
        axes_hdl = gca();
        set(axes_hdl, 'Position', [ 0.015 0.025 0.97 0.92] );
        title('Select the detection threshold to use in the similarity function')
        
        set(axes_hdl, 'Ytick', []);
        
        aux_val = sort([ break_sample+(0:(win_sample+break_sample):(n_excerpts-1)*(win_sample+break_sample)) break_sample+win_sample+(0:(win_sample+break_sample):(n_excerpts-1)*(win_sample+break_sample)) length(aux_val) ]);
        set(axes_hdl, 'Xtick', aux_val );
        
        aux_val = sort([ aux_windows (aux_windows + win_sample) ECG_struct.header.nsamp]);
        set(axes_hdl, 'XtickLabel', Seconds2HMS( aux_val ./ ECG_struct.header.freq ));
        
        legend(aux_hdls, {'ECG'; 'Similarity'} );

        figure(fig2_hdl);

        update_title_efimero('Select the threshold to use.', 5 );        
        
        [~, thr] = ginput(1);

        thr = thr * aux_thr_scale(2);
        
        commandwindow;
        
        key = input( 'Enter the minimum time in seconds between heartbeats (Default 0.3 s):\n' , 's');

        if( isempty(key) )
            detection_threshold = 0.3; % seconds
        else
            detection_threshold = str2double(key);
        end

        dt_samples = round(detection_threshold*ECG_struct.header.freq / ndown);
        
        xlims = get(axes_hdl, 'xlim');
        ylims = get(axes_hdl, 'ylim');
        
        dt_yloc = ylims(1) + 0.05*diff(ylims);
        dt_xloc = xlims(1) + 0.1*diff(xlims);
        
        hold(axes_hdl, 'on');
        arrow( [dt_xloc; dt_yloc], [dt_xloc + dt_samples; dt_yloc], 2, 0.5, [0 0 0], axes_hdl );
        text( dt_xloc, dt_yloc+0.01*diff(ylims), ['Min QRS sep '  Seconds2HMS(detection_threshold,2) ]);
        hold(axes_hdl, 'off');
        
        % limits to perform pattern match. Default whole recording.
        QRSxlims = 1;
        bContinue = true;

        
        while(bContinue)

            if( isempty(ECG_w) )
                % short or easy memory handlable signals
                ECG_struct.pattern_match.time = modmax(similarity, QRSxlims, thr, 1, round(detection_threshold*ECG_struct.header.freq) );
            else
                aux_w.ECGtaskHandle = 'arbitrary_function';
                aux_w.cacheResults = false;
                aux_w.ECGtaskHandle.lead_idx = 1;
                % generate QRS detections
                aux_w.ECGtaskHandle.signal_payload = false;
                aux_w.user_string = ['modmax_calc_for_leads_' num2str(sort(lead_idx)) ];
                aux_w.ECGtaskHandle.function_pointer = @arb_modmax;
                
                aux_payload.thr = thr;
                aux_payload.detection_threshold = detection_threshold;
                aux_w.ECGtaskHandle.payload = aux_payload;
                aux_w.Run
                
                % asume that the whole series keep in mem.
                aux_val = load(aux_w.Result_files{1});  
                ECG_struct.pattern_match.time = round(aux_val.result * ndown);
                
            end

    %         aux_idx = 1;
    %         while( ~isempty(aux_idx) )
    %             aux_idx = find(diff(ECG_struct.pattern_match.time) < round(0.15 * ECG_struct.header.freq));
    %             aux_idx2 = setdiff(1:length(ECG_struct.pattern_match.time), [colvec(aux_idx); colvec(aux_idx+1)]);
    %             [~, merged_times] = max([similarity(ECG_struct.pattern_match.time(aux_idx)) similarity(ECG_struct.pattern_match.time(aux_idx+1))],[],2);
    %             ECG_struct.pattern_match.time = sort([colvec(ECG_struct.pattern_match.time(aux_idx2)); colvec(ECG_struct.pattern_match.time(aux_idx( find(merged_times == 1) ))); colvec(ECG_struct.pattern_match.time(aux_idx( find(merged_times == 2) )+1 )) ]);
    %         end

    
            if( strcmp(AnnNames(end,1), cellstr('pattern_match') ) )
                aux_all_anns = all_annotations;
                aux_all_anns{end} = ECG_struct.pattern_match.time;
            else
                AnnNames = [AnnNames; cellstr('pattern_match') cellstr('time')];
                aux_all_anns = [all_annotations; {ECG_struct.pattern_match.time}];
            end

            if( isempty(ECG_w) )
                % only for short signals
                [ ratios, estimated_labs ] = CalcRRserieRatio(aux_all_anns, ECG_struct.header);
            else
                % ignore ratios and q measurements in long recordings.
                ratios = zeros(size(AnnNames,1),1);
                estimated_labs = cell(size(AnnNames,1),1);
                
            end
            
            all_annotations = aux_all_anns;
            
            AnnNames_idx = size(AnnNames,1);
            anns_under_edition = unique(round(colvec( ECG_struct.pattern_match.time )));

            hb_idx = 1;
            selected_hb_idx = [];

            undo_buffer_idx = 1;

            aux_val = anns_under_edition;

            bAnnsEdited = false;

            if( isempty(aux_val) )
                anns_under_edition = [];
                RRserie = {[]};
                RRserie_filt = {[]};
                all_annotations_selected_serie_location = [];
                serie_location_mask = [];
            else
                anns_under_edition = unique(round(colvec( aux_val )));
                
                [RRserie, RRserie_filt ] = RR_calculation(anns_under_edition, ECG_struct.header.freq, RRserie_filt{1});               
                RRserie = {RRserie * 1/ECG_struct.header.freq};          
                RRserie_filt = {RRserie_filt};
                
            end
            all_annotations_selected = {anns_under_edition};
            RR_idx = { find( anns_under_edition >= start_idx &  anns_under_edition <= end_idx ) };

            bSeriesChange = true;
            
            figure(prev_fig);
            Redraw();
%             figure(fig2_hdl);

            ocurrences = length(ECG_struct.pattern_match.time);
            
            update_title_efimero(['Threshold: ' num2str(thr) ' - found ' num2str(ocurrences) ' heartbeats with quality ' num2str(ratios(end)) ], 5 );        
            
            disp_string_framed('*[1,0.5,0]', ['Threshold: ' num2str(thr) ' - found ' num2str(ocurrences) ' heartbeats with quality ' num2str(ratios(end)) ] );        
            
            commandwindow;
            
            key = input(['[rt] to refine the QRS detection threshold.\n' ...
                         '[rx] to refine the time window to perform QRS detection.\n' ...
                         '[rh] to refine the minimum time between heartbeats.\n' ...
                         'any key to continue.\n' ...
                         ], 's');

            if( strcmp(key, 'rt') )
                figure(fig2_hdl);
                update_title_efimero('Select the threshold to use.', 5 );        
                
                [~, thr] = ginput(1);
                thr = thr * aux_thr_scale(2);
                
            elseif( strcmp(key, 'rh') )

                key = input( ['Enter the minimum time in seconds between heartbeats (Current ' sprintf('%3.3f', detection_threshold) ' s):\n'] , 's');

                if( isempty(key) )
                    detection_threshold = 0.3; % seconds
                else
                    detection_threshold = str2double(key);
                end

            elseif( strcmp(key, 'rx') )
                fig_hdl = figure(1);
                update_title_efimero('Click and drag the time interval to perform QRS detection in the RR series interval ...', 5 );        
                
                set(fig_hdl, 'CurrentAxes', RRserie_axes_hdl);
                waitforbuttonpress;
                point1 = get(RRserie_axes_hdl,'CurrentPoint');    % button down detected
                rbbox;
                point2 = get(RRserie_axes_hdl,'CurrentPoint');    % button up detected            
                QRSxlims = round(sort([ point1(1,1) point2(1,1) ]));
                QRSxlims = [ QRSxlims(1) QRSxlims(1) + max(diff(QRSxlims), 10 * ECG_struct.header.freq )];
            else
                bContinue = false;
            end

        end

        figure(prev_fig);
        update_title_efimero('Search pattern finished.', 5 );        

        cant_anns = size(AnnNames,1);
        aux_str = repmat( ' - ',cant_anns,1);

        [~, best_detections_idx] = sort(ratios, 'descend');
        
        aux_val = 1:length(ratios);
        [~, annotations_ranking] = sort(aux_val(best_detections_idx));
        
        set(annotation_list_control, 'string', [ char(cellstr(num2str((1:cant_anns)'))) aux_str char(AnnNames(:,1)) repmat( ' (',cant_anns,1) num2str(round(colvec(ratios * 1000))) aux_str num2str(colvec(annotations_ranking))  repmat( ')',cant_anns,1)  ] );
        set(annotation_under_edition_label, 'string', [ 'Annotation under edition: ' char(AnnNames( AnnNames_idx ,1)) ' (' num2str(ratios(AnnNames_idx)) ')' ])
        set(annotation_list_control, 'Value', AnnNames_idx);    

        close(fig2_hdl);

    end
    
    function UpdateRRserieZoom()

        if( isempty(hb_idx) )
            return
        end

        bRRempty = all(cellfun( @(a)(isempty(a)), RRserie));
        
        if( bRRempty )
            
            cla(RRserie_zoom_axes_hdl, 'reset');
            
        else
            
            aux_RR = RRserie{1};
            aux_idx = find( anns_under_edition(anns_under_edition_idx) >= (anns_under_edition(anns_under_edition_idx(hb_idx)) - round(win_size_zoom/2*ECG_struct.header.freq)) & anns_under_edition(anns_under_edition_idx) <= (anns_under_edition(anns_under_edition_idx(hb_idx)) + round(win_size_zoom/2*ECG_struct.header.freq)) );
            if(length(aux_idx) < 3 )
                min_idx = find(~isnan(anns_under_edition(anns_under_edition_idx)),1);
                max_idx = find(~isnan(anns_under_edition(anns_under_edition_idx)),1, 'last');
                aux_idx = max(min_idx, hb_idx - 1 );
                aux_idx = aux_idx:min( max_idx, max(hb_idx, aux_idx) + 1 );
            else
                hb_detail_window = round( length(aux_idx) / 2 );
            end
            
            aux_idx2 = cellfun( @(this_anns)( find( this_anns >= anns_under_edition(anns_under_edition_idx(aux_idx(1))) &  this_anns <= anns_under_edition(anns_under_edition_idx(aux_idx(end))) ) ), all_annotations_selected, 'UniformOutput', false);

            aux_idx3 = find(~all(cellfun( @(a)(isempty(a)), aux_idx2)));
            
            cla(RRserie_zoom_axes_hdl, 'reset');
                
            if( ~isempty(aux_idx3) )
            
        %         RRserie2 = cellfun( @(this_anns)( colvec(diff(this_anns)) ), all_annotations_selected, 'UniformOutput', false);
        %         RRserie2 = cellfun( @(this_rr_serie)( [this_rr_serie(1); this_rr_serie] ), RRserie, 'UniformOutput', false);

                hold(RRserie_zoom_axes_hdl, 'on')
                RRserie_zoom_hdl = cellfun( @(this_anns, this_rr_serie, this_idx, ii)( plot(RRserie_zoom_axes_hdl, this_anns(this_idx) , this_rr_serie(this_idx), 'LineStyle', ':', 'Marker', all_markers{ii}, 'MarkerEdgeColor', ColorOrder(ii,:), 'Color', ColorOrder(ii,:) )  ), all_annotations_selected(aux_idx3), RRserie(aux_idx3), aux_idx2(aux_idx3), num2cell(colvec(aux_idx3)), 'UniformOutput', false );

                this_ylims = get(RRserie_axes_hdl, 'Ylim' );
                this_xlims_orig = get(RRserie_zoom_axes_hdl, 'Xlim' );

                set(RRserie_zoom_axes_hdl, 'Ylim', this_ylims );

                % blue box around
                this_xlims = this_xlims_orig + 0.01*[ diff(this_xlims_orig) -diff(this_xlims_orig) ];
                this_ylims = this_ylims + 0.015*[ diff(this_ylims) -diff(this_ylims) ];

                set(fig_hdl,'CurrentAxes', RRserie_zoom_axes_hdl);
                aux_hdl = patch([this_xlims(1) this_xlims(1) this_xlims(2) this_xlims(2) this_xlims(1) ], [this_ylims(1) this_ylims(2) this_ylims(2) this_ylims(1) this_ylims(1)], [1 1 1], 'EdgeColor', [0 0 1], 'LineWidth', 1.5, 'ButtonDownFcn', @inspect_RRserie );
                set(aux_hdl, 'FaceColor', 'none')
                uistack(aux_hdl, 'bottom');

                if( isempty(RRserie_zoom_zoombars_hdl) || ~ishandle(RRserie_zoom_zoombars_hdl) )
                    set(fig_hdl, 'CurrentAxes', RRserie_axes_hdl);
                    RRserie_zoom_zoombars_hdl = patch([this_xlims(1) this_xlims(1) this_xlims(2) this_xlims(2) this_xlims(1) ], [this_ylims(1) this_ylims(2) this_ylims(2) this_ylims(1) this_ylims(1)], [190 238 238]/255, 'EdgeColor', [0 0 1], 'LineWidth', 0.5, 'ButtonDownFcn', @inspect_RRserie);
                    uistack(RRserie_zoom_zoombars_hdl, 'bottom');
                else
                    set(RRserie_zoom_zoombars_hdl, 'Xdata', [this_xlims(1) this_xlims(1) this_xlims(2) this_xlims(2) this_xlims(1) ]);
                end

                if( ~isempty(aux_RR) )
                    [~, aux_idx2] = intersect(selected_hb_idx, aux_idx);
                    RRserie_zoom_hdl = [RRserie_zoom_hdl; colvec(arrayfun(@(a,b)( plot(RRserie_zoom_axes_hdl, a, b, 'og')), anns_under_edition(anns_under_edition_idx(selected_hb_idx(aux_idx2))), aux_RR(anns_under_edition_idx(selected_hb_idx(aux_idx2))), 'UniformOutput', false ) ) ];
                end

                if( ~isempty(aux_RR) && hb_idx < length(aux_RR) )
                    RRserie_zoom_hdl = [RRserie_zoom_hdl; colvec(arrayfun(@(a,b)( plot(RRserie_zoom_axes_hdl, a, b, 'or')), anns_under_edition(anns_under_edition_idx(hb_idx)), aux_RR(anns_under_edition_idx(hb_idx)), 'UniformOutput', false ) )];
                end

                set(RRserie_zoom_axes_hdl, 'Xlim', this_xlims_orig );

                hold(RRserie_zoom_axes_hdl, 'off');

                xlabel(RRserie_zoom_axes_hdl, 'Time');
                ylabel(RRserie_zoom_axes_hdl, 'Serie value');

                set(RRserie_zoom_axes_hdl,'ButtonDownFcn',@inspect_RRserie);  
                cellfun( @(a)(set(a, 'ButtonDownFcn', @inspect_RRserie ) ), RRserie_zoom_hdl );  

                aux_hb_idx = find(hb_idx == aux_idx);

                if( length(aux_idx) > 5 )
                    cant_ticks = 5;
                    if( isempty(aux_hb_idx) )
                        aux_idx = round(linspace(aux_idx(1), aux_idx(end), cant_ticks));
                    else
                        cant_ticks = cant_ticks - 1;
                        aux_idx = sort(unique([ round(linspace(aux_idx(1), aux_idx(end), cant_ticks)) aux_idx(aux_hb_idx) ]));
                        aux_hb_idx = find(hb_idx == aux_idx);
                    end
                end

                set(RRserie_zoom_axes_hdl, 'XTick', rowvec(anns_under_edition(anns_under_edition_idx(aux_idx))) );
                aux_str = cellstr(num2str(colvec(anns_under_edition_idx(aux_idx))));
                if( ~isempty(anns_under_edition) && ~isnan(anns_under_edition(anns_under_edition_idx(hb_idx))) && hb_idx <= length(anns_under_edition_idx) )
                    aux_str{aux_hb_idx} = [aux_str{aux_hb_idx} ' (' Seconds2HMS((anns_under_edition(anns_under_edition_idx(hb_idx)) + base_start_time - 1)*1/ECG_struct.header.freq) ' - ' Seconds2HMS((anns_under_edition(anns_under_edition_idx(hb_idx)))*1/ECG_struct.header.freq) ')'];
                    set(RRserie_zoom_axes_hdl, 'XTickLabel', char(aux_str) );
                end
            end
        end
        
    end

    function DeleteAnnotations(obj,event_obj) 

        cant_anns = size(AnnNames,1);

        if( cant_anns == 0 )
            return;
        end

        if( strcmpi(questdlg('Are you sure ?', 'Delete annotations', 'No'), 'yes') )

            ann_idx = get(annotation_list_control, 'Value');
            ECG_struct = rmfield(ECG_struct, AnnNames{ann_idx,1});

            if( cant_anns == 1)
                % last annotation deleted
                ECG_struct.Default.time = [];
                AnnNames = {'Default' 'time'};
                set(annotation_list_control, 'string', '1 - Default' );
                set(annotation_under_edition_label, 'string', 'Annotation under edition: Default' )
                AnnNames_idx = 1;
                anns_under_edition = [];

            elseif( cant_anns > 1)
                aux_idx = 1:cant_anns;
                aux_idx(ann_idx) = [];
                AnnNames = AnnNames(aux_idx,:);
                ratios = ratios(aux_idx);
                annotations_ranking = annotations_ranking(aux_idx);

                cant_anns = size(AnnNames,1);
                aux_str = repmat( ' - ',cant_anns,1);

                set(annotation_list_control, 'string', [ char(cellstr(num2str((1:cant_anns)'))) aux_str char(AnnNames(:,1)) repmat( ' (',cant_anns,1) num2str(round(colvec(ratios * 1000))) aux_str num2str(colvec(annotations_ranking))  repmat( ')',cant_anns,1)  ] );

                AnnNames_idx = 1;
                set(annotation_under_edition_label, 'string', [ 'Annotation under edition: ' char(AnnNames( AnnNames_idx ,1)) ' (' num2str(ratios(AnnNames_idx)) ')' ])
                set(annotation_list_control, 'Value', AnnNames_idx);
                anns_under_edition = unique(round(colvec( ECG_struct.(AnnNames{AnnNames_idx,1}).(AnnNames{AnnNames_idx,2}) )));

            end

            if( isfield(ECG_struct, 'series_quality' ) ) 
                ECG_struct.series_quality.AnnNames = AnnNames;
                ECG_struct.series_quality.ratios = ratios;
            end
            
            undo_buffer_idx = 1;

            bRecEdited = true;
            bAnnsEdited = false;

            selected_hb_idx = [];        

            Redraw();

        end
        
    end

    function MergeAnnotations(obj,event_obj) 


        ann_idx = get(annotation_list_control, 'Value');
        
        if( length(ann_idx) > 1 )
        
            merged_detections = unique(round(colvec( ECG_struct.(AnnNames{ann_idx(1),1}).(AnnNames{ann_idx(1),2}) )));

            for ii = 2:length(ann_idx)

                detB = unique(round(colvec( ECG_struct.(AnnNames{ann_idx(ii),1}).(AnnNames{ann_idx(ii),2}) )));

                merged_detections = merge_QRS_detections(merged_detections, detB, ECG_struct.header);

            end
        
            aux_val = sum(~cellfun( @isempty, strfind( AnnNames(:,1), 'merged' )));
            
            aux_str = [ 'merged_' num2str(aux_val+1)];
            
            ECG_struct.(aux_str).time = merged_detections;
            
            AnnNames = [AnnNames; cellstr(aux_str) cellstr('time')];
            aux_all_anns = [all_annotations; {merged_detections}];

            if( isempty(ECG_w) )
                % only for short signals
                [ ratios, estimated_labs ] = CalcRRserieRatio(aux_all_anns, ECG_struct.header);
            else
                % ignore ratios and q measurements in long recordings.
                ratios = zeros(size(AnnNames,1),1);
                estimated_labs = cell(size(AnnNames,1),1);
            end

            [ratios, best_detections_idx] = sort(ratios, 'descend');

            aux_val = 1:length(ratios);

            [~, annotations_ranking] = sort(aux_val(best_detections_idx));
            
            all_annotations = aux_all_anns;

            AnnNames_idx = size(AnnNames,1);
            anns_under_edition = merged_detections;

            cant_anns = size(AnnNames,1);
            aux_str = repmat( ' - ',cant_anns,1);

            set(annotation_list_control, 'string', [ char(cellstr(num2str((1:cant_anns)'))) aux_str char(AnnNames(:,1)) repmat( ' (',cant_anns,1) num2str(round(colvec(ratios * 1000))) aux_str num2str(colvec(annotations_ranking))  repmat( ')',cant_anns,1)  ] );
            set(annotation_under_edition_label, 'string', [ 'Annotation under edition: ' char(AnnNames( AnnNames_idx ,1)) ' (' num2str(ratios(AnnNames_idx)) ')' ])

            set(annotation_list_control, 'Value', AnnNames_idx);
            
            hb_idx = 1;
            selected_hb_idx = [];

            undo_buffer_idx = 1;

            aux_val = anns_under_edition;

            bAnnsEdited = false;

            if( isempty(aux_val) )
                anns_under_edition = [];
                RRserie = {[]};
                RRserie_filt = {[]};
                all_annotations_selected_serie_location = [];
                serie_location_mask = [];
            else
                [RRserie, RRserie_filt ] = RR_calculation(anns_under_edition, ECG_struct.header.freq, RRserie_filt{1});               
                RRserie = {RRserie * 1/ECG_struct.header.freq};          
                RRserie_filt = {RRserie_filt};
            end
            all_annotations_selected = {anns_under_edition};
            RR_idx = { find( anns_under_edition >= start_idx &  anns_under_edition <= end_idx ) };

            bSeriesChange = true;        

            undo_buffer_idx = 1;

            bRecEdited = true;
            bAnnsEdited = false;

            selected_hb_idx = [];   

            Redraw();        

        end
            
    end

    function ChangeRecordingSelected(obj,event_obj) 


        if (strcmp(get(fig_hdl,'SelectionType'),'open'))
            %Double click        

        else
            %Single click

            rec_selected = get(obj, 'Value');

            if( rec_selected ~= rec_idx )
                % change of annotation

                if(bRecEdited)

                    update_annotations();

                    if( bLoadECG )
                        update_title_efimero('Saving data ...', 5 );
                        save(rec_path, '-struct', 'ECG_struct');        
                        update_title_efimero(['Saved ' rec_path], 5 );
                    else

                        if( isfield(ECG_struct, 'series_quality' ) ) 
                            anns_struct.series_quality = ECG_struct.series_quality;
                        end
                        for ii = 1:size(AnnNames,1)
                            anns_struct.(AnnNames{ii,1}) = ECG_struct.(AnnNames{ii,1});
                        end
                        assignin( 'caller', OutputVarName, anns_struct );
                        update_title_efimero(sprintf('Saving ''%s'' variable in caller workspace', OutputVarName), 5 );

                        bAnnsEdited = false;
                        bRecEdited = false;

                    end
                end

                rec_idx = get(recordings_control, 'Value' );

                DoRecording();


            end

        end

    end
    
    function ChangeLeadsSelected(obj,event_obj) 

        if (strcmp(get(fig_hdl,'SelectionType'),'open'))
            %Double click        

        else
            %Single click

            leads_selected = get(obj, 'Value');

            if(length(leads_selected ) == 1 && ( length(lead_idx) ~= length(leads_selected) || leads_selected ~= lead_idx) )
                % change of annotation

                standarize_ECG_view = false;

                lead_idx = leads_selected;

                cla(ECG_axes_hdl);

                if( bSeries )
                    this_all_anns = all_annotations_selected_serie_location;
                    aux_val = this_all_anns{1};
                    aux_val(serie_location_mask) = nan;
                    this_all_anns{1} = aux_val;
                else
                    this_all_anns = all_annotations_selected;
                end
                
                if( ~bLockECG )
                    ECG_limits = [];
                end            
                
                
                if( isempty(anns_under_edition_idx) )
                    [ECG_hdl, ECG_limits ] = plot_ecg_heartbeat(ECG_struct.signal, lead_idx, this_all_anns, start_idx, [] , hb_detail_window, ECG_struct.header, filtro, ECG_axes_hdl, ECG_limits);    
                else
                    [ECG_hdl, ECG_limits ] = plot_ecg_heartbeat(ECG_struct.signal, lead_idx, this_all_anns, start_idx, anns_under_edition_idx(hb_idx) , hb_detail_window , ECG_struct.header, filtro, ECG_axes_hdl, ECG_limits);    
                end

                title(ECG_axes_hdl, ['Heartbeat ' num2str(hb_idx) ' : Lead ' ECG_struct.header.desc(lead_idx,:)] )

                cellfun(@(a)( set(a,'ButtonDownFcn',@inspect_ECG)), ECG_hdl);            
                set(ECG_axes_hdl,'ButtonDownFcn',@inspect_ECG);   


            elseif(length(leads_selected ) > 1)

                lead_idx = leads_selected;

                cla(ECG_axes_hdl);

                if( bSeries )
                    this_all_anns = all_annotations_selected_serie_location;
                    aux_val = this_all_anns{1};
                    aux_val(serie_location_mask) = nan;
                    this_all_anns{1} = aux_val;
                else
                    this_all_anns = all_annotations_selected;
                end
                
                if( ~bLockECG )
                    ECG_limits = [];
                end            
                
                if( isempty(anns_under_edition_idx) )
                    [ECG_hdl, ECG_limits ] = plot_ecg_heartbeat(ECG_struct.signal, lead_idx, this_all_anns, start_idx, [] , hb_detail_window, ECG_struct.header, filtro, ECG_axes_hdl, ECG_limits);    
                else
                    [ECG_hdl, ECG_limits ] = plot_ecg_heartbeat(ECG_struct.signal, lead_idx, this_all_anns, start_idx, anns_under_edition_idx(hb_idx) , hb_detail_window , ECG_struct.header, filtro, ECG_axes_hdl, ECG_limits);    
                end
                

                aux_str = rowvec(colvec([repmat(',', length(lead_idx), 1) ECG_struct.header.desc(lead_idx,:) ]'));
                title(ECG_axes_hdl, ['Heartbeat ' num2str(hb_idx) ' : Leads ' aux_str(2:end) ] )

                cellfun(@(a)( set(a,'ButtonDownFcn',@inspect_ECG)), ECG_hdl);            
                set(ECG_axes_hdl,'ButtonDownFcn',@inspect_ECG);   


            end

        end

    end
    
    function update_annotations()
        
        corrected_prefix = 'corrected_';
        
        if( bSeries )
            aux_val = all_annotations_selected_serie_location{1};
            aux_val(serie_location_mask) = nan;
            aux_val = ( aux_val - anns_under_edition) / ECG_struct.header.freq;
            aux_val = [anns_under_edition aux_val];
        else
            aux_val = unique(round(colvec( anns_under_edition )));
        end

        if( isempty(strfind(AnnNames{AnnNames_idx,1}, corrected_prefix )) )
            % modifying an original annotation, duplicate annotation
            
            ii = 1;
            while( isfield(ECG_struct, [ corrected_prefix AnnNames{AnnNames_idx,1}]) )
                corrected_prefix = [corrected_prefix 'v' num2str(ii) '_' ];
                ii = ii+1;
            end
            
            ECG_struct.([ corrected_prefix AnnNames{AnnNames_idx,1}]).(AnnNames{AnnNames_idx,2}) = aux_val;
%             ECG_struct.([ corrected_prefix AnnNames{AnnNames_idx,1}]) = ECG_struct.(AnnNames{AnnNames_idx,1});
%             ECG_struct = rmfield(ECG_struct, AnnNames{AnnNames_idx,1});


            % add it to the ann list
%             AnnNames( AnnNames_idx , : ) = { [ corrected_prefix AnnNames{AnnNames_idx,1}] AnnNames{AnnNames_idx,2} };
            AnnNames = [{ [ corrected_prefix AnnNames{AnnNames_idx,1}] AnnNames{AnnNames_idx,2} }; AnnNames];
            ratios = [ ratios(AnnNames_idx); ratios ];
            AnnNames_idx = 1;
            annotations_ranking = [ 1; colvec(annotations_ranking+1) ];

            cant_anns = size(AnnNames,1);
            aux_str = repmat( ' - ',cant_anns,1);
            set(annotation_list_control, 'string', [ char(cellstr(num2str((1:cant_anns)'))) aux_str char(AnnNames(:,1)) repmat( ' (',cant_anns,1) num2str(round(colvec(ratios * 1000))) aux_str num2str(colvec(annotations_ranking))  repmat( ')',cant_anns,1)  ] );

            if( isfield(ECG_struct, 'series_quality' ) ) 
                ECG_struct.series_quality.AnnNames = AnnNames;
                ECG_struct.series_quality.ratios = ratios;
            end
            
        else
            % updating an already corrected ann.
            ECG_struct.(AnnNames{AnnNames_idx,1}).(AnnNames{AnnNames_idx,2}) = aux_val;
        end            

    end

    function bAux = IsClicked(this_hdl, this_xy )
       
        dfltUnits = get(this_hdl, 'Units');
        
        set(this_hdl, 'Units', 'pixels');
        
        aux_pos = get( this_hdl, 'Position');
        
        bAux = this_xy(1) >= aux_pos(1) & this_xy(1) <= (aux_pos(1) + aux_pos(3) ) & this_xy(2) >= aux_pos(2) & this_xy(2) <= (aux_pos(2) + aux_pos(4) );

        set(this_hdl, 'Units', dfltUnits);
        
    end

    function TimeOffsetClick(obj,event_obj) 

        if ( ~fIsDragTimeAllowed )

            [drag_timeScroll_start_x, drag_timeScroll_start_y ] = GetCursorCoordOnWindow();
           
            if( IsClicked(RRserie_global_axes_hdl, [drag_timeScroll_start_x, drag_timeScroll_start_y ]) )
                point = get(RRserie_global_axes_hdl,'CurrentPoint');
            elseif( IsClicked(RRserie_axes_hdl, [drag_timeScroll_start_x, drag_timeScroll_start_y ]) )
                point = get(RRserie_axes_hdl,'CurrentPoint');
            else
                return
            end
            
            x_timeScroll_units = point(1);
            
            if (strcmp(get(fig_hdl,'SelectionType'),'alt'))
                bChangeWin = true;
            else
                bChangeWin = false;
            end
            
            UpdateStartX( 0 );
            
            fIsDragTimeAllowed = true;
            
            PrevStateWindowButtonMotionFcn = get(fig_hdl, 'WindowButtonMotionFcn');
            set(fig_hdl, 'WindowButtonMotionFcn', @WindowButtonMotionCallback2D);
            
        end
        
    end


    function ChangeAnnotationsSelected(obj,event_obj) 


        if (strcmp(get(fig_hdl,'SelectionType'),'open'))
            %Double click        

            answer = char(inputdlg([ 'Enter the new name of the annotation ' char(AnnNames( AnnNames_idx ,1)) ], 'Change annotation name', 1, AnnNames( AnnNames_idx ,1)) );

            if( ~isempty(answer) && ischar(answer) )
                ann_idx = get(annotation_list_control, 'Value');
                ECG_struct.(answer) = ECG_struct.(AnnNames{ann_idx,1});
                ECG_struct = rmfield(ECG_struct, AnnNames{ann_idx,1});

                AnnNames{ ann_idx ,1} = answer;

                if( isfield(ECG_struct, 'series_quality' ) ) 
                    ECG_struct.series_quality.AnnNames = AnnNames;
                end
                cant_anns = size(AnnNames,1);
                aux_str = repmat( ' - ',cant_anns,1);

                if( isempty(ratios) )
                    set(annotation_list_control, 'string', [ char(cellstr(num2str((1:cant_anns)'))) aux_str char(AnnNames(:,1))   ] );
                    set(annotation_under_edition_label, 'string', [ 'Annotation under edition: ' char(AnnNames( AnnNames_idx ,1)) ])
                else
                    set(annotation_list_control, 'string', [ char(cellstr(num2str((1:cant_anns)'))) aux_str char(AnnNames(:,1)) repmat( ' (',cant_anns,1) num2str(round(colvec(ratios * 1000))) aux_str num2str(colvec(annotations_ranking))  repmat( ')',cant_anns,1)  ] );
                    set(annotation_under_edition_label, 'string', [ 'Annotation under edition: ' char(AnnNames( AnnNames_idx ,1)) ' (' num2str(ratios(AnnNames_idx)) ')' ])
                end
                
                bAnnsEdited = true;
                bRecEdited = true;   

            end

        else
            %Single click

            anns_selected = get(obj, 'Value');

            if(length(anns_selected ) == 1 && anns_selected ~= AnnNames_idx )
                % change of annotation

                if(bAnnsEdited)
                    update_annotations();
                end

                AnnNames_idx = anns_selected;

    %             disp( ['Using ' AnnNames{AnnNames_idx,1} ' annotations (' num2str(ratios(AnnNames_idx)) ')'] );

                set(annotation_under_edition_label, 'string', [ 'Annotation under edition: ' char(AnnNames( AnnNames_idx ,1)) ' (' num2str(ratios(AnnNames_idx)) ')' ])

                undo_buffer_idx = 1;

                aux_val = ECG_struct.(AnnNames{AnnNames_idx,1}).(AnnNames{AnnNames_idx,2});

                bAnnsEdited = false;

                selected_hb_idx = [];

                if( isempty(aux_val) )

                    anns_under_edition = [];
                    RRserie = {[]};
                    RRserie_filt = {[]};
                    all_annotations_selected_serie_location = [];
                    serie_location_mask = [];
                    
                else
                
                    if( bSeries )
                        [anns_under_edition, aux_idx ]= unique(round(colvec( aux_val(:,1) )));
                        RRserie = { aux_val(aux_idx,2) };

                        % absolute position
                        all_annotations_selected_serie_location = {anns_under_edition + round( aux_val(aux_idx,2) * ECG_struct.header.freq) };
                        serie_location_mask = false(size(anns_under_edition));
                    else
                        anns_under_edition = unique(round(colvec( aux_val )));
                        [RRserie, RRserie_filt ] = RR_calculation(anns_under_edition, ECG_struct.header.freq, RRserie_filt{1});               
                        RRserie = {RRserie * 1/ECG_struct.header.freq};          
                        RRserie_filt = {RRserie_filt};
                        
                    end
                
                end
                
                all_annotations_selected = {anns_under_edition};
                RR_idx = { find( anns_under_edition >= start_idx &  anns_under_edition <= end_idx ) };
                
                bSeriesChange = true;
                
                Redraw();


            elseif(length(anns_selected ) > 1)

                aux_anns = {anns_under_edition};
                if( bSeries )
                    aux_val = ECG_struct.(AnnNames{AnnNames_idx,1}).(AnnNames{AnnNames_idx,2});
                    aux_anns2 = { aux_val(:,2) };
                end
                anns_selected( anns_selected == AnnNames_idx) = [];

                for ii = rowvec(anns_selected)
                    aux_val = ECG_struct.(AnnNames{ii,1}).(AnnNames{ii,2});
                    aux_anns = [aux_anns; ... 
                                {aux_val(:,1)}
                                ];
                    if( bSeries )
                        aux_anns2 = [aux_anns2; ... 
                                    {aux_val(:,2)}
                                    ];
                        
                    end
                            
                end

                if( bSeries )
                    all_annotations_selected_serie_location = cellfun( @(a,b)(a + round( b * ECG_struct.header.freq) ), aux_anns, aux_anns2, 'UniformOutput', false);
                end
                
                all_annotations_selected = aux_anns;

                RR_idx = cellfun( @(this_anns)( find( this_anns >= start_idx &  this_anns <= end_idx ) ), all_annotations_selected, 'UniformOutput', false);
                
                if( bSeries )
                    RRserie = aux_anns2;
                else
                    [RRserie, RRserie_filt ]= cellfun( @(this_anns)( RR_calculation(this_anns, ECG_struct.header.freq) * 1/ECG_struct.header.freq ), all_annotations_selected, 'UniformOutput', false);
                    RRserie_filt = RRserie_filt * ECG_struct.header.freq;
                end
                
                bSeriesChange = true;
                
                Redraw();

            end

        end
    end

    function update_title_efimero( strTitle, delay )
       
        if( ishandle(title_efimero_hdl) )
            set(title_efimero_hdl, 'String',strTitle )
            set(title_efimero_hdl, 'Visible', 'on');
        else

            title_efimero_hdl = annotation( 'textbox', [0 0.98 0.3 0.01 ], ...
                                  'String', strTitle, ... 
                                  'Tag', 'title_efimero', ...
                                  'FontSize', 8, ...
                                  'Interpreter', 'none', ...
                                  'FitBoxToText', 'on', ...
                                  'HorizontalAlignment', 'left', ...
                                  'EdgeColor', 'none', ...
                                  'BackgroundColor', 'r' );

        end
        
        if( ~isinf(delay) && strcmpi(my_timer.Running, 'off') )
            my_timer.StartDelay = delay;
            start(my_timer)
        end
        
    end

    function timer_fcn(obj,event_obj)
        
        set(title_efimero_hdl, 'Visible', 'off');

%         if(~bPreserveFix)
%             % allow edition of the closer wave
%             bFixedWave = false;
%         end
        
    end

end    
    
        