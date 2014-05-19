function ann_output = QRScorrector(varargin)

%%
% Description: 
% This function implements a graphical user interface (GUI) to correct annotations
% possibly previous performed by an automatic algorithm. The idea is to
% easily visualize, compare and correct annotations. For this purpose the
% GUI presents several representations of the time evolution of the events
% (possibly but not necessarily heartbeats).
% 
% Arguments:
%     
%     +ECG: [numeric or cell] REQUIRED
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
% 
% Limits and Known bugs:
%   Probably a lot :( ... but dont panic! send me feedback if you need help.
% 
% Example:
% 
% % This example ...
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate  : 16/2/2012
% Last update: 7/2/2014

%% Constants
cached_filename = 'tmp_QRScorrector_cache.mat';

user_data.inv_dower = [   0.156 -0.010  -0.172 -0.074  0.122  0.231 0.239 0.194 ;
                -0.227  0.887 0.057 -0.019 -0.106 -0.022 0.041 0.048;
                0.022  0.102  -0.229 -0.310 -0.246 -0.063 0.055 0.108]';

cHeaderFieldNamesRequired = {'freq' 'nsamp' 'nsig' 'gain' 'adczero' };
cAnnotationsFieldNamesRequired = {'time' };
   
ann_output = [];

%% Argument parsing

%argument definition
p = inputParser;   % Create instance of inputParser class.
p.addParamValue('FilterSignals', true, @(x)( islogical(x)));
p.addParamValue('BuildArtificial', true, @(x)( islogical(x)));
p.addParamValue('recording', [], @(x)( ischar(x)));
p.addParamValue('recording_path', [], @(x)( ischar(x)));
p.addParamValue('recording_indexes', [], @(x)( isnumeric(x) && all(x > 0) ) );
p.addParamValue('ECG', [], @(x)(isnumeric(x)) );
p.addParamValue('ECG_header', [], @(x)(isstruct(x)) );
p.addParamValue('QRS_annotations', [], @(x)( (isnumeric(x) && all(x > 0)) || isstruct(x) ) );
p.addParamValue('tmp_path', [], @(x)(ischar(x)) );

try
    p.parse( varargin{:} );
catch MyError
    rethrow(MyError);
end


bFilter = p.Results.FilterSignals;
bBuildArtificial = p.Results.BuildArtificial;
recording = p.Results.recording;
recording_path = p.Results.recording_path;
recording_indexes = p.Results.recording_indexes;
ECG = p.Results.ECG;
ECG_header = p.Results.ECG_header;
ECG_annotations = p.Results.QRS_annotations;
tmp_path = p.Results.tmp_path;

% Dont know why this variable uses a lot of bytes to store at disk.
clear p

%% Argument parsing

user_data.bLoadECG = true;

if( ~isempty(ECG) )
    %% ECG already read

    ECG_struct.signal = ECG;
    
    if( isempty(ECG_header))
       error( 'QRScorrector:ArgCheck:InvalidHeader', 'Please provide the ECG header.\n\n' );
    else
        if( ~isfield(ECG_header, cHeaderFieldNamesRequired ) )
            strAux = [ repmat(' + ', length(cHeaderFieldNamesRequired), 1) char(cHeaderFieldNamesRequired) repmat('\n', length(cHeaderFieldNamesRequired), 1 ) ];
            error( 'QRScorrector:ArgCheck:InvalidHeader', ['Please provide the following fields in the header struct:\n ' rowvec(strAux') ] );
        end
    end

    ECG_struct.header = ECG_header;
    
    if( isempty(ECG_annotations))
        disp( 'No annotations provided.' )
    else
        if(isstruct(ECG_annotations))        
            if( ~isfield(ECG_annotations, cAnnotationsFieldNamesRequired ) )
                strAux = [ repmat(' + ', length(cAnnotationsFieldNamesRequired), 1) char(cAnnotationsFieldNamesRequired) repmat('\n', length(cAnnotationsFieldNamesRequired), 1 ) ];
                error( 'QRScorrector:ArgCheck:InvalidAnnotations', ['Please provide the following fields in the annotations struct:\n ' rowvec(strAux') ] );
            end
            ECG_struct.provided_anns = ECG_annotations;
        else
            ECG_struct.provided_anns.time = ECG_annotations;
        end
    end
    
    clear ECG ECG_header ECG_annotations
    
    recording_indexes = 1;
    rec_names.name = 'ECG_already_loaded.ecg';
    
    if( isempty(recording_path) )
        recording_path = [fileparts(mfilename('fullpath')) filesep 'tmp' filesep ];
    end   

    user_data.bLoadECG = false;
    
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
            fprintf(1, 'Loading cached ratios from %s.\n', filename)
            aux_val = load(filename);
            ratios = aux_val.ratios;
        end
            
        if( ~bAux || length(ratios) ~= lrec_names)
            
            ratios = nan(lrec_names,1);

            pb = progress_bar( 'Processing ratios', sprintf( '%d recordings found', lrec_names ) );

            pb.Loop2do = lrec_names;

            for ii = recording_indexes
                pb.start_loop();

                pb.checkpoint(sprintf( 'Loading recording %d', ii ));
                pepe = load([recording_path rec_names(ii).name ], 'series_quality');
                ratios(ii) = max(pepe.series_quality.ratios);

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

user_data.side_plot = [];

user_data.filtro = [];

user_data.rec_names = rec_names(recording_indexes,:);
user_data.rec_idx = 1;
user_data.recording_ratios = ratios;
user_data.recording_indexes = recording_indexes;

user_data.bFilter = bFilter;
user_data.bBuildArtificial = bBuildArtificial;
user_data.Xoffset = 0.09;
user_data.Yoffset = 0;
user_data.Xscale = 1.05;
user_data.Yscale = 1.05;

user_data.tmp_path = tmp_path;

user_data.hb_detail_window = 3;

user_data = DoRecording(user_data);

set(user_data.fig_hdl, 'UserData', user_data);


function user_data = DoRecording(user_data)

    [~, rec_name] = fileparts(user_data.rec_names(user_data.rec_idx).name);

    user_data.rec_path = [user_data.tmp_path rec_name ];

    disp([num2str(user_data.recording_indexes(user_data.rec_idx)) ' : ' rec_name])

    if(user_data.bLoadECG)
        disp('Loading ...')
        ECG_struct = load(user_data.rec_path);
    end

    user_data.bAnnsEdited = false;
    user_data.bRecEdited = false;
    user_data.ECG_struct = ECG_struct;

    user_data.ECG_struct.header.recname = rec_name;
    user_data.hb_idx = 1;
    user_data.lead_idx = 1;
    user_data.bLockRRserie = false;
    user_data.bLockScatter = false;
    user_data.selected_hb_idx = [];
    user_data.RRscatter_hb_idx_hdl = [];   
    user_data.RRserie_hb_idx_hdl = [];   
    user_data.undo_buffer = [];
    user_data.undo_buffer_idx = 1;
    user_data.Pattern_hdl = [];

    [~, ~, user_data.inv_dower_idx] = intersect({'I','II','V1','V2','V3','V4','V5','V6'}, cellstr(ECG_struct.header.desc));
    if( length(user_data.inv_dower_idx) == 8 )
        user_data.bUseDower = true;
    else
        user_data.inv_dower_idx = 1:ECG_struct.header.nsig;
        user_data.bUseDower = false;
    end

    user_data.ratios = [];
    user_data.AnnNames = [];

    if( user_data.bFilter && isempty(user_data.filtro) )
        user_data.filtro = bandpass_filter_design( ECG_struct.header.freq );
    end
    
    if( isfield(ECG_struct, 'noise_power') )
        noise_power = ECG_struct.noise_power;
    else
        noise_power = [];
    end

    if( isfield(user_data.ECG_struct, 'series_quality' ) ) 
        
        user_data.AnnNames = user_data.ECG_struct.series_quality.AnnNames;
        user_data.all_annotations = user_data.ECG_struct.series_quality.all_annotations;
        user_data.ratios = user_data.ECG_struct.series_quality.ratios;
        user_data.estimated_labs = user_data.ECG_struct.series_quality.estimated_labs;
        
    else
        
        disp('calculating series quality ...')

        user_data = calc_ratios(user_data);
        
    end
    
    [user_data.ratios, user_data.best_detections_idx] = sort(user_data.ratios, 'descend');
    
    if( user_data.bBuildArtificial )

        % generate artificial annotations combining K best annotations
        aux_idx = user_data.best_detections_idx(1:10);
        artificial_annotations = combine_anns(user_data.all_annotations(aux_idx), user_data.estimated_labs(aux_idx), user_data.ECG_struct.header);

        for ii = length(artificial_annotations):-1:1
            aux_str = ['artificial_' num2str(ii)];
            user_data.ECG_struct.(aux_str) = artificial_annotations(ii);
        end

        user_data = calc_ratios(user_data);

        [user_data.ratios, user_data.best_detections_idx] = sort(user_data.ratios, 'descend');
    
    end
    
    aux_val = 1:length(user_data.ratios);
    
    [~, user_data.annotations_ranking] = sort(aux_val(user_data.best_detections_idx));
    
    %     if( isfield(user_data.ECG_struct, 'corrected_anns' ) ) 
    %         user_data.AnnNames_idx = find(user_data.best_detections_idx == find(strcmp(user_data.AnnNames(:,1), 'corrected_anns') ));
    %         if( isfield(user_data.ECG_struct.corrected_anns, 'time' ) ) 
    %             user_data.anns_under_edition = unique(colvec(user_data.ECG_struct.corrected_anns.time));
    %         else
    %             user_data.anns_under_edition = unique(colvec(user_data.ECG_struct.corrected_anns));
    %         end
    %     else
    %         user_data.AnnNames_idx = 1;
    %         user_data.anns_under_edition = unique(round(colvec( user_data.ECG_struct.(user_data.AnnNames{user_data.AnnNames_idx,1}).(user_data.AnnNames{user_data.AnnNames_idx,2}) )));
    %     end

    user_data.AnnNames_idx = 1;
    if( isempty(user_data.AnnNames) )
        user_data.anns_under_edition = [];
    else
        user_data.anns_under_edition = unique(round(colvec( user_data.ECG_struct.(user_data.AnnNames{user_data.AnnNames_idx,1}).(user_data.AnnNames{user_data.AnnNames_idx,2}) )));
        disp( ['Using ' user_data.AnnNames{user_data.AnnNames_idx,1} ' annotations (' num2str(user_data.ratios(user_data.AnnNames_idx)) ')'] );
    end
    
    disp('drawing ...')

    user_data.bFirstLoad = true;
    
    user_data = Redraw(user_data);

    set(user_data.annotation_list_control, 'Value', 1);
    set(user_data.leads_control, 'Value', 1);

    
    function update_q_ratios(obj,event_obj) 

        fig_hnd = gcbf; user_data = get(fig_hnd,'userdata');

        user_data = calc_ratios(user_data);

        cant_anns = size(user_data.AnnNames,1);
        aux_str = repmat( ' - ',cant_anns,1);
        
        user_data.annotation_list_control = uicontrol( ... 
                      'style','listbox', ...
                      'units','normalized', ...
                      'string', [ char(cellstr(num2str((1:cant_anns)'))) aux_str char(user_data.AnnNames(:,1)) repmat( ' (',cant_anns,1) num2str(round(colvec(user_data.ratios * 1000))) aux_str num2str(colvec(user_data.annotations_ranking))  repmat( ')',cant_anns,1)  ] , ...
                      'position', [0.865 0.11 0.13 0.18] , ...
                      'min', 2, ...
                      'max', 4, ...
                      'callback', @ChangeAnnotationsSelected);
    
    function user_data = calc_ratios(user_data)
        
        user_data.AnnNames = [];
        
        for fname = rowvec(fieldnames(user_data.ECG_struct))
            if( isfield(user_data.ECG_struct.(fname{1}), 'time') )
                user_data.AnnNames = [user_data.AnnNames; cellstr(fname{1}) cellstr('time')];
            end
            if( isfield(user_data.ECG_struct.(fname{1}), 'qrs') )
                user_data.AnnNames = [user_data.AnnNames; cellstr(fname{1}) cellstr('qrs')];
            end
        end

        cant_anns = size(user_data.AnnNames,1);

        aux_val = cell(cant_anns,1);
        for ii = 1:cant_anns
            aux_val{ii} = user_data.ECG_struct.(user_data.AnnNames{ii,1}).(user_data.AnnNames{ii,2});
        end
        user_data.all_annotations = aux_val;

    %     user_data.ratios = CalcRRserieQuality(user_data.ECG_struct.signal, user_data.ECG_struct.header, user_data.all_annotations);
        [ user_data.ratios, user_data.estimated_labs ] = CalcRRserieRatio(user_data.all_annotations, user_data.ECG_struct.header);
    
        [user_data.ratios, user_data.best_detections_idx] = sort(user_data.ratios, 'descend');
        
        aux_val = 1:cant_anns;

        [~, user_data.annotations_ranking] = sort(aux_val(user_data.best_detections_idx));
    
        user_data.AnnNames = user_data.AnnNames(user_data.best_detections_idx,:);
        user_data.all_annotations = user_data.all_annotations(user_data.best_detections_idx);
        user_data.estimated_labs = user_data.estimated_labs(user_data.best_detections_idx);
    
        

function user_data = Redraw(user_data)

    user_data.RRserie = colvec(diff(user_data.anns_under_edition));
    if( isempty(user_data.RRserie) )
        limits = [0.6 0.7];
    else        
        user_data.RRserie = [user_data.RRserie(1); user_data.RRserie] * 1/user_data.ECG_struct.header.freq;
        limits = prctile(user_data.RRserie, [2.5 97.5]);
    end

    aux_val = diff(limits);
    aux_val = max(0.01, 0.2*aux_val);
    limits(1) = limits(1) - aux_val;
    limits(2) = limits(2) + aux_val;
        
    if( user_data.bLockScatter )
        x_lims_scatter = get(user_data.Scatter_axes_hdl, 'Xlim');
        y_lims_scatter = get(user_data.Scatter_axes_hdl, 'Ylim');
    end
    
    if( user_data.bLockRRserie )
        x_lims_RRserie = get(user_data.RRserie_axes_hdl, 'Xlim');
        y_lims_RRserie = get(user_data.RRserie_axes_hdl, 'Ylim');
    end
    
    user_data.fig_hdl = figure(1);
    
    if( ~isfield(user_data, 'maximized_size') )
        maximize(user_data.fig_hdl);
        user_data.maximized_size = get(user_data.fig_hdl, 'Position');
    end
    
    set(user_data.fig_hdl, 'Position', [ user_data.maximized_size(3:4) user_data.maximized_size(3:4) ] .* [ 0.05 0.13 0.95 0.9] );
    
%     set(user_data.fig_hdl, 'Toolbar', 'none');

    %% Axis de Poincaré

    if( ~isfield(user_data, 'Scatter_axes_hdl') )
        user_data.Scatter_axes_hdl = axes('Position', ([0.55 0.41 0.355 0.52 ] - [ user_data.Xoffset user_data.Yoffset 0 0 ]) .* [ 1 1 user_data.Xscale user_data.Yscale ] );
    end
    
    if( isempty(user_data.RRserie) )
        cla(user_data.Scatter_axes_hdl)
        RRscatter_hdl = [];
    else
        RRscatter_hdl = plot(user_data.Scatter_axes_hdl, user_data.RRserie, [ user_data.RRserie(2:end); user_data.RRserie(end) ], 'bo' );
    end

    if( ~isempty(user_data.RRserie) )
        hold(user_data.Scatter_axes_hdl, 'on')
        RRnext = [user_data.RRserie(2:end); user_data.RRserie(end)];
        user_data.RRscatter_selection_hdl = plot(user_data.Scatter_axes_hdl, user_data.RRserie(user_data.selected_hb_idx), RRnext(user_data.selected_hb_idx), 'xg' );
        hold(user_data.Scatter_axes_hdl, 'off')
    end    
    
    if( user_data.bLockScatter )
        xlim(user_data.Scatter_axes_hdl, x_lims_scatter);
        ylim(user_data.Scatter_axes_hdl, y_lims_scatter);
    else
        xlim(user_data.Scatter_axes_hdl, limits);
        ylim(user_data.Scatter_axes_hdl, limits);
%         zoom reset
    end
    
    title(user_data.Scatter_axes_hdl, ['Poincaré plot - Registro ' user_data.ECG_struct.header.recname ]);
    xlabel(user_data.Scatter_axes_hdl, 'RR current')
    ylabel(user_data.Scatter_axes_hdl, 'RR next')

    
    %% Axis de RR serie
    
    if( ~isfield(user_data, 'RRserie_axes_hdl') )
        user_data.RRserie_axes_hdl = axes('Position', ([0.13 0.709 0.375 0.216 ] - [ user_data.Xoffset user_data.Yoffset 0 0 ]) .* [ 1 1 user_data.Xscale user_data.Yscale ] );
    end
    
    if( isempty(user_data.RRserie) )
        cla(user_data.RRserie_axes_hdl)
        RRserie_hdl = [];
    else
        RRserie_hdl = plot(user_data.RRserie_axes_hdl, user_data.anns_under_edition, user_data.RRserie, '-xb' );
    end
    
    if( ~isempty(user_data.RRserie) )
        hold(user_data.RRserie_axes_hdl, 'on')
        user_data.RRserie_selection_hdl = plot(user_data.RRserie_axes_hdl, user_data.anns_under_edition(user_data.selected_hb_idx), user_data.RRserie(user_data.selected_hb_idx), 'og' );
        hold(user_data.RRserie_axes_hdl, 'off')
    end
    
    if( user_data.bLockRRserie )
        xlim(user_data.RRserie_axes_hdl, x_lims_RRserie);
        ylim(user_data.RRserie_axes_hdl, y_lims_RRserie);
    else
        if( ~isempty(user_data.RRserie) )
            ylim(user_data.RRserie_axes_hdl, limits);
            x_range = user_data.anns_under_edition(end) - user_data.anns_under_edition(1);
            xlim(user_data.RRserie_axes_hdl, [user_data.anns_under_edition(1) - 0.05 * x_range user_data.anns_under_edition(end) + 0.05 * x_range ] );
%         zoom reset
        end
    end
    
    if( length(user_data.RRserie) > 5 )
        cant_ticks = 5;
        aux_idx = round(linspace(0, user_data.anns_under_edition(end), cant_ticks));
        set(user_data.RRserie_axes_hdl, 'XTick', rowvec(aux_idx) );
        aux_str = cellstr(Seconds2HMS(aux_idx*1/user_data.ECG_struct.header.freq));    
        set(user_data.RRserie_axes_hdl, 'XTickLabel', char(aux_str) );
    end
    xlabel(user_data.RRserie_axes_hdl, 'Time')
    ylabel(user_data.RRserie_axes_hdl, 'RR interval')

    %% Axis de RR serie ZOOM

    if( ~isfield(user_data, 'RRserie_zoom_axes_hdl') )
        user_data.RRserie_zoom_axes_hdl = axes('Position', ([0.13 0.41 0.375 0.216 ]- [ user_data.Xoffset user_data.Yoffset 0 0 ]) .* [ 1 1 user_data.Xscale user_data.Yscale ] );
    end

    if( isempty(user_data.RRserie) )
        cla(user_data.RRserie_zoom_axes_hdl)
    else
        UpdateRRserieZoom(user_data);
    end
    
    %% Axis de ECG
    
    if( ~isfield(user_data, 'ECG_axes_hdl') )
        user_data.ECG_axes_hdl = axes('Position', ([0.13 0.11 0.775 0.216] - [ user_data.Xoffset user_data.Yoffset 0 0 ]) .* [ 1 1 user_data.Xscale user_data.Yscale ] );
    end

    if( isfield(user_data, 'annotation_list_control') )
        anns_selected = get(user_data.annotation_list_control, 'Value');
    else
        anns_selected = [];
    end
    
    if( length(anns_selected) > 1 )
        
        aux_anns = {user_data.anns_under_edition};
        anns_selected( anns_selected == user_data.AnnNames_idx) = [];
        
        for ii = anns_selected
            aux_anns = [aux_anns; ... 
                        {user_data.ECG_struct.(user_data.AnnNames{ii,1}).(user_data.AnnNames{ii,2})}
                        ];
        end
        
    else
        aux_anns = user_data.anns_under_edition;
    end
    
    user_data.ECG_hdl = plot_ecg_heartbeat(user_data.ECG_struct.signal(:,user_data.lead_idx), aux_anns, user_data.hb_idx , user_data.hb_detail_window, user_data.ECG_struct.header, user_data.filtro, user_data.ECG_axes_hdl);    
    
    if( length(user_data.lead_idx) > 1 )
        title(user_data.ECG_axes_hdl, ['Heartbeat ' num2str(user_data.hb_idx) ' : All Leads' ] )
    else
        title(user_data.ECG_axes_hdl, ['Heartbeat ' num2str(user_data.hb_idx) ' : Lead ' user_data.ECG_struct.header.desc(user_data.lead_idx,:)] )
    end
    
    xlabel(user_data.ECG_axes_hdl, 'Sample #')
    ylabel(user_data.ECG_axes_hdl, 'ECG')
%     zoom reset

    %% Controls
    
    if( ~isfield(user_data, 'annotation_list_control') || user_data.bFirstLoad )

        user_data.bFirstLoad = false;
        
        %% Recording list
                
        aux_str = repmat( ' - ', length(user_data.recording_indexes),1);
        
        user_data.recordings_control = uicontrol( ... 
                      'style','listbox', ...
                      'units','normalized', ...
                      'string', [ char(cellstr(num2str(colvec(user_data.recording_indexes)))) aux_str num2str(round(user_data.recording_ratios*1000)) aux_str char(user_data.rec_names.name) ] , ...
                      'position', [0.865 0.75 0.13 0.2] , ...
                      'min', 1, ...
                      'max', 1, ...
                      'Value', user_data.rec_idx, ...
                      'callback', @ChangeRecordingSelected);
                  
        user_data.recordings_label = uicontrol( ...
                      'style','text', ...
                      'string', 'Available Recordings', ...
                      'units','normalized', ...
                      'position', [0.865 0.96 0.13 0.025] );  
                  
        %% Signal list
        user_data.leads_control = uicontrol( ... 
                      'style','listbox', ...
                      'units','normalized', ...
                      'string', [ char(cellstr(num2str((1:user_data.ECG_struct.header.nsig)'))) repmat( ' - ',user_data.ECG_struct.header.nsig,1) user_data.ECG_struct.header.desc ] , ...
                      'position', [0.865 0.4 0.13 0.3] , ...
                      'min', 2, ...
                      'max', 4, ...
                      'Value', user_data.lead_idx, ...
                      'callback', @ChangeLeadsSelected);
                  
        user_data.leads_label = uicontrol( ...
                      'style','text', ...
                      'string', 'Available signals', ...
                      'units','normalized', ...
                      'position', [0.865 0.71 0.13 0.025] );        
                          
        %% Annotation list
                  
        cant_anns = size(user_data.AnnNames,1);
        
        aux_str = repmat( ' - ',cant_anns,1);
        
        user_data.annotation_list_control = uicontrol( ... 
                      'style','listbox', ...
                      'units','normalized', ...
                      'string', [ char(cellstr(num2str((1:cant_anns)'))) aux_str char(user_data.AnnNames(:,1)) repmat( ' (',cant_anns,1) num2str(round(colvec(user_data.ratios * 1000))) aux_str num2str(colvec(user_data.annotations_ranking))  repmat( ')',cant_anns,1)  ] , ...
                      'position', [0.865 0.11 0.13 0.18] , ...
                      'min', 2, ...
                      'max', 4, ...
                      'callback', @ChangeAnnotationsSelected);
                  
        user_data.annotation_list_label = uicontrol( ...
                      'style','text', ...
                      'string', 'Annotations available', ...
                      'units','normalized', ...
                      'position', [0.865 0.30 0.13 0.025] );
                  
        user_data.annotation_under_edition_label = uicontrol( ...
                      'style','text', ...
                      'string', [ 'Annotation under edition: ' char(user_data.AnnNames( user_data.AnnNames_idx ,1)) ' (' num2str(user_data.ratios(user_data.AnnNames_idx)) ')' ], ...
                      'units','normalized', ...
                      'position', [0.865 0.35 0.13 0.05] );
                  
        user_data.DeleteAnnotation = uicontrol( ... 
                        'style','pushbutton', ...
                        'string', 'Delete Annotations', ...
                        'units','normalized', ...
                        'position', [0.865 0.07 0.063 0.03], ...
                        'callback', @DeleteAnnotations);
                  
        user_data.DeleteAnnotation = uicontrol( ... 
                        'style','pushbutton', ...
                        'string', 'Upd Q ratios', ...
                        'units','normalized', ...
                        'position', [0.93 0.07 0.063 0.03], ...
                        'callback', @update_q_ratios);
                  
                  
    end
    

    user_data.KeyPressed = '';
    
    user_data.CurrAxesHdl = gca;

    if( ~isempty(RRscatter_hdl) )
        set(RRscatter_hdl,'ButtonDownFcn',@inspect_scatter);            
    end
    set(user_data.Scatter_axes_hdl,'ButtonDownFcn',@inspect_scatter);            
    if( ~isempty(RRserie_hdl) )
        set(RRserie_hdl, 'ButtonDownFcn',@inspect_RRserie);            
    end
    set(user_data.RRserie_axes_hdl,'ButtonDownFcn',@inspect_RRserie);            
    set(user_data.ECG_hdl,'ButtonDownFcn',@inspect_ECG);            
    set(user_data.ECG_axes_hdl,'ButtonDownFcn',@inspect_ECG);            
%     set(user_data.fig_hdl,'WindowScrollWheelFcn',@ScrollWheel);            
    set(user_data.fig_hdl,'WindowKeyPressFcn',@KeyPress);            

    
function inspect_scatter(obj,event_obj)

	fig_hnd = gcbf; 
    user_data = get(fig_hnd, 'UserData');
    
    if( strcmp(get(fig_hnd,'SelectionType'),'alt'))
        % Delete annotation
        
        point = get(gca,'CurrentPoint');

        [~, user_data.hb_idx] = min(sum( bsxfun(@minus, point(1,1:2), [user_data.RRserie, [user_data.RRserie(2:end); user_data.RRserie(end)]]).^2,2));
        
        user_data = PushUndoAction(user_data);
        
        user_data.anns_under_edition(user_data.hb_idx) = [];
        
        user_data.selected_hb_idx = [];
        
        user_data = Redraw(user_data);
        
        user_data.bAnnsEdited = true;
        user_data.bRecEdited = true;
        
    else
        
        if (strcmp(get(fig_hnd,'SelectionType'),'extend'))
            point1 = get(gca,'CurrentPoint');    % button down detected
            rbbox;
            point2 = get(gca,'CurrentPoint');    % button up detected            

            if( ~isempty(user_data.selected_hb_idx) )
                delete(user_data.RRserie_selection_hdl)
                delete(user_data.RRscatter_selection_hdl)
            end
            
            xlims = sort([ point1(1,1) point2(1,1) ]);
            ylims = sort([ point1(1,2) point2(1,2) ]);
            aux_RRcurr = user_data.RRserie;
            aux_RRnext = [user_data.RRserie(2:end); user_data.RRserie(end)];
            user_data.selected_hb_idx = find( aux_RRcurr > xlims(1) & aux_RRcurr < xlims(2) & aux_RRnext > ylims(1) & aux_RRnext < ylims(2) );

            disp([num2str(length(user_data.selected_hb_idx)) ' heartbeats selected.'])

            hold(user_data.Scatter_axes_hdl, 'on')
            RRnext = [user_data.RRserie(2:end); user_data.RRserie(end)];
            user_data.RRscatter_selection_hdl = plot(user_data.Scatter_axes_hdl, user_data.RRserie(user_data.selected_hb_idx), RRnext(user_data.selected_hb_idx), 'xg' );
            hold(user_data.Scatter_axes_hdl, 'off')

            hold(user_data.RRserie_axes_hdl, 'on')
            user_data.RRserie_selection_hdl = plot(user_data.RRserie_axes_hdl,user_data.anns_under_edition(user_data.selected_hb_idx) , user_data.RRserie(user_data.selected_hb_idx), 'og' );
            hold(user_data.RRserie_axes_hdl, 'off')
            
            min_hb_idx = min(user_data.selected_hb_idx);
            cant_hb_idx = max(user_data.selected_hb_idx) - min_hb_idx + 1;
            
            if( isfield(user_data, 'annotation_list_control') )
                anns_selected = get(user_data.annotation_list_control, 'Value');
            else
                anns_selected = [];
            end

            if( length(anns_selected) > 1 )

                aux_anns = {user_data.anns_under_edition};
                anns_selected( anns_selected == user_data.AnnNames_idx) = [];

                for ii = anns_selected
                    aux_anns = [aux_anns; ... 
                                {user_data.ECG_struct.(user_data.AnnNames{ii,1}).(user_data.AnnNames{ii,2})}
                                ];
                end

            else
                aux_anns = user_data.anns_under_edition;
            end

            user_data.ECG_hdl = plot_ecg_heartbeat(user_data.ECG_struct.signal(:,user_data.lead_idx), aux_anns, min_hb_idx , cant_hb_idx, user_data.ECG_struct.header, user_data.filtro, user_data.ECG_axes_hdl);    

            if( length(user_data.lead_idx) > 1 )
                title(user_data.ECG_axes_hdl, ['Heartbeat ' num2str(min_hb_idx) ' : All Leads' ] )
            else
                title(user_data.ECG_axes_hdl, ['Heartbeat ' num2str(min_hb_idx) ' : Lead ' user_data.ECG_struct.header.desc(user_data.lead_idx,:)] )
            end
            
%             zoom reset
            
            set(user_data.ECG_hdl,'ButtonDownFcn',@inspect_ECG);            
            set(user_data.ECG_axes_hdl,'ButtonDownFcn',@inspect_ECG);            
            
            % Zoom RR serie
            UpdateRRserieZoom(user_data);
            
        else
            % Show ECG
            point = get(gca,'CurrentPoint');

            if( ishandle(user_data.RRscatter_hb_idx_hdl) )
                delete(user_data.RRscatter_hb_idx_hdl)
            end
            if( ishandle(user_data.RRserie_hb_idx_hdl) )
                delete(user_data.RRserie_hb_idx_hdl)
            end
            
            [~, user_data.hb_idx] = min(sum( bsxfun(@minus, point(1,1:2), [user_data.RRserie(1:end-1), user_data.RRserie(2:end)]).^2,2));
            
            hold(user_data.Scatter_axes_hdl, 'on')
            user_data.RRscatter_hb_idx_hdl = plot(user_data.Scatter_axes_hdl, user_data.RRserie(user_data.hb_idx), user_data.RRserie(user_data.hb_idx+1), 'xr' );
            hold(user_data.Scatter_axes_hdl, 'off')

            hold(user_data.RRserie_axes_hdl, 'on')
            user_data.RRserie_hb_idx_hdl = plot(user_data.RRserie_axes_hdl, user_data.anns_under_edition(user_data.hb_idx) , user_data.RRserie(user_data.hb_idx), 'or' );
            hold(user_data.RRserie_axes_hdl, 'off')
            
            if( isfield(user_data, 'annotation_list_control') )
                anns_selected = get(user_data.annotation_list_control, 'Value');
            else
                anns_selected = [];
            end

            if( length(anns_selected) > 1 )

                aux_anns = {user_data.anns_under_edition};
                anns_selected( anns_selected == user_data.AnnNames_idx) = [];

                for ii = anns_selected
                    aux_anns = [aux_anns; ... 
                                {user_data.ECG_struct.(user_data.AnnNames{ii,1}).(user_data.AnnNames{ii,2})}
                                ];
                end

            else
                aux_anns = user_data.anns_under_edition;
            end

            user_data.ECG_hdl = plot_ecg_heartbeat(user_data.ECG_struct.signal(:,user_data.lead_idx), aux_anns, user_data.hb_idx , user_data.hb_detail_window , user_data.ECG_struct.header, user_data.filtro, user_data.ECG_axes_hdl);    

            if( length(user_data.lead_idx) > 1 )
                title(user_data.ECG_axes_hdl, ['Heartbeat ' num2str(user_data.hb_idx) ' : All Leads' ] )
            else
                title(user_data.ECG_axes_hdl, ['Heartbeat ' num2str(user_data.hb_idx) ' : Lead ' user_data.ECG_struct.header.desc(user_data.lead_idx,:)] )
            end
            
%             zoom reset
            
            set(user_data.ECG_hdl,'ButtonDownFcn',@inspect_ECG);            
            set(user_data.ECG_axes_hdl,'ButtonDownFcn',@inspect_ECG);            
            
            % Zoom RR serie
            UpdateRRserieZoom(user_data);
            
        end
        
    end

    user_data.KeyPressed = '';
    
    set(fig_hnd, 'UserData', user_data);

function inspect_ECG(obj,event_obj)


	fig_hnd = gcbf; 
    
%     disp(get(fig_hnd,'SelectionType'))
    
    if (strcmp(get(fig_hnd,'SelectionType'),'alt'))
        % Delete annotation
        user_data = get(fig_hnd, 'UserData');
        
        point = get(gca,'CurrentPoint');

        [~, user_data.hb_idx] = min(abs( point(1) - user_data.anns_under_edition));
        
        user_data = PushUndoAction(user_data);
        
        user_data.anns_under_edition(user_data.hb_idx) = [];
        
        user_data.selected_hb_idx( user_data.selected_hb_idx > length(user_data.anns_under_edition) | user_data.selected_hb_idx == user_data.hb_idx ) = [];
        
        if( isempty(user_data.anns_under_edition) )
            user_data.hb_idx = [];
        else
            user_data.hb_idx = max(1, user_data.hb_idx - 1); 
        end
        
        user_data = Redraw(user_data);

        user_data.bAnnsEdited = true;
        user_data.bRecEdited = true;
        
        
        set(fig_hnd, 'UserData', user_data);
        
    elseif (strcmp(get(fig_hnd,'SelectionType'),'normal'))

        % Show ECG
        user_data = get(fig_hnd, 'UserData');
        
        point = get(gca,'CurrentPoint');

        user_data = PushUndoAction(user_data);
        
        aux_val = round(point(1));
        
        user_data.anns_under_edition = sort([user_data.anns_under_edition; aux_val ]);
        
        user_data.selected_hb_idx = find(aux_val == user_data.anns_under_edition, 1);
        user_data.hb_idx = user_data.selected_hb_idx;
        
        user_data = Redraw(user_data);
        
        user_data.bAnnsEdited = true;
        user_data.bRecEdited = true;
        
        set(fig_hnd, 'UserData', user_data);

    elseif (strcmp(get(fig_hnd,'SelectionType'),'extend'))
    
        % Show ECG
        
        user_data = get(fig_hnd, 'UserData');
        point = get(gca,'CurrentPoint');
        aux_val = round(point(1));
        
        if( isempty(user_data.side_plot) || ~ishandle(user_data.side_plot) )
            user_data.side_plot = figure();
            set(user_data.side_plot, 'Position', [ user_data.maximized_size(3:4) user_data.maximized_size(3:4) ] .* [ 0 0 0.95 0.9] );
            set(fig_hnd, 'UserData', user_data);
        else
            figure(user_data.side_plot);
        end
        
        plot_ecg_strip(user_data.ECG_struct.signal, ...
                        'ECG_header', user_data.ECG_struct.header, ...
                        'QRS_annotations', user_data.anns_under_edition, ...
                        'Start_time', aux_val/user_data.ECG_struct.header.freq - 1 , ...
                        'End_time', aux_val/user_data.ECG_struct.header.freq + 1 );

        figure(fig_hnd);
                    
    end
    
function ScrollWheel(obj,event_obj)

    user_data = get(obj, 'UserData');

    set(user_data.fig_hdl, 'CurrentAxes', user_data.CurrAxesHdl)
    
    point = get(user_data.CurrAxesHdl,'CurrentPoint');

    x_lims = get(user_data.CurrAxesHdl, 'Xlim');
    y_lims = get(user_data.CurrAxesHdl, 'Ylim');

    movement = point(1,1:2) - [ mean(x_lims) mean(y_lims) ];

    set(user_data.CurrAxesHdl, 'Xlim', x_lims + movement(1));
    set(user_data.CurrAxesHdl, 'Ylim', y_lims + movement(2));
    
    z_hdl = zoom(obj);
%     setAllowAxesZoom(z_hdl,[user_data.Scatter_axes_hdl user_data.RRserie_axes_hdl user_data.ECG_axes_hdl ], false)
%     
%     setAllowAxesZoom(z_hdl,user_data.CurrAxesHdl, true)
    
    if( strcmp(user_data.KeyPressed, 'control' ) )
        set(z_hdl, 'Motion', 'horizontal')
        zoom( 1-(event_obj.VerticalScrollCount*event_obj.VerticalScrollAmount*0.1) );
    elseif( strcmp(user_data.KeyPressed, 'shift' ) )
        set(z_hdl, 'Motion', 'vertical')
        zoom( 1-(event_obj.VerticalScrollCount*event_obj.VerticalScrollAmount*0.1) );
    else
        set(z_hdl, 'Motion', 'both')
        zoom( 1-(event_obj.VerticalScrollCount*event_obj.VerticalScrollAmount*0.1) );
    end
 
function KeyPress(obj,event_obj)

    user_data = get(obj, 'UserData');
    
    if (any(strcmp(event_obj.Key,{'control' 'alt' 'shift'})))
        user_data.KeyPressed = event_obj.Key;
        
    elseif (strcmp(event_obj.Key,'c') && strcmp(event_obj.Modifier,'control'))
        % copy selection
        if( ~isempty(user_data.selected_hb_idx) )
            user_data.copy_paste_buffer = user_data.anns_under_edition(user_data.selected_hb_idx);
            disp('Selection copied')
        end
        
    elseif (strcmp(event_obj.Key,'v') && strcmp(event_obj.Modifier,'control'))
        % paste selection
        if( ~isempty(user_data.copy_paste_buffer) )
            user_data.bAnnsEdited = true;
            user_data.bRecEdited = true;
            user_data = PushUndoAction(user_data);
            user_data.anns_under_edition = sort( unique([user_data.anns_under_edition; colvec(user_data.copy_paste_buffer) ]) );
            user_data.selected_hb_idx = [];
            user_data = Redraw(user_data);
        end
        
    elseif (strcmp(event_obj.Key,'delete'))
        % delete selection
        if( ~isempty(user_data.selected_hb_idx) )
            user_data.bAnnsEdited = true;
            user_data.bRecEdited = true;
            user_data = PushUndoAction(user_data);
            if( isempty(user_data.hb_idx) )
                user_data.hb_idx = 1;
            else
                user_data.hb_idx = find(user_data.anns_under_edition < user_data.anns_under_edition(user_data.hb_idx), 1, 'first'  );
            end
            user_data.anns_under_edition(user_data.selected_hb_idx) = [];
            user_data.selected_hb_idx = [];
            user_data = Redraw(user_data);
        end
        
    elseif (strcmp(event_obj.Key,'y') && ~isempty(event_obj.Modifier) && strcmp(event_obj.Modifier,'control'))
        % redo last change
        if( user_data.undo_buffer_idx < length(user_data.undo_buffer) )
            user_data.undo_buffer_idx = user_data.undo_buffer_idx + 1;
            user_data.anns_under_edition = user_data.undo_buffer{user_data.undo_buffer_idx};
            user_data.selected_hb_idx = [];
            user_data = Redraw(user_data);
        end
        
    elseif (strcmp(event_obj.Key,'z') && strcmp(event_obj.Modifier,'control'))
        % undo last change
        if( user_data.undo_buffer_idx > 1 )
            user_data.undo_buffer_idx = user_data.undo_buffer_idx - 1;
            user_data.anns_under_edition = user_data.undo_buffer{user_data.undo_buffer_idx};
            user_data.selected_hb_idx = [];
            user_data = Redraw(user_data);
        end

    elseif (strcmp(event_obj.Key,'rightarrow') && ~isempty(event_obj.Modifier) && strcmp(event_obj.Modifier,'control'))
    
        user_data = get(user_data.fig_hdl, 'UserData');

        if(user_data.bRecEdited)

            if( user_data.bLoadECG )
                disp('Saving data ...')
                user_data.ECG_struct.(user_data.AnnNames{user_data.AnnNames_idx,1}).(user_data.AnnNames{user_data.AnnNames_idx,2}) = unique(round(colvec( user_data.anns_under_edition )));
                ECG_struct = user_data.ECG_struct;
                save(user_data.rec_path, '-struct', 'ECG_struct');        
                disp(['Saved ' user_data.rec_path])
            else
                disp('Returning data ...')
                ann_output = user_data.anns_under_edition; 
            end
        end

        if( user_data.rec_idx >= 1 && user_data.rec_idx < length(user_data.recording_indexes) )
            user_data.rec_idx = user_data.rec_idx + 1;
        else
            user_data.rec_idx = 1;
        end
        
        set(user_data.recordings_control, 'Value', user_data.rec_idx );
        
        user_data = DoRecording(user_data);
        
    elseif (strcmp(event_obj.Key,'leftarrow') && ~isempty(event_obj.Modifier) && strcmp(event_obj.Modifier,'control'))

        user_data = get(user_data.fig_hdl, 'UserData');

        if(user_data.bRecEdited)

            if( user_data.bLoadECG )
                disp('Saving data ...')
                user_data.ECG_struct.(user_data.AnnNames{user_data.AnnNames_idx,1}).(user_data.AnnNames{user_data.AnnNames_idx,2}) = unique(round(colvec( user_data.anns_under_edition )));
                ECG_struct = user_data.ECG_struct;
                save(user_data.rec_path, '-struct', 'ECG_struct');        
                disp(['Saved ' user_data.rec_path])
            else
                disp('Returning data ...')
                ann_output = user_data.anns_under_edition; 
            end
        end

        if( user_data.rec_idx > 1 && user_data.rec_idx <= length(user_data.recording_indexes) )
            user_data.rec_idx = user_data.rec_idx - 1;
        else
            user_data.rec_idx = length(user_data.recording_indexes);
        end
        
        set(user_data.recordings_control, 'Value', user_data.rec_idx );
        
        user_data = DoRecording(user_data);

    elseif (strcmp(event_obj.Key,'p') )
        
        disp('Click and drag the time interval where the pattern is ...')
        set(obj,'CurrentAxes', user_data.ECG_axes_hdl);
        waitforbuttonpress;
        user_data = SearchPattern(user_data);
        
    elseif ( ~isempty(event_obj.Modifier) && strcmp(event_obj.Key,'g') && strcmp(event_obj.Modifier,'control'))
        
        if( user_data.bLoadECG )
            disp('Saving data ...')
            user_data = get(obj, 'UserData');
            user_data.ECG_struct.(user_data.AnnNames{user_data.AnnNames_idx,1}).(user_data.AnnNames{user_data.AnnNames_idx,2}) = unique(round(colvec( user_data.anns_under_edition )));
            ECG_struct = user_data.ECG_struct;
            save(user_data.rec_path, '-struct', 'ECG_struct');        
            disp(['Saved ' user_data.rec_path])
        else
            disp('Cant save in this mode.')
        end
        
    elseif (strcmp(event_obj.Key,'1'))
        user_data.CurrAxesHdl = user_data.Scatter_axes_hdl;
        set(user_data.fig_hdl, 'CurrentAxes', user_data.CurrAxesHdl)
        user_data.bLockScatter = true;
    elseif (strcmp(event_obj.Key,'2'))
        user_data.CurrAxesHdl = user_data.RRserie_axes_hdl;
        set(user_data.fig_hdl, 'CurrentAxes', user_data.CurrAxesHdl)
        user_data.bLockRRserie = true;
    elseif (strcmp(event_obj.Key,'3'))
        user_data.CurrAxesHdl = user_data.ECG_axes_hdl;
        set(user_data.fig_hdl, 'CurrentAxes', user_data.CurrAxesHdl)
    elseif (strcmp(event_obj.Key,'t'))
        %% Toggle ECG lead
        if( length(user_data.lead_idx) > 1 )
            user_data.lead_idx = 1;
        else
            if( user_data.lead_idx < user_data.ECG_struct.header.nsig )
                user_data.lead_idx = user_data.lead_idx + 1;
            else
                user_data.lead_idx = 1:user_data.ECG_struct.header.nsig;
            end
        end
        cla(user_data.ECG_axes_hdl);
        
        if( isfield(user_data, 'annotation_list_control') )
            anns_selected = get(user_data.annotation_list_control, 'Value');
        else
            anns_selected = [];
        end

        if( length(anns_selected) > 1 )

            aux_anns = {user_data.anns_under_edition};
            anns_selected( anns_selected == user_data.AnnNames_idx) = [];

            for ii = anns_selected
                aux_anns = [aux_anns; ... 
                            {user_data.ECG_struct.(user_data.AnnNames{ii,1}).(user_data.AnnNames{ii,2})}
                            ];
            end

        else
            aux_anns = user_data.anns_under_edition;
        end
            
        user_data.ECG_hdl = plot_ecg_heartbeat(user_data.ECG_struct.signal(:,user_data.lead_idx), aux_anns, user_data.hb_idx , user_data.hb_detail_window , user_data.ECG_struct.header, user_data.filtro, user_data.ECG_axes_hdl);    
        
        if( length(user_data.lead_idx) > 1 )
            title(user_data.ECG_axes_hdl, ['Heartbeat ' num2str(user_data.hb_idx) ' : All Leads' ] )
        else
            title(user_data.ECG_axes_hdl, ['Heartbeat ' num2str(user_data.hb_idx) ' : Lead ' user_data.ECG_struct.header.desc(user_data.lead_idx,:)] )
        end
        
%         zoom reset
        set(user_data.ECG_hdl,'ButtonDownFcn',@inspect_ECG);            
        set(user_data.ECG_axes_hdl,'ButtonDownFcn',@inspect_ECG);            
        
    elseif (strcmp(event_obj.Key,'rightarrow'))
        %% right arrow

        if(user_data.bAnnsEdited)
            user_data = update_annotations(user_data);
        end

        user_data.AnnNames_idx = user_data.AnnNames_idx + 1;
        if( user_data.AnnNames_idx > length(user_data.AnnNames) )
            user_data.AnnNames_idx = 1;
        end

        disp( ['Using ' user_data.AnnNames{user_data.AnnNames_idx,1} ' annotations (' num2str(user_data.ratios(user_data.AnnNames_idx)) ')'] );

        user_data.undo_buffer_idx = 1;

        user_data.anns_under_edition = unique(round(colvec( user_data.ECG_struct.(user_data.AnnNames{user_data.AnnNames_idx,1}).(user_data.AnnNames{user_data.AnnNames_idx,2}) )));

        user_data.bAnnsEdited = false;
        
        user_data.selected_hb_idx = [];

        user_data = Redraw(user_data);
            
        
    elseif (strcmp(event_obj.Key,'leftarrow'))
        %% left arrow

        if(user_data.bAnnsEdited)
            user_data = update_annotations(user_data);
        end
        
        user_data.AnnNames_idx = user_data.AnnNames_idx - 1;
        if( user_data.AnnNames_idx < 1  )
            user_data.AnnNames_idx = length(user_data.AnnNames);
        end

        disp( ['Using ' user_data.AnnNames{user_data.AnnNames_idx,1} ' annotations (' num2str(user_data.ratios(user_data.AnnNames_idx)) ')'] );

        user_data.undo_buffer_idx = 1;
        
        user_data.anns_under_edition = unique(round(colvec( user_data.ECG_struct.(user_data.AnnNames{user_data.AnnNames_idx,1}).(user_data.AnnNames{user_data.AnnNames_idx,2}) )));

        user_data.bAnnsEdited = false;
        
        user_data.selected_hb_idx = [];

        user_data = Redraw(user_data);
        
    else
        user_data.KeyPressed = '';
        user_data.bLockRRserie = false;
        user_data.bLockScatter = false;
        user_data.lead_idx = 1:user_data.ECG_struct.header.nsig;
    end
    
    set(obj, 'UserData', user_data);
    
function user_data = PushUndoAction(user_data)

    if( user_data.undo_buffer_idx < length(user_data.undo_buffer) )
        %me cargo todo lo que está por rehacerse
        user_data.undo_buffer = user_data.undo_buffer(1:user_data.undo_buffer_idx);
    end
    
    user_data.undo_buffer{user_data.undo_buffer_idx} = user_data.anns_under_edition;
    user_data.undo_buffer_idx = user_data.undo_buffer_idx + 1;

    
function inspect_RRserie(obj,event_obj)

	fig_hnd = gcbf; 
    
    if (strcmp(get(fig_hnd,'SelectionType'),'extend'))
        % Show ECG
        user_data = get(fig_hnd, 'UserData');
        
        point1 = get(gca,'CurrentPoint');    % button down detected
        rbbox;
        point2 = get(gca,'CurrentPoint');    % button up detected            

        if( ~isempty(user_data.selected_hb_idx) )
            delete(user_data.RRserie_selection_hdl)
            delete(user_data.RRscatter_selection_hdl)
        end
        
        [~, aux_val] = sort( abs(point1(1,1) - user_data.anns_under_edition) );
        xlims = aux_val(1);
        [~, aux_val] = sort( abs(point2(1,1) - user_data.anns_under_edition) );
        xlims = sort([xlims aux_val(1)]);
        
        ylims = sort([ point1(1,2) point2(1,2) ]);
        lanns_under_ed = length(user_data.anns_under_edition);
        bAux = false(lanns_under_ed,1);
        bAux(max(1,xlims(1)):min(lanns_under_ed,xlims(2))) = true;
        user_data.selected_hb_idx = find( bAux & user_data.RRserie > ylims(1) & user_data.RRserie < ylims(2) );

        if( isempty(user_data.selected_hb_idx) )
            % probar en el rango de x solamente
            user_data.selected_hb_idx = find( bAux );
        end
        
        if( ~isempty(user_data.selected_hb_idx) )
            
            disp([num2str(length(user_data.selected_hb_idx)) ' heartbeats selected.'])

            hold(user_data.Scatter_axes_hdl, 'on')
            RRnext = [user_data.RRserie(2:end); user_data.RRserie(end)];
            user_data.RRscatter_selection_hdl = plot(user_data.Scatter_axes_hdl, user_data.RRserie(user_data.selected_hb_idx), RRnext(user_data.selected_hb_idx), 'xg' );
            hold(user_data.Scatter_axes_hdl, 'off')

            hold(user_data.RRserie_axes_hdl, 'on')
            user_data.RRserie_selection_hdl = plot(user_data.RRserie_axes_hdl, user_data.anns_under_edition(user_data.selected_hb_idx) , user_data.RRserie(user_data.selected_hb_idx), 'og' );
            hold(user_data.RRserie_axes_hdl, 'off')

            min_hb_idx = min(user_data.selected_hb_idx);
            cant_hb_idx = max(user_data.selected_hb_idx) - min_hb_idx + 1;
            
            if( isfield(user_data, 'annotation_list_control') )
                anns_selected = get(user_data.annotation_list_control, 'Value');
            else
                anns_selected = [];
            end

            if( length(anns_selected) > 1 )

                aux_anns = {user_data.anns_under_edition};
                anns_selected( anns_selected == user_data.AnnNames_idx) = [];

                for ii = anns_selected
                    aux_anns = [aux_anns; ... 
                                {user_data.ECG_struct.(user_data.AnnNames{ii,1}).(user_data.AnnNames{ii,2})}
                                ];
                end

            else
                aux_anns = user_data.anns_under_edition;
            end
            
            user_data.ECG_hdl = plot_ecg_heartbeat(user_data.ECG_struct.signal(:,user_data.lead_idx), aux_anns, min_hb_idx , cant_hb_idx, user_data.ECG_struct.header, user_data.filtro, user_data.ECG_axes_hdl);    

            if( length(user_data.lead_idx) > 1 )
                title(user_data.ECG_axes_hdl, ['Heartbeat ' num2str(min_hb_idx) ' : All Leads' ] )
            else
                title(user_data.ECG_axes_hdl, ['Heartbeat ' num2str(min_hb_idx) ' : Lead ' user_data.ECG_struct.header.desc(user_data.lead_idx,:)] )
            end
            
%             zoom reset
            
            set(user_data.ECG_hdl,'ButtonDownFcn',@inspect_ECG);            
            set(user_data.ECG_axes_hdl,'ButtonDownFcn',@inspect_ECG);            
            
            % Zoom RR serie
            UpdateRRserieZoom(user_data);
            
            set(fig_hnd, 'UserData', user_data);
        
        end
        
    elseif (strcmp(get(fig_hnd,'SelectionType'),'normal'))
        
        % Show ECG
        user_data = get(fig_hnd, 'UserData');

        if( ishandle(user_data.RRscatter_hb_idx_hdl) )
            delete(user_data.RRscatter_hb_idx_hdl)
        end
        if( ishandle(user_data.RRserie_hb_idx_hdl) )
            delete(user_data.RRserie_hb_idx_hdl)
        end
        
        point = get(gca,'CurrentPoint');
        
        lanns_under_ed = length(user_data.anns_under_edition);
        [~, aux_val] = sort( abs(point(1) - user_data.anns_under_edition) );
        aux_val = aux_val(1);
        
        user_data.hb_idx = max(1, min(lanns_under_ed, aux_val ));

        hold(user_data.Scatter_axes_hdl, 'on')
        RRnext = [user_data.RRserie(2:end); user_data.RRserie(end)];
        user_data.RRscatter_selection_hdl = plot(user_data.Scatter_axes_hdl, user_data.RRserie(user_data.selected_hb_idx), RRnext(user_data.selected_hb_idx), 'xg' );
        hold(user_data.Scatter_axes_hdl, 'off')

        hold(user_data.RRserie_axes_hdl, 'on')
        user_data.RRserie_hb_idx_hdl = plot(user_data.RRserie_axes_hdl, user_data.anns_under_edition(user_data.hb_idx), user_data.RRserie(user_data.hb_idx), 'or' );
        hold(user_data.RRserie_axes_hdl, 'off')        

        %Zoom RR serie
        UpdateRRserieZoom(user_data);

        if( isfield(user_data, 'annotation_list_control') )
            anns_selected = get(user_data.annotation_list_control, 'Value');
        else
            anns_selected = [];
        end

        if( length(anns_selected) > 1 )

            aux_anns = {user_data.anns_under_edition};
            anns_selected( anns_selected == user_data.AnnNames_idx) = [];

            for ii = anns_selected
                aux_anns = [aux_anns; ... 
                            {user_data.ECG_struct.(user_data.AnnNames{ii,1}).(user_data.AnnNames{ii,2})}
                            ];
            end

        else
            aux_anns = user_data.anns_under_edition;
        end

        user_data.ECG_hdl = plot_ecg_heartbeat(user_data.ECG_struct.signal(:,user_data.lead_idx), aux_anns, user_data.hb_idx , user_data.hb_detail_window , user_data.ECG_struct.header, user_data.filtro, user_data.ECG_axes_hdl);    
        
        if( length(user_data.lead_idx) > 1 )
            title(user_data.ECG_axes_hdl, ['Heartbeat ' num2str(user_data.hb_idx) ' : All Leads' ] )
        else
            title(user_data.ECG_axes_hdl, ['Heartbeat ' num2str(user_data.hb_idx) ' : Lead ' user_data.ECG_struct.header.desc(user_data.lead_idx,:)] )
        end

%         zoom reset
        
        user_data.KeyPressed = '';
        
        set(user_data.ECG_hdl,'ButtonDownFcn',@inspect_ECG);            
        set(user_data.ECG_axes_hdl,'ButtonDownFcn',@inspect_ECG);            

        set(fig_hnd, 'UserData', user_data);
        
    end
        
    
function user_data = SearchPattern(user_data)

    point1 = get(user_data.ECG_axes_hdl,'CurrentPoint');    % button down detected
    rbbox;
    point2 = get(user_data.ECG_axes_hdl,'CurrentPoint');    % button up detected            

    xlims = round(sort([ point1(1,1) point2(1,1) ]));
    aux_seq = xlims(1):xlims(2);

    if( ishandle(user_data.Pattern_hdl) )
        delete(user_data.Pattern_hdl)
    end
    disp( 'Filtering ECG ...')
    if( isempty(user_data.filtro) )
        ECG = user_data.ECG_struct.signal(:,user_data.lead_idx);
    else
        ECG = filter(user_data.filtro, flipud(user_data.ECG_struct.signal(:,user_data.lead_idx)) );
        ECG = filter(user_data.filtro, flipud(ECG) );
    end
    
    hold(user_data.ECG_axes_hdl, 'on')
    user_data.Pattern_hdl = plot(user_data.ECG_axes_hdl, aux_seq, ECG(aux_seq,:), '.g' );
    hold(user_data.ECG_axes_hdl, 'off')        
    drawnow;

    pattern2detect = user_data.ECG_struct.signal(aux_seq,user_data.lead_idx);
    pattern2detect = bsxfun( @minus, pattern2detect, mean(pattern2detect));

    [nsamp_pattern nsig_pattern] = size(pattern2detect);
    max_idx = max_index(pattern2detect);
    
%     cross_corr = 1/sqrt(var(ECG)*var(pattern2detect) )*conv(ECG, flipud(pattern2detect) , 'valid');
%     pattern2detect_var = var(pattern2detect);
    disp( 'Looking for the pattern ...')
%     similarity = arrayfun( @(a)( 1/sqrt(var(ECG(a-nsamp_pattern+1:a, :))* pattern2detect_var)*(rowvec(ECG(a-nsamp_pattern+1:a, :)-mean(ECG(a-nsamp_pattern+1:a, :))) * colvec(pattern2detect) ) ), colvec(nsamp_pattern:user_data.ECG_struct.header.nsamp));    

%     similarity = arrayfun( @(a)( 1/sqrt(var(ECG(a-nsamp_pattern+1:a))* pattern2detect_var)*(rowvec(ECG(a-nsamp_pattern+1:a)-mean(ECG(a-nsamp_pattern+1:a))) * colvec(pattern2detect) ) ), colvec(nsamp_pattern:user_data.ECG_struct.header.nsamp));    
%     similarity = [repmat(similarity(1,:),nsamp_pattern-1,1); similarity];

    similarity = cellfun( @(a,b)( diff(conv( a, b, 'same' )) ), mat2cell(user_data.ECG_struct.signal(:,user_data.lead_idx), user_data.ECG_struct.header.nsamp, ones(1,length(user_data.lead_idx))), mat2cell(pattern2detect, diff(xlims)+1, ones(1,length(user_data.lead_idx))), 'UniformOutput', false);
    similarity = cell2mat(cellfun( @(a,b)( diff(conv( a, b, 'same' )) ), similarity, mat2cell(flipud(pattern2detect), diff(xlims)+1, ones(1,length(user_data.lead_idx))), 'UniformOutput', false));
%     similarity = [repmat(similarity(1,:),round((nsamp_pattern-1)/2),1); similarity];
    similarity = abs(mean(similarity,2));
    
    figure(2)
    aux_idx = 1:2*user_data.ECG_struct.header.freq;
    aux_windows = 0:round(user_data.ECG_struct.header.nsamp/4):user_data.ECG_struct.header.nsamp*4/5;
    aux_idx = repmat(aux_idx,length(aux_windows),1);
    aux_idx = bsxfun(@plus, aux_idx',aux_windows);
    aux_idx = colvec(aux_idx);
    aux_val = [(ECG(aux_idx)-mean(ECG(aux_idx))) (similarity(aux_idx)-mean(similarity(aux_idx)))];
    aux_val = bsxfun(@times, aux_val, 1./max(abs(aux_val)));
    plot(aux_val)
    
    disp('Select the threshold to use.')
    [~, thr] = ginput(1);
    
    QRSxlims = 1;
    bContinue = true;
    
    detection_threshold = 0.3;
    
    while(bContinue)
        
        user_data.ECG_struct.pattern_match.time = modmax(similarity, QRSxlims, thr, +1, round(detection_threshold*user_data.ECG_struct.header.freq) );

%         aux_idx = 1;
%         while( ~isempty(aux_idx) )
%             aux_idx = find(diff(user_data.ECG_struct.pattern_match.time) < round(0.15 * user_data.ECG_struct.header.freq));
%             aux_idx2 = setdiff(1:length(user_data.ECG_struct.pattern_match.time), [colvec(aux_idx); colvec(aux_idx+1)]);
%             [~, merged_times] = max([similarity(user_data.ECG_struct.pattern_match.time(aux_idx)) similarity(user_data.ECG_struct.pattern_match.time(aux_idx+1))],[],2);
%             user_data.ECG_struct.pattern_match.time = sort([colvec(user_data.ECG_struct.pattern_match.time(aux_idx2)); colvec(user_data.ECG_struct.pattern_match.time(aux_idx( find(merged_times == 1) ))); colvec(user_data.ECG_struct.pattern_match.time(aux_idx( find(merged_times == 2) )+1 )) ]);
%         end

        if( strcmp(user_data.AnnNames(end,1), cellstr('pattern_match') ) )
            user_data.ratios(end) = CalcRRserieQuality(user_data.ECG_struct.signal, user_data.ECG_struct.header,  user_data.all_annotations, user_data.ECG_struct.noise_power, length(user_data.ratios) );
        else
            user_data.AnnNames = [user_data.AnnNames; cellstr('pattern_match') cellstr('time')];
            user_data.ratios = [user_data.ratios; CalcRRserieQuality(user_data.ECG_struct.signal, user_data.ECG_struct.header, user_data.all_annotations, user_data.ECG_struct.noise_power, size(user_data.AnnNames,1))];
        end
        
        user_data.AnnNames_idx = size(user_data.AnnNames,1);
        user_data.anns_under_edition = user_data.ECG_struct.pattern_match.time;

        user_data.hb_idx = 1;
        user_data.selected_hb_idx = [];
        
        user_data = Redraw(user_data);
        
        ocurrences = length(user_data.ECG_struct.pattern_match.time);
        disp(['Threshold: ' num2str(thr) ' - found ' num2str(ocurrences) ' heartbeats with quality ' num2str(user_data.ratios(end)) ]);

        set(user_data.fig_hdl, 'UserData', user_data);
        
        key = input(['[rt] to refine the QRS detection threshold.\n' ...
                     '[rx] to refine the time window to perform QRS detection.\n' ...
                     '[rh] to refine the minimum time between heartbeats.\n' ...
                     'any key to continue.\n' ...
                     ], 's');
                 
        if( strcmp(key, 'rt') )
            figure(2)
            disp('Select the threshold to use.')
            [~, thr] = ginput(1);
        
        elseif( strcmp(key, 'rh') )
            
            key = input( 'Enter the minimum time between heartbeats:\n' , 's');
                 
            detection_threshold = str2double(key);
            
        elseif( strcmp(key, 'rx') )
            fig_hdl = figure(1);
            disp('Click and drag the time interval to perform QRS detection in the RR series interval ...')
            set(fig_hdl, 'CurrentAxes', user_data.RRserie_axes_hdl);
            waitforbuttonpress;
            point1 = get(user_data.RRserie_axes_hdl,'CurrentPoint');    % button down detected
            rbbox;
            point2 = get(user_data.RRserie_axes_hdl,'CurrentPoint');    % button up detected            
            QRSxlims = round(sort([ point1(1,1) point2(1,1) ]));
            QRSxlims = [ QRSxlims(1) QRSxlims(1) + max(diff(QRSxlims), 10 * user_data.ECG_struct.header.freq )];
        else
            bContinue = false;
        end
        
    end
    
    disp('Search pattern finished.')

    cant_anns = size(user_data.AnnNames,1);
    aux_str = repmat( ' - ',cant_anns,1);
    
    set(user_data.annotation_list_control, 'string', [ char(cellstr(num2str((1:cant_anns)'))) aux_str char(user_data.AnnNames(:,1)) repmat( ' (',cant_anns,1) num2str(round(colvec(user_data.ratios * 1000))) aux_str num2str(colvec(user_data.annotations_ranking))  repmat( ')',cant_anns,1)  ] );
    set(user_data.annotation_under_edition_label, 'string', [ 'Annotation under edition: ' char(user_data.AnnNames( user_data.AnnNames_idx ,1)) ' (' num2str(user_data.ratios(user_data.AnnNames_idx)) ')' ])
    set(user_data.annotation_list_control, 'Value', user_data.AnnNames_idx);    
    
    [~, user_data.best_detections_idx] = sort(user_data.ratios, 'descend');
    
    aux_val = 1:cant_anns;
    
    [~, user_data.annotations_ranking] = sort(aux_val(user_data.best_detections_idx));
    
    close(2)
    

function UpdateRRserieZoom(user_data)

    if( isempty(user_data.hb_idx) )
        return
    end
    
    aux_idx = max(1, user_data.hb_idx - user_data.hb_detail_window):min( length(user_data.RRserie), user_data.hb_idx + user_data.hb_detail_window);

    RRserie_zoom_hdl = plot(user_data.RRserie_zoom_axes_hdl, user_data.anns_under_edition(aux_idx) , user_data.RRserie(aux_idx), '-xb' );

    hold(user_data.RRserie_zoom_axes_hdl, 'on')

    if( ~isempty(user_data.RRserie) )
        [~, aux_idx2] = intersect(user_data.selected_hb_idx, aux_idx);
        RRserie_zoom_hdl = [RRserie_zoom_hdl; colvec(plot(user_data.RRserie_zoom_axes_hdl, user_data.anns_under_edition(user_data.selected_hb_idx(aux_idx2)), user_data.RRserie(user_data.selected_hb_idx(aux_idx2)), '-og' ))];
    end

    if( ~isempty(user_data.RRserie) && user_data.hb_idx < length(user_data.RRserie) )
        RRserie_zoom_hdl = [RRserie_zoom_hdl; colvec(plot(user_data.RRserie_zoom_axes_hdl, user_data.anns_under_edition(user_data.hb_idx), user_data.RRserie(user_data.hb_idx), '-or' ))];
    end
    
    hold(user_data.RRserie_zoom_axes_hdl, 'off')

    xlabel(user_data.RRserie_zoom_axes_hdl, 'Time')
    ylabel(user_data.RRserie_zoom_axes_hdl, 'RR interval')
    
    set(user_data.RRserie_zoom_axes_hdl,'ButtonDownFcn',@inspect_RRserie);  
    set(RRserie_zoom_hdl,'ButtonDownFcn',@inspect_RRserie);  
    
    set(user_data.RRserie_zoom_axes_hdl, 'XTick', rowvec(user_data.anns_under_edition(aux_idx)) );
    aux_str = cellstr(num2str(colvec(aux_idx)));
    aux_hb_idx = find(user_data.hb_idx == aux_idx);
    if( ~isempty(user_data.anns_under_edition) && user_data.hb_idx <= length(user_data.anns_under_edition) )
        aux_str{aux_hb_idx} = [aux_str{aux_hb_idx} ' (' Seconds2HMS(user_data.anns_under_edition(user_data.hb_idx)*1/user_data.ECG_struct.header.freq) ')'];
        set(user_data.RRserie_zoom_axes_hdl, 'XTickLabel', char(aux_str) );
    end
    
function DeleteAnnotations(obj,event_obj) 

	fig_hnd = gcbf; user_data = get(fig_hnd,'userdata');

    cant_anns = size(user_data.AnnNames,1);
    
    if( cant_anns == 0 )
        return;
    end
    
    if( strcmpi(questdlg('Are you sure ?', 'Delete annotations', 'No'), 'yes') )
        
        ann_idx = get(user_data.annotation_list_control, 'Value');
        user_data.ECG_struct = rmfield(user_data.ECG_struct, user_data.AnnNames{ann_idx,1});
        
        if( cant_anns == 1)
            % last annotation deleted
            user_data.ECG_struct.Default.time = [];
            user_data.AnnNames = {'Default' 'time'};
            set(user_data.annotation_list_control, 'string', '1 - Default' );
            set(user_data.annotation_under_edition_label, 'string', 'Annotation under edition: Default' )
            user_data.AnnNames_idx = 1;
            user_data.anns_under_edition = [];
            
        elseif( cant_anns > 1)
            aux_idx = 1:cant_anns;
            aux_idx(ann_idx) = [];
            user_data.AnnNames = user_data.AnnNames(aux_idx,:);
            user_data.ratios = user_data.ratios(aux_idx);
            cant_anns = cant_anns -1;
            set(user_data.annotation_list_control, 'Value', 1);
            
            cant_anns = size(user_data.AnnNames,1);
            aux_str = repmat( ' - ',cant_anns,1);

            set(user_data.annotation_list_control, 'string', [ char(cellstr(num2str((1:cant_anns)'))) aux_str char(user_data.AnnNames(:,1)) repmat( ' (',cant_anns,1) num2str(round(colvec(user_data.ratios * 1000))) aux_str num2str(colvec(user_data.annotations_ranking))  repmat( ')',cant_anns,1)  ] );
            set(user_data.annotation_under_edition_label, 'string', [ 'Annotation under edition: ' char(user_data.AnnNames( user_data.AnnNames_idx ,1)) ' (' num2str(user_data.ratios(user_data.AnnNames_idx)) ')' ])
            
            user_data.AnnNames_idx = 1;
            user_data.anns_under_edition = unique(round(colvec( user_data.ECG_struct.(user_data.AnnNames{user_data.AnnNames_idx,1}).(user_data.AnnNames{user_data.AnnNames_idx,2}) )));
            
        end
        
        user_data.undo_buffer_idx = 1;

        user_data.bAnnsEdited = false;

        user_data.selected_hb_idx = [];        
        
        user_data = Redraw(user_data);

        set(user_data.fig_hdl, 'UserData', user_data);    
        
    end
    

function ChangeRecordingSelected(obj,event_obj) 

	fig_hnd = gcbf; user_data = get(fig_hnd,'userdata');

    if (strcmp(get(fig_hnd,'SelectionType'),'open'))
        %Double click        

    else
        %Single click

        rec_selected = get(obj, 'Value');
        
        if( rec_selected ~= user_data.rec_idx )
            % change of annotation

            if(user_data.bRecEdited)

                if( user_data.bLoadECG )
                    disp('Saving data ...')
                    user_data = update_annotations(user_data);
                    ECG_struct = user_data.ECG_struct;
                    save(user_data.rec_path, '-struct', 'ECG_struct');        
                    disp(['Saved ' user_data.rec_path])
                else
                    disp('Returning data ...')
                    ann_output = user_data.anns_under_edition; 
                end
            end

            user_data.rec_idx = get(user_data.recordings_control, 'Value' );

            user_data = DoRecording(user_data);
        
            set(user_data.fig_hdl, 'UserData', user_data);
            
        end
        
    end


function ChangeLeadsSelected(obj,event_obj) 

	fig_hnd = gcbf; user_data = get(fig_hnd,'userdata');

    if (strcmp(get(fig_hnd,'SelectionType'),'open'))
        %Double click        

    else
        %Single click

        leads_selected = get(obj, 'Value');
        
        if(length(leads_selected ) == 1 && ( length(user_data.lead_idx) ~= length(leads_selected) || leads_selected ~= user_data.lead_idx) )
            % change of annotation

            user_data.standarize_ECG_view = false;
            
            user_data.lead_idx = leads_selected;
                
            cla(user_data.ECG_axes_hdl);
            
            if( isfield(user_data, 'annotation_list_control') )
                anns_selected = get(user_data.annotation_list_control, 'Value');
            else
                anns_selected = [];
            end

            if( length(anns_selected) > 1 )

                aux_anns = {user_data.anns_under_edition};
                anns_selected( anns_selected == user_data.AnnNames_idx) = [];

                for ii = anns_selected
                    aux_anns = [aux_anns; ... 
                                {user_data.ECG_struct.(user_data.AnnNames{ii,1}).(user_data.AnnNames{ii,2})}
                                ];
                end

            else
                aux_anns = user_data.anns_under_edition;
            end
            
            user_data.ECG_hdl = plot_ecg_heartbeat(user_data.ECG_struct.signal(:,user_data.lead_idx), aux_anns, user_data.hb_idx , user_data.hb_detail_window , user_data.ECG_struct.header, user_data.filtro, user_data.ECG_axes_hdl);    

            title(user_data.ECG_axes_hdl, ['Heartbeat ' num2str(user_data.hb_idx) ' : Lead ' user_data.ECG_struct.header.desc(user_data.lead_idx,:)] )

            set(user_data.ECG_hdl,'ButtonDownFcn',@inspect_ECG);            
            set(user_data.ECG_axes_hdl,'ButtonDownFcn',@inspect_ECG);   

            set(user_data.fig_hdl, 'UserData', user_data);
            
            
        elseif(length(leads_selected ) > 1)
                        
            user_data.lead_idx = leads_selected;

            cla(user_data.ECG_axes_hdl);
             
            if( isfield(user_data, 'annotation_list_control') )
                anns_selected = get(user_data.annotation_list_control, 'Value');
            else
                anns_selected = [];
            end

            if( length(anns_selected) > 1 )

                aux_anns = {user_data.anns_under_edition};
                anns_selected( anns_selected == user_data.AnnNames_idx) = [];

                for ii = anns_selected
                    aux_anns = [aux_anns; ... 
                                {user_data.ECG_struct.(user_data.AnnNames{ii,1}).(user_data.AnnNames{ii,2})}
                                ];
                end

            else
                aux_anns = user_data.anns_under_edition;
            end
            
            user_data.ECG_hdl = plot_ecg_heartbeat(user_data.ECG_struct.signal(:,user_data.lead_idx), aux_anns, user_data.hb_idx , user_data.hb_detail_window , user_data.ECG_struct.header, user_data.filtro, user_data.ECG_axes_hdl);    

            title(user_data.ECG_axes_hdl, ['Heartbeat ' num2str(user_data.hb_idx) ' : Leads ' rowvec(colvec([repmat(' ', length(user_data.lead_idx), 1) num2str(colvec(user_data.lead_idx))]')) ] )

            set(user_data.ECG_hdl,'ButtonDownFcn',@inspect_ECG);            
            set(user_data.ECG_axes_hdl,'ButtonDownFcn',@inspect_ECG);   
            
            set(user_data.fig_hdl, 'UserData', user_data);
            
        end
        

        
    end

function user_data = update_annotations(user_data)
    
    user_data.ECG_struct.(user_data.AnnNames{user_data.AnnNames_idx,1}).(user_data.AnnNames{user_data.AnnNames_idx,2}) = unique(round(colvec( user_data.anns_under_edition )));

    if( isempty(strfind(user_data.AnnNames{user_data.AnnNames_idx,1}, 'corrected_' )) )
        user_data.ECG_struct.([ 'corrected_' user_data.AnnNames{user_data.AnnNames_idx,1}]) = user_data.ECG_struct.(user_data.AnnNames{user_data.AnnNames_idx,1});
        user_data.ECG_struct = rmfield(user_data.ECG_struct, user_data.AnnNames{user_data.AnnNames_idx,1});

        user_data.AnnNames{ user_data.AnnNames_idx ,1} = [ 'corrected_' user_data.AnnNames{user_data.AnnNames_idx,1}];

        if( isfield(user_data.ECG_struct, 'series_quality' ) ) 
            user_data.ECG_struct.series_quality.AnnNames = user_data.AnnNames;
        end
    end            


function ChangeAnnotationsSelected(obj,event_obj) 
        
	fig_hnd = gcbf; user_data = get(fig_hnd,'userdata');

    if (strcmp(get(fig_hnd,'SelectionType'),'open'))
        %Double click        
        
        answer = char(inputdlg([ 'Enter the new name of the annotation ' char(user_data.AnnNames( user_data.AnnNames_idx ,1)) ], 'Change annotation name', 1, user_data.AnnNames( user_data.AnnNames_idx ,1)) );
        
        if( ~isempty(answer) && ischar(answer) )
            ann_idx = get(user_data.annotation_list_control, 'Value');
            user_data.ECG_struct.(answer) = user_data.ECG_struct.(user_data.AnnNames{ann_idx,1});
            user_data.ECG_struct = rmfield(user_data.ECG_struct, user_data.AnnNames{ann_idx,1});
            
            user_data.AnnNames{ ann_idx ,1} = answer;
            
            if( isfield(user_data.ECG_struct, 'series_quality' ) ) 
                user_data.ECG_struct.series_quality.AnnNames = user_data.AnnNames;
            end
            cant_anns = size(user_data.AnnNames,1);
            aux_str = repmat( ' - ',cant_anns,1);

            set(user_data.annotation_list_control, 'string', [ char(cellstr(num2str((1:cant_anns)'))) aux_str char(user_data.AnnNames(:,1)) repmat( ' (',cant_anns,1) num2str(round(colvec(user_data.ratios * 1000))) aux_str num2str(colvec(user_data.annotations_ranking))  repmat( ')',cant_anns,1)  ] );
            set(user_data.annotation_under_edition_label, 'string', [ 'Annotation under edition: ' char(user_data.AnnNames( user_data.AnnNames_idx ,1)) ' (' num2str(user_data.ratios(user_data.AnnNames_idx)) ')' ])
            
            user_data.bAnnsEdited = true;
            user_data.bRecEdited = true;   
            
            set(user_data.fig_hdl, 'UserData', user_data);
            
        end

    else
        %Single click

        anns_selected = get(obj, 'Value');
        
        if(length(anns_selected ) == 1 && anns_selected ~= user_data.AnnNames_idx )
            % change of annotation
            
            if(user_data.bAnnsEdited)
                user_data = update_annotations(user_data);
            end

            user_data.AnnNames_idx = anns_selected;

%             disp( ['Using ' user_data.AnnNames{user_data.AnnNames_idx,1} ' annotations (' num2str(user_data.ratios(user_data.AnnNames_idx)) ')'] );

            set(user_data.annotation_under_edition_label, 'string', [ 'Annotation under edition: ' char(user_data.AnnNames( user_data.AnnNames_idx ,1)) ' (' num2str(user_data.ratios(user_data.AnnNames_idx)) ')' ])
            
            user_data.undo_buffer_idx = 1;

            user_data.anns_under_edition = unique(round(colvec( user_data.ECG_struct.(user_data.AnnNames{user_data.AnnNames_idx,1}).(user_data.AnnNames{user_data.AnnNames_idx,2}) )));

            user_data.bAnnsEdited = false;

            user_data.selected_hb_idx = [];

            user_data = Redraw(user_data);

            set(user_data.fig_hdl, 'UserData', user_data);
            
        elseif(length(anns_selected ) > 1)
            
            user_data = Redraw(user_data);

            set(user_data.fig_hdl, 'UserData', user_data);
            
        end
        
    end
    
    
    
        