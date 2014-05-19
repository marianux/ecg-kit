function [ ECG_hdl axes_hdl fig_hdl anns_hdl all_yranges ] = plot_ecg_mosaic( ECG, varargin )

%%
% Description: 
% This function plots several subplots in the same figure in order to do a
% mosaic with the different leads available in ECG. Annotations can be
% provided to each or for all the mosaic.
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
%     +WinSize: [numeric] OPTIONAL. Default values enclosed in ()
%           Width of the window around each fiducial point provided in
%           QRS_locations. (empty)
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
%     +MaxECGrange: [numeric or string] OPTIONAL. Force a vertial range in order to
%                     ease visual comparison of signals in the mosaic. 
%                   [string] 
%                     'max': force the maximum range to be the range
%                            for all mosaics.
%                     'min', 'mean', 'median': are also available options.
%                     'none': Each mosaic with a different range. (Default).
%     
%     +RowsCols: [numeric] OPTIONAL. Number of rows and columns of the
%                 mosaic. If ommited or if rows * cols ~= ECG_header.nsig, these values are
%                 automatically adapted to the best fit mosaic in relation to the
%                 aspect ratio of the screen.   
%     
%     +FigureHdl: [figure handle] OPTIONAL. Choose the figure to be produced
%                the mosaic. (handle produced by gcf)
% 
%     +ECG_delineation: [struct] OPTIONAL. Default values enclosed in ()
%               Annotation struct generated with wavedet. 
% 
%     +ECG_annotations: [cell] OPTIONAL. Default values enclosed in ()
%               Annotations to be included in the mosaic. The funcion
%               accepts 2 type of annotations: points and lines. 
%             
% 
% Limits and Known bugs:
%   Probably a lot :( ... but dont panic! send me feedback if you need help.
% 
% Example:
% 
% This example makes a synchronized plot of several random events (+1 0 -1)
% ocurryng randomly with additive white gaussian noise.  
% win_size = 100; 
% sig_samp = 10000;
% sig_size = 12;
% event_size = 50;
% x = 0.1*randn(sig_samp,sig_size); 
% event_locations = randsample(win_size:sig_samp-win_size, event_size);
% x(event_locations-1,:) = x(event_locations-1,:) + 1;
% x(event_locations+1,:) = x(event_locations+1,:) - 1;
% x_packed = pack_signal(x, event_locations, win_size);    
% 
% figure(1)
% % estimation of the signal averaged event
% plot_ecg_mosaic( mean(x_packed,3) );
% figure(2)
% % visualization of all events. In this case previous pack_signal call is
% % not needed. 
% plot_ecg_mosaic(x, 'QRS_locations', event_locations, 'WinSize', win_size);
% 
% figure(3)
% % introducing several kind of marks to the plot
% 
% h_line = cell(sig_size,7);
% h2_line = cell(sig_size,7);
% v_line = cell(sig_size,7);
% v2_line = cell(sig_size,7);
% point = cell(sig_size,7);
% a_line = cell(sig_size,7);
% 
% h_line(:,1) = {'line'};
% h_line(:,2) = { [ { 'String'                    'HeadWidth' 'HeadLength' 'LineStyle' 'LineWidth' 'Color' 'TextColor' }; ...
%                   { 'horizontal line text'      5           5            '--'         1.5         'r'     'r'        } ]'}; 
% h_line(1:sig_size, [6 7]) = num2cell( repmat(-0.5,sig_size,2) );
% 
% h2_line(:,1) = {'line'};
% h2_line(:,2) = { [ { 'String'           'HeadWidth' 'HeadLength' 'LineStyle' 'LineWidth' 'Color' 'TextColor' }; ...
%                   { 'other h-line'      5           5            '--'         1.5         'm'     'm'        } ]'}; 
% h2_line(1:sig_size, 4:7) = num2cell( [ repmat(60,sig_size,1) repmat(70,sig_size,1) repmat(0.5,sig_size,2) ] );
% 
% v_line(:,1) = {'line'};
% v_line(:,2) = { [ { 'String'             'HeadWidth' 'HeadLength' 'LineStyle' 'LineWidth' 'Color' 'TextColor' }; ...
%                   { 'vertical line text' 5           5            '--'         1.5         'g'     'g'        } ]'}; 
% v_line(1:sig_size, [4 5]) = num2cell( repmat(20,sig_size,2) );
% 
% v2_line(:,1) = {'line'};
% v2_line(:,2) = { [ { 'String'      'HeadWidth' 'HeadLength' 'LineStyle' 'LineWidth' 'Color' 'TextColor' }; ...
%                   { 'other v-line' 5           5            '--'         1.5         'b'     'b'        } ]'}; 
% v2_line(1:sig_size, 4:7) = num2cell( [ repmat(80,sig_size,2) repmat(-0.8,sig_size,1) repmat(-0.3,sig_size,1) ] );
% 
% point(:,1) = {'point'};
% point(:,2) = { [ { 'String'    'HeadWidth' 'HeadLength' 'Color'       'TextColor'   }; ...
%                  { 'one-point' 5            5           [0.2 0.3 0.4] [0.2 0.3 0.4] } ]'}; 
% point(1:sig_size,4) = num2cell( repmat(50,sig_size,1) );
% 
% a_line(:,1) = {'line'};                    
% a_line(:,2) = { [ { 'String'    'HeadWidth' 'HeadLength' 'LineStyle' 'LineWidth' 'Color' 'TextColor' }; ...
%                   { 'line text' 5           5            '--'         1.5         'k'     'k'        } ]'}; 
% a_line(1:sig_size,4:7) = num2cell( [ repmat(30,sig_size,1) repmat(40,sig_size,1) repmat(0.5,sig_size,1) repmat(-0.5,sig_size,1) ] );
%        
% aux_anns = cat(3,h_line,v_line,h2_line,v2_line,point,a_line);
% 
% plot_ecg_mosaic(mean(x_packed,3), 'ECG_annotations', aux_anns );
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 16/2/2012

%% constants

cModesMaxYrange = {'max','min','mean','median','none'};

% bgTextColor = [251 248 230]/255;
ColorMinorGrid = [0.8 0.8 0.8];

ann_typ_idx = 1;
plot_options_str_idx = 2;
x1_idx = 4;
x2_idx = 5;
y1_idx = 6;
y2_idx = 7;
n_cant_idx = 7;

lmarker_colors = 12;
marker_colors = my_colormap( lmarker_colors );

marker_types = '+o*xsd^v<>ph';
% marker_colors = 'rgbcmykw';
k_magnify = 1.3;

PwaveColor =   [210 255 189]/255;
QRScplxColor = [255 200 255]/255;
TwaveColor =   [255 240 170]/255;

%% parse input

p = inputParser;   % Create instance of inputParser class.
p.addParamValue('ECG_header', [], @(x)(isstruct(x)) );
p.addParamValue('QRS_locations', [], @(x)( isnumeric(x) && all(x > 0) ) );
p.addParamValue('WinSize', [], @(x)( all(isnumeric(x)) ) );
p.addParamValue('ECG_annotations', [], @(x)( isnumeric(x) || iscell(x) ) );
p.addParamValue('ECG_delineation', [], @(x)(isstruct(x)) );
p.addParamValue('MaxECGrange', [], @(x)(isnumeric(x) || ischar(x) && any(strcmpi(cModesMaxYrange, x)) ) );
p.addParamValue('RowsCols', [], @(x)( all(isnumeric(x)) && all(x > 0) ) );
p.addParamValue('FigureHdl', [], @(x)( ishandle(x) ) );

try
    p.parse( varargin{:} );
catch MyError
    fprintf(2, 'Incorrect argval/argvalue, use:\n' );
    fprintf(2, 'plot_ecg_mosaic( ECG, ''arg_name'', ''arg_value'');\n\n' );
    fprintf(2, 'Valid arg_name are:\n' );
    fprintf(2, '+ ECG_header\n' );
    fprintf(2, '+ QRS_locations\n' );
    fprintf(2, '+ ECG_annotations\n' );
    fprintf(2, '+ ECG_delineation\n' );
    fprintf(2, '+ WinSize\n' );
    fprintf(2, '+ MaxECGrange\n' );
    fprintf(2, '+ RowsCols\n' );
    fprintf(2, '+ FigureHdl\n' );
    fprintf(2, '\nOr run ''doc plot_ecg_mosaic'' for details.\n\n' );
    rethrow(MyError);    
end

QRS_locations = p.Results.QRS_locations;
winSize = p.Results.WinSize;
ECG_header = p.Results.ECG_header;
annotations = p.Results.ECG_annotations;
positions = p.Results.ECG_delineation;
all_yranges = p.Results.MaxECGrange;
rowscols = p.Results.RowsCols;
fig_hdl = p.Results.FigureHdl;

delete(p)

if( isempty(fig_hdl) )
    fig_hdl = gcf;
else
    figure(fig_hdl)
end
clf('reset')

if( iscell(ECG) )
    
    [lECG, cant_sig ] = size(ECG{1});
    cant_realizations = length(ECG);

    ECG = cell2mat(reshape(ECG, 1, 1, cant_realizations));
    
elseif( isnumeric(ECG) )

    if( ~isempty(QRS_locations) )
        ECG = pack_signal(ECG, QRS_locations, winSize);    
    end
    
    [lECG, cant_sig, cant_realizations] = size(ECG);

end

if( isempty(all_yranges) )
    bYrange_given = false;
    all_yranges = nan(cant_sig,2);
    rangeMode = 'none';
else
    if(ischar(all_yranges))
        rangeMode = lower(all_yranges);
        all_yranges = nan(cant_sig,2);
        bYrange_given = false;
    else
        bYrange_given = true;
    end
end

% plot arrows when range was not given.
bArrow = ~bYrange_given;

if( isempty(ECG_header) )
    ECG_header.freq = nan;
    ECG_header.nsamp = lECG;
    ECG_header.nsig = cant_sig;
    ECG_header.desc = num2str(colvec(1:cant_sig));
    ECG_header.units = repmat('#', ECG_header.nsig , 1);
    ECG_header.gain = repmat(1, ECG_header.nsig , 1);
    ECG_header.adczero = repmat(0, ECG_header.nsig , 1);
    
    bGridScale = false;
else
    if( ECG_header.nsig ~= cant_sig )
        ECG_header.nsig = cant_sig;
        ECG_header.units = repmat(ECG_header.units(1,:), ECG_header.nsig , 1);
        ECG_header.gain = repmat(ECG_header.gain(1), ECG_header.nsig , 1);
        ECG_header.adczero = repmat(ECG_header.adczero(1), ECG_header.nsig , 1);
    end
    bGridScale = true;
    
    [ECG ECG_header] = ADC2units(ECG, ECG_header, 'uv');
end

if( ~isempty(positions) )

    [min_vals_pw_x, min_vals_pw_y, max_vals_pw_x, max_vals_pw_y]     = CalcMinMax( positions, {'Pon' 'P' 'Poff'} );
    [min_vals_qrs_x, min_vals_qrs_y, max_vals_qrs_x, max_vals_qrs_y] = CalcMinMax( positions, {'QRSon' 'qrs' 'QRSoff'} );
    [min_vals_tw_x, min_vals_tw_y, max_vals_tw_x, max_vals_tw_y]     = CalcMinMax( positions, {'Ton' 'T' 'Toff'} );

    Pwave = cell(cant_sig,7);
    QRScplx = cell(cant_sig,7);
    Twave = cell(cant_sig,7);

    Pwave(:,1) = {'patch'};
    Pwave(:,2) = { [ { 'FaceColor' }; ...   
                     { PwaveColor  } ]'}; 
    Pwave(:,[x1_idx y1_idx x2_idx y2_idx]) = [ num2cell(min_vals_pw_x), num2cell(min_vals_pw_y), num2cell(max_vals_pw_x), num2cell(max_vals_pw_y)];

    QRScplx(:,1) = {'patch'};
    QRScplx(:,2) = { [ { 'FaceColor' }; ...   
                     { QRScplxColor  } ]'}; 
    QRScplx(:,[x1_idx y1_idx x2_idx y2_idx]) = [ num2cell(min_vals_qrs_x), num2cell(min_vals_qrs_y), num2cell(max_vals_qrs_x), num2cell(max_vals_qrs_y)];

    Twave(:,1) = {'patch'};
    Twave(:,2) = { [ { 'FaceColor' }; ...   
                     { TwaveColor  } ]'}; 
    Twave(:,[x1_idx y1_idx x2_idx y2_idx]) = [ num2cell(min_vals_tw_x), num2cell(min_vals_tw_y), num2cell(max_vals_tw_x), num2cell(max_vals_tw_y)];

    annotations = cat(3, Pwave, QRScplx, Twave);
    
end

if( isempty(annotations) )
    
    annotations = cell(cant_sig,1);
    anns_provided = false;
    
else
    if(iscell(annotations))
        
        anns_provided = true;

        [lannotations, cant_params, cant_type_anns] = size(annotations);
        
        if(cant_params ~= n_cant_idx )
            error( 'Invalid number of params: \n\nCol1: Ann point location\nCol2: Ann line start\nCol3: Ann line end\nCol4: plot LineSpec or Prop Name/Val\nCol5: Ann legend\nCol6: Annotation type\n' );
        end
        
        if(lannotations == 1)
            %same ann for all signals
            annotations = repmat(annotations, cant_sig, 1);
        elseif( lannotations > 1 && lannotations < cant_sig )
            error([ 'Invalid number of annotations: ' num2str(lannotations) ' for ' num2str(cant_sig) ' signals.'  ])
        end
        
    elseif(isnumeric(annotations))

        anns_provided = true;
        
        [lannotations, cant_type_anns, cant_anns] = size(annotations);
        
        if( lannotations > 1 && lannotations < cant_sig )
            error([ 'Invalid number of annotations: ' num2str(lannotations) ' for ' num2str(cant_sig) ' signals.'  ])
        end
        
        lmarker_types = length(marker_types);
        aux_ind_colors = repmat(1:lmarker_colors, lmarker_types, 1);
        aux_ind_types = repmat(1:lmarker_types, lmarker_colors, 1);
        aux_ind_coltype = [ colvec(aux_ind_colors') colvec(aux_ind_types') ];

        ann_aux = cell(cant_sig,3,cant_type_anns);
        cell_anntypes = cell(cant_type_anns,2);
        
        for ii = 1:cant_type_anns
            cell_anntypes{ii,1} = [ ':' marker_colors(aux_ind_coltype(ii,1),:) marker_types(aux_ind_coltype(ii,2))];
            cell_anntypes{ii,2} = num2str(ii);
        end
        
        if(lannotations == 1)
            
            for ii = 1:cant_type_anns
                ann_aux(1,:,ii) = { squeeze(annotations(1, ii, :)) cell_anntypes{ii,1} cell_anntypes{ii,2} };
            end
            annotations = repmat(ann_aux(1,:,:), cant_sig, 1);
            
        else
            for jj = 1:cant_sig
                for ii = 1:cant_type_anns
                    ann_aux(jj,:,ii) = { squeeze(annotations(jj, ii, :)) cell_anntypes{ii,1} cell_anntypes{ii,2} };
                end
            end
            annotations = ann_aux;
        end
        
    else
        error('Invalid annotation format.')
    end
end

if( isempty(rowscols) || prod(rowscols) < cant_sig )
    
    prev_units = get(fig_hdl,'units');
    prev_size = get(fig_hdl,'outerposition');
    
    set(fig_hdl,'Visible', 'off');
    set(fig_hdl,'units','normalized');
    set(fig_hdl,'outerposition',[0 0.03 1 0.97]);
    set(fig_hdl,'units','pixels');
    screen_size = get(fig_hdl,'outerposition');
    
    % restore originals
    set(fig_hdl,'outerposition',prev_size);
    set(fig_hdl,'units',prev_units);
    set(fig_hdl,'Visible', 'on');
    
    target_aspect = screen_size(3)/screen_size(4);
    aux_val = [1:cant_sig; ceil(cant_sig./(1:cant_sig)) ] ;
    bExact = (cant_sig - prod(aux_val)) == 0;
    if( any(bExact(2:end-1)) )
        % non-trivial solution
        aux_val = aux_val(:, bExact);
    end
    aux_aspects = aux_val(1,:)./aux_val(2,:);
    [~, aux_idx] = sort(abs(aux_aspects - target_aspect));
    n_rows = aux_val(1,aux_idx(1));
    n_cols = aux_val(2,aux_idx(1));
else
    n_rows = rowscols(1);
    n_cols = rowscols(2);
end

%%

db_status = dbstatus();
bRestoreErrorStatus = false;
if( length(db_status) > 1 && strcmpi(db_status(end).cond, 'caught error')  )
    dbclear if caught error
    bRestoreErrorStatus = true;
end

ECG_hdl = [];
axes_hdl = nan(cant_sig,1);
anns_hdl = [];

if( isnan(ECG_header.freq) )
    xTimeGridSpacing = round(lECG/5); % mseg
    xTimeGrid = 0:xTimeGridSpacing:lECG;
else
    xTimeGridSpacing = 50; % mseg
    xTimeGrid = 0:(xTimeGridSpacing/1e3*ECG_header.freq):lECG;
end
    

plotXmin = 1;
plotXmax = lECG;
plotXrange = plotXmax - plotXmin;
xTextOffset = 0.005*plotXrange;

if(anns_provided)
    lineWidth = 1.3;
else
    lineWidth = 1;
end

v_space = 0.02/(n_rows+1);
h_space = 0.02/(n_cols+2);
axe_width = 0.98/n_rows;
axe_height = 0.98/n_cols;
axe_x = h_space;
axe_y = 1-axe_height-v_space;

for ii = 1:cant_sig

    txthdl = [];

    axes_hdl(ii) = axes();
    
    % AUX ECG plot
    
    % ECG plot
    ECG_hdl = [ ECG_hdl; plot( squeeze(ECG(:,ii,:)), 'LineWidth', lineWidth )];
    
    if(~bYrange_given)
        all_yranges(ii,:) = rowvec(get(axes_hdl(ii), 'ylim'));
    end

    set(axes_hdl(ii), 'xlim', [1 lECG]);
    
    prev_units = get(axes_hdl(ii),'units');
    set(axes_hdl(ii),'units','pixels');
    axis_position = get(axes_hdl(ii),'Position');
    set(axes_hdl(ii),'units',prev_units);
    
    set(axes_hdl(ii), 'Box', 'off' );
    set(axes_hdl(ii), 'Xtick', [] );
    set(axes_hdl(ii), 'Ytick', [] );
    set(axes_hdl(ii), 'Xcolor', [1 1 1] );
    set(axes_hdl(ii), 'Ycolor', [1 1 1] );
    
    [x_loc, y_loc] = ind2sub([n_rows n_cols], ii);
    axe_x = h_space * x_loc + axe_width * (x_loc-1);
    axe_y = 1-( v_space * y_loc + axe_height * y_loc);
    
%     aux_val = get(axes_hdl(ii), 'Position' );
%     set(axes_hdl(ii), 'Position', [ aux_val(1:2) - (k_magnify-1)/2*aux_val(3:4) k_magnify*aux_val(3:4) ] );
    set(axes_hdl(ii), 'Position', [ axe_x axe_y axe_width axe_height ] );
    
end

max_yrange = [min(all_yranges(:,1)) max(all_yranges(:,2)) ];
% min_range = min(all_yranges(:,2) - all_yranges(:,1));
    
for ii = 1:cant_sig
    
    set(fig_hdl, 'CurrentAxes', axes_hdl(ii));
    hold(axes_hdl(ii), 'on');
    
    if(bGridScale) 
        % grid plot

        yTimeGridSpacing = (all_yranges(ii,2) - all_yranges(ii,1)) / 5; % uV
        
        yTimeGrid = max_yrange(1):yTimeGridSpacing:max_yrange(2);
        
        lyTimeGrid = length(yTimeGrid);
        if( lyTimeGrid > 10 )
            yTimeGrid = yTimeGrid( yTimeGrid >= all_yranges(ii,1) & yTimeGrid <= all_yranges(ii,2) );
        end
        
        plot( repmat([1; lECG], 1, length(yTimeGrid)), repmat(yTimeGrid,2,1), 'Color', ColorMinorGrid, 'LineStyle', ':' );
        plot( repmat(xTimeGrid,2,1), repmat(colvec(max_yrange), 1, length(xTimeGrid) ), 'Color', ColorMinorGrid, 'LineStyle', ':' );

    end

    % zero line
    plot( [ xTimeGrid(1) xTimeGrid(end)], [0 0], '-k' );
    
    hold(axes_hdl(ii), 'off');
   
end

if(~bYrange_given)
    switch(rangeMode)
        case 'max'
            all_yranges = repmat([ min(all_yranges(:,1)) max(all_yranges(:,2)) ], cant_sig, 1   );
        case 'min'
            all_yranges = repmat([ max(all_yranges(:,1)) min(all_yranges(:,2)) ], cant_sig, 1   );
        case 'mean'
            all_yranges = repmat([ mean(all_yranges(:,1)) mean(all_yranges(:,2)) ], cant_sig, 1   );
        case 'median'
            all_yranges = repmat([ median(all_yranges(:,1)) median(all_yranges(:,2)) ], cant_sig, 1   );
    end
end


for ii = 1:cant_sig
%     aux_val = get(axes_hdl(ii), 'ylim');
%     aux_val2 = (all_yranges - abs(diff(aux_val)))/2;
%     set(axes_hdl(ii), 'ylim', aux_val + [-aux_val2 aux_val2] );
    set(axes_hdl(ii), 'ylim', all_yranges(ii,:) );
    
    set(fig_hdl, 'CurrentAxes', axes_hdl(ii));
    y_lims = get(axes_hdl(ii), 'ylim');
    x_lims = get(axes_hdl(ii), 'xlim');
    
    if( ii <= size(ECG_header.desc,1) )
        text(x_lims(1) + 0.01*abs(diff(x_lims)), y_lims(2) - 0.05*abs(diff(y_lims)), ECG_header.desc(ii,:), 'EdgeColor', [0 0 0] );
    else
        text(x_lims(1) + 0.01*abs(diff(x_lims)), y_lims(2) - 0.05*abs(diff(y_lims)), ['Sig ' num2str(ii)], 'EdgeColor', [0 0 0] );
    end
    
    plotXmin = x_lims(1);
    plotXmax = x_lims(2);
    plotXrange = plotXmax - plotXmin;
    plotYmin = y_lims(1);
    plotYmax = y_lims(2);
    plotYrange = plotYmax - plotYmin;
    yTextOffset = 0.01*plotYrange;

    yTimeGridSpacing = (all_yranges(ii,2) - all_yranges(ii,1)) / 5; % uV
    yTimeGrid = max_yrange(1):yTimeGridSpacing:max_yrange(2);

    aux_idx = find(yTimeGrid > plotYmin & yTimeGrid < plotYmax);

    width_legend = 0.18*plotXrange;
    height_legend = 0.1*plotYrange;
    left_legend = plotXmax - width_legend - 5*xTextOffset;
    bottom_legend = plotYmax - height_legend - 2*yTextOffset;
    
    % scale X
%         plot( [-xTimeGridSpacing/2 -xTimeGridSpacing/2 xTimeGridSpacing/2 xTimeGridSpacing/2] + 5/7*lECG, plotYmin + 0.05*plotYrange + [yTextOffset 0 0 yTextOffset], 'Color', [0 0 0], 'LineWidth', 1 );
    if( bArrow )
        arrow( [ xTimeGrid(2) plotYmin + 0.05*plotYrange], [ xTimeGrid(3) plotYmin + 0.05*plotYrange], 2, 0.5, [0 0 0], axes_hdl(ii) );
    end
    
    if( isnan(ECG_header.freq) )
        aux_str =  [ num2str(xTimeGridSpacing) ' #' ];
    else
        aux_str =  [ num2str(xTimeGridSpacing) ' ms' ];
    end
    
    text( mean(xTimeGrid(2:3)),  plotYmin + 0.05*plotYrange + 5*yTextOffset, aux_str, 'FontSize', 8, 'HorizontalAlignment', 'center');

    % scale Y
%         plot( [ -xTextOffset 0 0 -xTextOffset] + 0.95*lECG, plotYmin + 5/7*plotYrange + [-yTimeGridSpacing/2 -yTimeGridSpacing/2 yTimeGridSpacing/2 yTimeGridSpacing/2], 'Color', [0 0 0], 'LineWidth', 1 );
    if( bArrow )
        arrow( [0.95*lECG - xTextOffset; yTimeGrid(aux_idx(2))], [0.95*lECG - xTextOffset; yTimeGrid(aux_idx(3))], 2, 0.5, [0 0 0], axes_hdl(ii) );
    end
    
    text( 0.95*lECG - 10*xTextOffset, mean(yTimeGrid(aux_idx(2:3))), [ num2str(yTimeGridSpacing) ' ' ECG_header.units(ii,:) ], 'FontSize', 8, 'HorizontalAlignment', 'center', 'Rotation', 90);
    
    txthdl = [];
    
    if( ~isempty(annotations{ii}) )
        
        x_lims = get(axes_hdl(ii), 'xlim');
        x_range = abs(diff(x_lims));
        y_lims = all_yranges(ii,:);
        y_range = abs(diff(y_lims));
        
        for jj = 1:cant_type_anns
            
            switch( annotations{ii,ann_typ_idx,jj} )
                    
                case 'point'
                    
                    if( ~isnan(annotations{ii,x1_idx,jj}) )

                        txthdl = [ txthdl text_arrow(annotations{ii,x1_idx,jj}, ECG(annotations{ii,x1_idx,jj}, ii), '', annotations{ii,plot_options_str_idx,jj}, axes_hdl(ii) )];
                        
%                         aux_hdl = eval( [ 'plot( annotations{' num2str(ii) ',' num2str(x1_idx) ',' num2str(jj) '}, ECG(annotations{' num2str(ii) ',' num2str(x1_idx) ',' num2str(jj) '},' num2str(ii) ',1)  ' annotations{ii,plot_options_str_idx,jj} ')']);
%                         anns_hdl = [ anns_hdl; colvec(aux_hdl)];
%                         txthdl = [ txthdl; text( annotations{ii,x1_idx,jj} + 0.03*x_range, ECG(annotations{ii,x1_idx,jj}, ii, 1), annotations{ii,text_idx,jj}, 'BackgroundColor', bgTextColor )];
                    end
                    
                case 'line'
                    
                    if( isempty(annotations{ii,x1_idx,jj}) )
                        annotations{ii,x1_idx,jj} = repmat(x_lims(1),1,length(annotations{ii,y1_idx,jj}));
                    end
                    if( isempty(annotations{ii,x2_idx,jj}) )
                        annotations{ii,x2_idx,jj} = repmat(x_lims(2),1,length(annotations{ii,y1_idx,jj}));
                    end
                    
                    if(cant_realizations > 1)
                        if( isempty(annotations{ii,y1_idx,jj}) )
                            annotations{ii,y1_idx,jj} = min(squeeze(ECG(annotations{ii,x1_idx,jj},ii,:)),[],2);
                        end
                        if( isempty(annotations{ii,y2_idx,jj}) )
                            annotations{ii,y2_idx,jj} = max(squeeze(ECG(annotations{ii,x1_idx,jj},ii,:)),[],2);
                        end
                    else
                        if( isempty(annotations{ii,y1_idx,jj}) )
                            annotations{ii,y1_idx,jj} = repmat(min(ECG(:,ii)),1,length(annotations{ii,x1_idx,jj}));
                        end
                        if( isempty(annotations{ii,y2_idx,jj}) )
                            annotations{ii,y2_idx,jj} = repmat(max(ECG(:,ii)),1,length(annotations{ii,x1_idx,jj}));
                        end
                    end
                    
                    if( ~isempty(cell2mat(annotations(ii,4:7,jj)))   )
                        
                        txthdl = [ txthdl text_line(cell2mat(annotations(ii,[x1_idx x2_idx],jj)), cell2mat(annotations(ii,[y1_idx y2_idx],jj)), '', annotations{ii,plot_options_str_idx,jj}, axes_hdl(ii) )];
                        
%                         aux_hdl = eval( [ 'plot( cell2mat(annotations(' num2str(ii) ', [' num2str(x1_idx) ' ' num2str(x2_idx) '],' num2str(jj) ')), cell2mat(annotations(' num2str(ii) ', [' num2str(y1_idx) ' ' num2str(y2_idx) '],' num2str(jj) '))'  annotations{ii,plot_options_str_idx,jj} ')']);
%                         anns_hdl = [ anns_hdl; colvec(aux_hdl)];
%                         if( ~strcmpi('', annotations{ii,text_idx,jj}) )
%                             txthdl = [ txthdl; text(annotations{ii,x1_idx,jj} + 0.03*x_range, annotations{ii,y1_idx,jj} + 0.03*y_range, annotations{ii,text_idx,jj}, 'VerticalAlignment', 'bottom','BackgroundColor', bgTextColor )];
% %                             set(txthdl(end), 'Rotation', round(180/pi*atan2( (annotations{ii,y2_idx,jj} - annotations{ii,y1_idx,jj}) * axis_position(4) / y_range , (annotations{ii,x2_idx,jj} - annotations{ii,x1_idx,jj}) * axis_position(3) / x_range )) );
% %                             set(txthdl(end), 'Rotation', 90 );
%                         end                    
                    end
                    
                case 'patch'
                    
                    if( isempty(annotations{ii,x1_idx,jj}) )
                        annotations{ii,x1_idx,jj} = repmat(x_lims(1),1,length(annotations{ii,y1_idx,jj}));
                    end
                    if( isempty(annotations{ii,x2_idx,jj}) )
                        annotations{ii,x2_idx,jj} = repmat(x_lims(2),1,length(annotations{ii,y1_idx,jj}));
                    end
                    
                    if(cant_realizations > 1)
                        if( isempty(annotations{ii,y1_idx,jj}) )
                            annotations{ii,y1_idx,jj} = min(squeeze(ECG(annotations{ii,x1_idx,jj},ii,:)),[],2);
                        end
                        if( isempty(annotations{ii,y2_idx,jj}) )
                            annotations{ii,y2_idx,jj} = max(squeeze(ECG(annotations{ii,x1_idx,jj},ii,:)),[],2);
                        end
                    else
                        if( isempty(annotations{ii,y1_idx,jj}) )
                            annotations{ii,y1_idx,jj} = repmat(min(ECG(:,ii)),1,length(annotations{ii,x1_idx,jj}));
                        end
                        if( isempty(annotations{ii,y2_idx,jj}) )
                            annotations{ii,y2_idx,jj} = repmat(max(ECG(:,ii)),1,length(annotations{ii,x1_idx,jj}));
                        end
                    end                    
                    
                    if( ~isempty(cell2mat(annotations(ii,4:7,jj)))   )
                        
                        aux_hdl = patch( [annotations{ii,x1_idx,jj} annotations{ii,x1_idx,jj} annotations{ii,x2_idx,jj} annotations{ii,x2_idx,jj} ], [ annotations{ii,y1_idx,jj} annotations{ii,y2_idx,jj} annotations{ii,y2_idx,jj} annotations{ii,y1_idx,jj} ], [1 1 1], 'EdgeColor', 'none');
                        aux_prop_vals = annotations{ii,plot_options_str_idx,jj};
                        set( aux_hdl, aux_prop_vals(:,1)', aux_prop_vals(:,2)' );
                        
                        txthdl = [ txthdl; aux_hdl];

                    end                    
                    
            end
            
        end
        
    end
   
    uistack(txthdl,'bottom');

    patch([left_legend left_legend [left_legend left_legend]+width_legend left_legend ], [bottom_legend [bottom_legend bottom_legend]+height_legend bottom_legend bottom_legend], [1 1 1], 'EdgeColor', [0 0 0]);
    text( left_legend + width_legend/4, bottom_legend + height_legend/2, 'P', 'FontSize', 8, 'HorizontalAlignment', 'center', 'BackgroundColor', PwaveColor);
    text( left_legend + width_legend/2, bottom_legend + height_legend/2, 'QRS', 'FontSize', 8, 'HorizontalAlignment', 'center', 'BackgroundColor', QRScplxColor );
    text( left_legend + width_legend*3/4, bottom_legend + height_legend/2, 'T', 'FontSize', 8, 'HorizontalAlignment', 'center', 'BackgroundColor', TwaveColor );
    
end

if(bRestoreErrorStatus)
    dbstop if caught error
end

if( nargout == 0 )
    ECG_hdl = [];
end


    function [min_vals_x, min_vals_y, max_vals_x, max_vals_y] = CalcMinMax( this_annotation, field_names )

        if( isfield(this_annotation, field_names{1} ) )
            aux_on = cell2mat({this_annotation(:).(field_names{1})});
        else
            aux_on = [];
        end

        if( isfield(this_annotation, field_names{2} ) )
            aux_peak = cell2mat({this_annotation(:).(field_names{2})});
        else
            aux_peak = [];
        end

        if( isfield(this_annotation, field_names{3} ) )
            aux_off = cell2mat({this_annotation(:).(field_names{3})});
        else
            aux_off = [];
        end

        start_sample = 1;

        % wave start
        bOn = ~isnan(aux_on) & aux_on >= start_sample & aux_on <= start_sample+lECG;
        % wave end
        bOff = ~isnan(aux_off) & aux_off >= start_sample & aux_off <= start_sample+lECG;

        % wave peak
        bPeak = ~isnan(aux_peak) & aux_peak >= start_sample & aux_peak <= start_sample+lECG;

        % wave conection between start-end
        bOnOff = bOn & bOff;
        aux_complete_idx = find(bOnOff);

        aux_complete_idxx = num2cell(nan(ECG_header.nsig,1));
        % normal sampled version
        aux_complete_idxx(aux_complete_idx) = arrayfun( @(a)( max(1, aux_on(a)):min(lECG,aux_off(a)) ),aux_complete_idx, 'UniformOutput', false);
        if( any(~bOnOff & bOn & bPeak) )
            % on-peak
            aux_complete_idx = find( ~bOnOff & bOn & bPeak);
            aux_complete_idxx(aux_complete_idx) = arrayfun( @(a)( aux_on(a):aux_peak(a) ),aux_complete_idx, 'UniformOutput', false);
        end

        if( any(~bOnOff & bOff & bPeak) )
            % peak-off
            aux_complete_idx = find(~bOnOff & bPeak & bOff);
            aux_complete_idxx(aux_complete_idx) = arrayfun( @(a)( aux_peak(a):aux_off(a) ),aux_complete_idx, 'UniformOutput', false);
        end         

        max_vals_x = cellfun( @(a,lead)( a(end) ), aux_complete_idxx);
        max_vals_y = cellfun( @(a,lead)( max(dummyfun(ECG,a, lead))), aux_complete_idxx, num2cell((1:ECG_header.nsig)'));
        min_vals_x = cellfun( @(a,lead)( a(1) ), aux_complete_idxx);
        min_vals_y = cellfun( @(a,lead)( min(dummyfun(ECG,a, lead))), aux_complete_idxx, num2cell((1:ECG_header.nsig)'));

    end

    function fun_val = dummyfun(a,b,c)
        
        fun_val = nan;
        if( ~any(isnan(b)) && ~any(isnan(c)) )
            fun_val = a(b,c);
        end
        
    end


end
