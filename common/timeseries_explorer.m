%% (Internal) GUI for inspecting and correcting timeseries
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
% Birthdate  : 16/08/2017
% Last update: 16/08/2017
% Copyright 2008-2017
% 
function ts_out = timeseries_explorer(ts_in)

    ts_out

    user_data.RRserie = colvec(diff(user_data.anns_under_edition));
    user_data.RRserie = [user_data.RRserie(1); user_data.RRserie] * 1/user_data.ECG_struct.header.freq;

    limits = prctile(user_data.RRserie, [5 95]);
    
    if(limits(1) > 0) 
        limits(1) = limits(1) * 0.9;
    else
        limits(1) = limits(1) * 1.1;
    end
    limits(2) = limits(2) * 1.1;
        
    if( user_data.bLockScatter )
        x_lims_scatter = get(user_data.Scatter_axes_hdl, 'Xlim');
        y_lims_scatter = get(user_data.Scatter_axes_hdl, 'Ylim');
    end
    
    if( user_data.bLockRRserie )
        x_lims_RRserie = get(user_data.RRserie_axes_hdl, 'Xlim');
        y_lims_RRserie = get(user_data.RRserie_axes_hdl, 'Ylim');
    end
    
    user_data.fig_hdl = figure(1);
%     set(user_data.fig_hdl, 'Toolbar', 'none');
    clf

    %% Axis de Poincar�

    user_data.Scatter_axes_hdl = axes('Position', [0.55 0.41 0.355 0.52 ]);

    RRscatter_hdl = scatter(user_data.Scatter_axes_hdl, user_data.RRserie, [ user_data.RRserie(2:end); user_data.RRserie(end) ] );

    if( ~isempty(user_data.selected_hb_idx) )
        hold(user_data.Scatter_axes_hdl, 'on')
        RRnext = [user_data.RRserie(2:end); user_data.RRserie(end)];
        user_data.RRscatter_selection_hdl = scatter(user_data.Scatter_axes_hdl, user_data.RRserie(user_data.selected_hb_idx), RRnext(user_data.selected_hb_idx), 'xg' );
        hold(user_data.Scatter_axes_hdl, 'off')
    end    
    
    if( user_data.bLockScatter )
        xlim(user_data.Scatter_axes_hdl, x_lims_scatter);
        ylim(user_data.Scatter_axes_hdl, y_lims_scatter);
    else
        xlim(user_data.Scatter_axes_hdl, limits);
        ylim(user_data.Scatter_axes_hdl, limits);
        zoom reset
    end
    
    title(user_data.Scatter_axes_hdl, ['Poincar� plot - Registro ' user_data.ECG_struct.header.recname ]);
    xlabel(user_data.Scatter_axes_hdl, 'RR current')
    ylabel(user_data.Scatter_axes_hdl, 'RR next')

    
    %% Axis de RR serie
    user_data.RRserie_axes_hdl = axes('Position', [0.13 0.709 0.375 0.216 ]);
    
    RRserie_hdl = plot(user_data.RRserie_axes_hdl, user_data.RRserie, '-xb' );

    if( ~isempty(user_data.selected_hb_idx) )
        hold(user_data.RRserie_axes_hdl, 'on')
        user_data.RRserie_selection_hdl = plot(user_data.RRserie_axes_hdl, user_data.selected_hb_idx, user_data.RRserie(user_data.selected_hb_idx), 'og' );
        hold(user_data.RRserie_axes_hdl, 'off')
    end
    
    if( user_data.bLockRRserie )
        xlim(user_data.RRserie_axes_hdl, x_lims_RRserie);
        ylim(user_data.RRserie_axes_hdl, y_lims_RRserie);
    else
        ylim(user_data.RRserie_axes_hdl, limits);
        zoom reset
    end
    
    cant_ticks = 5;
    aux_idx = round(linspace(1, length(user_data.RRserie), cant_ticks));
    set(user_data.RRserie_axes_hdl, 'XTick', rowvec(aux_idx) );
    aux_str = strcat(cellstr(num2str(colvec(aux_idx))), cellstr(repmat(' (',cant_ticks,1)), cellstr(Seconds2HMS(user_data.anns_under_edition(aux_idx)*1/user_data.ECG_struct.header.freq)), cellstr(repmat(')',cant_ticks,1)) );    set(user_data.RRserie_axes_hdl, 'XTickLabel', char(aux_str) );    
    set(user_data.RRserie_axes_hdl, 'XTickLabel', char(aux_str) );
    xlabel(user_data.RRserie_axes_hdl, 'heartbeat #')
    ylabel(user_data.RRserie_axes_hdl, 'RR interval')


