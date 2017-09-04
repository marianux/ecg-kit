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
function ts_out = timeseries_explorer(tsc_in)

    tsc_out

    
    %% Global variables
    
    fig_hdl = figure();
    
    title_efimero_hdl = [];

    
    
    %%
    
    set(fig_hdl, 'WindowButtonDownFcn', @WindowButtonDownCallback);            
    set(fig_hdl, 'WindowButtonUpFcn',   @WindowButtonUpCallback);            
    set(fig_hdl, 'WindowKeyPressFcn',   @WindowKeyPress);            
    
    maximize(fig_hdl);
    maximized_size = get(fig_hdl, 'Position');

    set(fig_hdl, 'Position', [ maximized_size(3:4) maximized_size(3:4) ] .* [ 0.05 0.13 0.95 0.9] );
    
    set(fig_hdl,'CloseRequestFcn',@my_closefcn)

    
    
    
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

    
%% global functions    

    function update_title_efimero( strTitle, delay )
       
        if( ~isempty(title_efimero_hdl) && ~isempty(ancestor(title_efimero_hdl,'figure')) && gcf == ancestor(title_efimero_hdl,'figure') && ishandle(title_efimero_hdl) )
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
        
    end

%% Window callbacks

    function WindowButtonDownCallback(obj,event_obj)
        
        if (strcmp(get(fig_hdl,'SelectionType'),'alt'))

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

    function WindowButtonUpCallback(obj,event_obj)

        DragMouseEnd()
            
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

            point1 = get(ECG_axes_hdl,'CurrentPoint');    % button down detected
            rbbox;
            point2 = get(ECG_axes_hdl,'CurrentPoint');    % button up detected            

            pattern_match_xlims = round(sort([ point1(1,1) point2(1,1) ]));
            aux_seq = pattern_match_xlims(1):pattern_match_xlims(2);

            update_title_efimero('Filtering ECG ...', 5 );

            if( isempty(filtro) )
                ECG = ECG_struct.signal(:,lead_idx);
            else
                ECG = filter(filtro, flipud(ECG_struct.signal(:,lead_idx)) );
                ECG = filter(filtro, flipud(ECG) );
            end

            if( ishandle(Pattern_hdl) )
                delete(Pattern_hdl)
            end

            hold(ECG_axes_hdl, 'on')
            Pattern_hdl = plot(ECG_axes_hdl, aux_seq, ECG(aux_seq,:), '.g' );
            hold(ECG_axes_hdl, 'off')        
            drawnow;            
            
            try
                SearchPattern();
                set(fig_hdl, 'WindowButtonDownFcn', @WindowButtonDownCallback2D);            
                set(ECG_axes_hdl,'ButtonDownFcn',@inspect_ECG);            
                
            catch ME
                
                set(fig_hdl, 'WindowButtonDownFcn', @WindowButtonDownCallback2D);            
                set(ECG_axes_hdl,'ButtonDownFcn',@inspect_ECG);            
                rethrow(ME);
            end
            
            % store the axes_hdl created in SearchPattern();
            axes_hdl( RR_global_PM_k, 1:2) = [RR_global_PM_hdl, figPatternMatch_hdl]; 
            axes_hdl( proximity_k, 1:2) = [proximity_hdl, figPatternMatch_hdl]; 
            axes_hdl( similarity_hdl_k, 1:2) = [similarity_hdl, figPatternMatch_hdl]; 
            axes_hdl( pattMatch_hdl_k, 1:2) = [pattMatch_hdl, figPatternMatch_hdl]; 
            

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
    


end


