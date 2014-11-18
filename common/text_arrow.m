function text_arrow_handles = text_arrow( x, y, text_str, text_arrow_prop_vals, axes_hdl )

    %% interface
    k_length = 0.06;
    k_text = 0.08;
    
    if(nargin < 3 || isempty(text_str) || ~ischar(text_str) )
        text_str = '';
    end

    if(nargin < 4 || ~iscell(text_arrow_prop_vals) )
        text_arrow_prop_vals = [];
    end
    
    if(nargin < 5 || isempty(axes_hdl) )
        axes_hdl = gca;
    end
    
    x_lims = get(axes_hdl, 'xlim');
    x_range = abs(diff(x_lims));
    y_lims = get(axes_hdl, 'ylim');
    y_range = abs(diff(y_lims));
    
    if( x - x_lims(1) > x_lims(2) - x  )
        x_text = x - k_text*x_range;
        x = [(x - k_length*x_range), x];
    else
        x_text = x + k_text*x_range;
        x = [(x + k_length*x_range), x];
    end
    
    if( y - y_lims(1) > y_lims(2) - y  )
        y_text = y - k_text*y_range;
        y = [(y - k_length*y_range), y];
    else
        y_text = y + k_text*y_range;
        y = [(y + k_length*y_range), y];
    end
    
%     [xaf,yaf] = ds2nfu(x, y);
%     
%     if( any(xaf < 0 | xaf > 1) || any(yaf < 0 | yaf > 1)  )
%         return
%     end
%     
% %     dbclear if caught error
%     text_arrow_handles = annotation('textarrow', xaf, yaf, 'String' , text_str);
% 
%     aux_prop_vals = text_arrow_prop_vals;
% 
%     set( text_arrow_handles, aux_prop_vals(:,1)', aux_prop_vals(:,2)' );

    hold(axes_hdl, 'on');

    aux_hdl = plot( x, y );

    
    [~, aux_idx] = setdiff(text_arrow_prop_vals(:,1), {'String' 'TextColor'});
    set( aux_hdl, text_arrow_prop_vals(aux_idx,1)', text_arrow_prop_vals(aux_idx,2)' );
    this_colour = get(aux_hdl, 'Color');
    set( aux_hdl, 'MarkerEdgeColor', 1-this_colour );
    
    bAux = strcmpi('String', text_arrow_prop_vals(:,1));
    if(any(bAux))
        text_str = text_arrow_prop_vals{bAux,2};
    end
    
    txthdl = text( x_text, y_text, text_str, 'Color', 1-this_colour, 'BackgroundColor', this_colour );
    
    hold(axes_hdl, 'off');

    text_arrow_handles = [aux_hdl txthdl];
    
