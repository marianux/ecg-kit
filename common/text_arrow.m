function text_arrow_handles = text_arrow( x, y, text_str, text_arrow_prop_vals, axes_hdl )

    %% interface
    text_arrow_handles = [];
    k_length = 0.1;
    
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
        x = [(x + k_length*x_range), x];
    else
        x = [(x - k_length*x_range), x];
    end
    
    if( y - y_lims(1) > y_lims(2) - y  )
        y = [(y + k_length*y_range), y];
    else
        y = [(y - k_length*y_range), y];
    end
    
    [xaf,yaf] = ds2nfu(x, y);
    
    if( any(xaf < 0 | xaf > 1) || any(yaf < 0 | yaf > 1)  )
        return
    end
    
%     dbclear if caught error
    text_arrow_handles = annotation('textarrow', xaf, yaf, 'String' , text_str);

    aux_prop_vals = text_arrow_prop_vals;

    set( text_arrow_handles, aux_prop_vals(:,1)', aux_prop_vals(:,2)' );
    
