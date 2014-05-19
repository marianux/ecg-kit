function text_line_handles = text_line( x, y, text_str, text_line_prop_vals, axes_hdl )

    %% interface
    k_length = 0.1;
    
    if(nargin < 3 || isempty(text_str) || ~ischar(text_str) )
        text_str = '';
    end

    if(nargin < 4 || ~iscell(text_line_prop_vals) )
        text_line_prop_vals = [];
    end
    
    if(nargin < 5 || isempty(axes_hdl) )
        axes_hdl = gca;
    end

    text_line_handles = text_arrow(x(1), y(1), '', text_line_prop_vals, axes_hdl);

    aux_prop_vals = text_line_prop_vals;

    set( text_line_handles, aux_prop_vals(:,1)', aux_prop_vals(:,2)' );
    
    aux_hdl = arrow( [x(1) y(1)], [x(2) y(2)], 2, 0.5, [0 0 0], axes_hdl );

    bStr = strcmpi('String', aux_prop_vals(:,1)) | strcmpi('TextColor', aux_prop_vals(:,1));
    
    set( aux_hdl, aux_prop_vals(~bStr,1)', aux_prop_vals(~bStr,2)' );
    
    set( aux_hdl, {'Head1Style' 'Head2Style' }, {'diamond' 'diamond'} );
    set( aux_hdl, {'Head1Width' 'Head2Width' 'Head1Length' 'Head2Length' }, {2 2 2 2} );
    
    text_line_handles = [ text_line_handles aux_hdl ];

