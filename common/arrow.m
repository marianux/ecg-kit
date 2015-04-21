%% (Internal) Creates arrows in grapics
%   
%   arrow_handles = arrow( arrow_start, arrow_end, num_of_arrows, arrow_size, arrow_color, axes_hdl, arrow_style )
% 
% Arguments:
% 
%      + arrow_start/end: arrow limits
% 
%      + num_of_arrows: number of arrows in 1 or 2 extremes.
% 
%      + arrow_size: positive scalar > 0.01. default: 1
% 
%      + arrow_color: [R G B]
% 
%      + axes_hdl: an axes handle to the place to plot.
% 
%      + arrow_style: Some string {'none' 'star4' 'plain' 'rectangle'
%      'ellipse' 'diamond' 'vback1' 'rose' 'vback2' '(Default)'
%      'hypocycloid' 'vback3' 'astroid' 'cback1' 'deltoid' 'cback2'
%      'cback3'   
% 
% Output:
% 
%      + arrow_handles: handler to the created handles.
% 
% Example:
% 
% See also addpath
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function arrow_handles = arrow( arrow_start, arrow_end, num_of_arrows, arrow_size, arrow_color, axes_hdl, arrow_style )

    %% interface
    
    if(nargin < 3 || isempty(num_of_arrows) || ~isnumeric(num_of_arrows) )
        num_of_arrows = 1;
    end

    if(nargin < 4 || isempty(arrow_size) || ~isnumeric(arrow_size) )
        arrow_size = 1;
    end
    arrow_size = max(0.01, arrow_size);
    
    if(nargin < 5 || isempty(arrow_color) || ~isnumeric(arrow_color) )
        arrow_color = [0 0 0];
        line_color = [0 0 0];
    else
        if( size(arrow_color,1) > 1 )
            line_color = arrow_color(2,:);
        else
            line_color = arrow_color(1,:);
        end
        arrow_color = arrow_color(1,:);
    end
    
    if(nargin < 6 || isempty(axes_hdl) )
        axes_hdl = gca;
    end
    
    %% original implementation

%     extra_arguments = varargin;

% 
%     arrow_start = colvec(arrow_start);
%     arrow_end = colvec(arrow_end);
%     
%     arrow_handles = [];
%     arrow_angle_max = pi/6;
%     
%     
%     % esto se hace para que la flecha no se distorcione cuando los ejes x e
%     % y difieren mucho en escala.
% 
%     if( ishold(axes_hdl) )
%         bHoldAxes = false;
%     else
%         bHoldAxes = true;
%     end
%     
%     if( bHoldAxes )
%         hold(axes_hdl, 'on');
%     end
%     
%     prev_state = get(axes_hdl, 'XlimMode');
%     set(axes_hdl, 'XlimMode', 'manual');
%     
%     prev_units = get(axes_hdl, 'Units');
%     set(axes_hdl, 'Units', 'pixels');
%     axes_pos = get(axes_hdl, 'Position');
%     set(axes_hdl, 'Units', prev_units);
%     
%     x_lims = get(axes_hdl, 'Xlim');
%     x_range = abs(diff(x_lims));
%     max_arrow_size_x = 10; % px
%     
%     y_lims = get(axes_hdl, 'Ylim');
%     y_range = abs(diff(y_lims));
%     max_arrow_size_y = max_arrow_size_x  * axes_pos(4) / axes_pos(3);
%     
%     template_arrow = [ -2 -2 0; 0.5 -0.5 0 ];
%     % lo adaptamos al aspect ratio del axis.
%     diff_vecs = arrow_end - arrow_start;
%     
%     rotated_arrow = rotateArrow(-atan2(diff_vecs(2),diff_vecs(1)) );
%     
%     % draw arrow lines
%     arrow_handles = [arrow_handles; plot( axes_hdl, [arrow_start(1,:); arrow_end(1,:)] - 0.9*[-diff_vecs(1,:); diff_vecs(1,:)], [arrow_start(2,:); arrow_end(2,:)] - 0.9*[-diff_vecs(2,:); diff_vecs(2,:)], 'Color', line_color, varargin{:})];
% 
%     % draw end arrow
%     arrow_handles = [arrow_handles;  plotArrow(rotated_arrow, arrow_end)];
%     
%     if(num_of_arrows == 2)
%         rotated_arrow = rotateArrow(pi-atan2(diff_vecs(2),diff_vecs(1)) );
%         % draw start arrow
%         arrow_handles = [arrow_handles;  plotArrow(rotated_arrow, arrow_start)];
%     end
%     
%     if( bHoldAxes )
%         hold(axes_hdl, 'off');
%     end
% 
%     set(axes_hdl, 'XlimMode', prev_state);
%     
%     function rot_mat = rotate_matrix(angle)
%         
%         rot_mat = [ cos(angle)  -sin(angle) ; ... 
%                     sin(angle)  cos(angle)  ; ...
%                     ]';
%     end
% 
%     function [v1,v2,v3] = arrow2vertices(rotated_arrow)
%         v1 = rotated_arrow(:,1);
%         v2 = rotated_arrow(:,2);
%         v3 = rotated_arrow(:,3);
%     end
% 
%     function arrow_handles = plotArrow(rotated_arrow, arrow_location)
%         % draw end arrow
%         [v1,v2,v3] = arrow2vertices(rotated_arrow);
%         cant_arrows = size(arrow_location,2);
%         arrow_handles = colvec(arrayfun(@(a)(fill([v1(1,a) v2(1,a) v3(1,a)] + arrow_location(1,a),[v1(2,a) v2(2,a) v3(2,a)]+ arrow_location(2,a), arrow_color, extra_arguments{:})), 1:cant_arrows ));
%     end
% 
%     function arrow_rot = rotateArrow(rot_angle)
%        
%         aux_norm_ratio = sqrt(max_arrow_size_x^2+max_arrow_size_y^2) / sqrt(sum(template_arrow(:,1).^2));
%         arrow_rot = template_arrow * aux_norm_ratio;
%         arrow_rot = rotate_matrix(rot_angle) * arrow_rot;
%         arrow_rot = [arrow_rot(1,:) * x_range / axes_pos(3); arrow_rot(2,:) * y_range / axes_pos(4)];
%         arrow_rot = arrow_size * arrow_rot;
%         
%     end

    %% Annotation based implementation
    
    arrow_handles = [];
    
    arrow_styles = {'none' 'star4' 'plain' 'rectangle' 'ellipse' 'diamond' 'vback1' 'rose' 'vback2' '(Default)' 'hypocycloid' 'vback3' 'astroid' 'cback1' 'deltoid' 'cback2' 'cback3'};

    if(nargin < 7 || isempty(arrow_style) || any(strcmpi(arrow_styles, arrow_style)) )
        arrow_style = 'vback2';
    end
    
    if(num_of_arrows == 2)
        arrow_str = 'doublearrow';
    elseif(num_of_arrows == 1)
        arrow_str = 'arrow';
    end

    [xaf,yaf] = ds2nfu( [arrow_start(1) arrow_end(1)], [arrow_start(2) arrow_end(2)]);

    if( any(xaf < 0 | xaf > 1) || any(yaf < 0 | yaf > 1)  )
        return
    end
    
    arrow_handles = annotation(arrow_str, xaf, yaf );
    
%     harrow_handle = handle(arrow_handles);
%     try
%         harrow_handle.pinAtAffordance(1);
%         harrow_handle.pinAtAffordance(2);
%         harrow_handle.Pin(1).DataPosition = [arrow_start(1), arrow_end(1), 0];
%         harrow_handle.Pin(2).DataPosition = [arrow_start(2), arrow_end(2), 0];
%     catch
%         % never mind - ignore (no error)
%     end
    
    
    size_props = {'Head1Width' 'Head2Width' 'Head1Length' 'Head2Length' 'LineWidth'};
    aux_size = get( arrow_handles, size_props);
    set( arrow_handles,  size_props, num2cell(arrow_size * cell2mat(aux_size)) );

    set( arrow_handles,  'Color', arrow_color );
    
    set( arrow_handles,  {'Head1Style' 'Head2Style'}, {arrow_style arrow_style} );
    
end