function reportECG(fig_hdl, filename, bVectorize, bAppend, bHideAxes, bWithTitle)

if( nargin < 6)
    bWithTitle = true;
end

if( nargin < 5)
    bHideAxes = true;
end

if( nargin < 4)
    bAppend = true;
end

if( nargin < 3)
    bVectorize = true;
end

fig_hdl = gcf;

set(fig_hdl, 'Visible', 'off');

paperorient = 'landscape';
paperpos = [0 0 29.7 21.0 ];

set(fig_hdl, 'NumberTitle', 'off');
% set(fig_hdl, 'Position', screenpos);
set(fig_hdl, 'PaperUnits', 'centimeters');
set(fig_hdl, 'PaperType', 'A4');
set(fig_hdl, 'PaperOrientation', paperorient);
set(fig_hdl, 'PaperPosition', paperpos);
set(fig_hdl, 'Units', 'centimeters');
set(fig_hdl, 'Position', paperpos);

axes_hdl = gca;

offset = -10;
ytop = 0.78;

if( bWithTitle )
    set(axes_hdl, 'OuterPosition', [ 0.15 0.05  ytop 0.75]);
else
    set(axes_hdl, 'Position', [ 0.15 0.05  ytop 0.75]);
end

if(bHideAxes)
    set(axes_hdl,'Visible','off')
end


if( bVectorize )
    if ( bAppend )
        print('-dpsc2',sprintf('-f%d',gcf),[filename '.eps'], '-append');
    else
        print('-dpsc2',sprintf('-f%d',gcf),[filename '.eps']);
    end
else
    print('-djpeg', '-r300', sprintf('-f%d',gcf),[filename '.jpg']);
end

if(bHideAxes)
    set(axes_hdl,'Visible','on')
end
