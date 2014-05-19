function plot_auc( areas, pvalues, strTitle, strLegend, bPrint, fig_hdl)

if( nargin < 4  ) 
    strLegend = [];
end

if( nargin < 3  ) 
    strTitle = [];
end

if( nargin < 5 || isempty(bPrint) ) 
    bPrint = false;
end

if( size(areas,1) ~= 3 )
    error('areas should have lower CI - mean Area - upper CI in the rows')
end

if( nargin < 6 || isempty(fig_hdl) ) 
    fig_hdl = figure();
else
    fig_hdl = figure(fig_hdl);
    clf(fig_hdl);
end


areas2plot = size(areas,2);
Xvalues = [ ];
Yvalues = [ ];
significance_levels = [0.05 0.01 0.001];
significance_levels_str = {'0.05' '0.01' '0.001'};

box_dist = 1;
box_width = box_dist* 2/3;
box_loc = 1:box_dist:box_dist*areas2plot;
box_height_k = 0.75;
colors = my_colormap(areas2plot);
areas_min = min(min(areas));
areas_max = max(max(areas));
areas_range = areas_max - areas_min;

if(bPrint)
    strFontName = 'LMRoman12';
    FontSizeTitle = 22;
    FontSizeLabel = 20;
    FontSizeTick = 18;
else
    strFontName = 'Helvetica';
    FontSizeTitle = 12;
    FontSizeLabel = 12;
    FontSizeTick = 10;
end

axes_hdl = gca();
cla(axes_hdl);

xlim([0 (box_dist*areas2plot)+box_dist])
ylim([areas_min - 0.1*areas_range  areas_max + 0.2*areas_range ])

hold on

for ii = 1:areas2plot
    
    box_height = box_height_k * (areas(3,ii) - areas(1,ii));
    this_colour = colors( 1+rem(ii-1, size(colors,1)),:);
    plot( box_loc(ii) + [-box_width/4 box_width/4], [areas(3,ii) areas(3,ii)], 'Color', this_colour, 'LineWidth', 2 );
    plot( box_loc(ii) + [-box_width/4 box_width/4], [areas(1,ii) areas(1,ii)], 'Color', this_colour, 'LineWidth', 2 )
    plot( [box_loc(ii) box_loc(ii)], [areas(3,ii) areas(1,ii)], 'Color', this_colour, 'LineWidth', 2 )
    fill( box_loc(ii) + [-box_width/2 -box_width/2 box_width/2 box_width/2], areas(2,ii) + [ -box_height/2 box_height/2 box_height/2 -box_height/2], this_colour, 'EdgeColor',  this_colour );
    text( box_loc(ii), areas(2,ii), num2str(areas(2,ii), '%3.2f'), 'HorizontalAlignment', 'center', 'FontName', strFontName, 'FontSize', FontSizeTick, 'Color', (1-this_colour));
    
end

for ii = 1:(areas2plot-1)
    for jj = ii+1:areas2plot
        if( pvalues(ii,jj) < 0.05 ) 
            y_loc = max(areas(3,ii), areas(3,jj)) + (0.05*areas_range);
            if( bPrint )
                plot( [box_loc(ii) box_loc(ii) box_loc(jj) box_loc(jj)], y_loc + (0.02*areas_range) * [0 1 1 0], 'LineWidth', 1 , 'Color', [0 0 0] );
            else
                plot( [box_loc(ii) box_loc(ii) box_loc(jj) box_loc(jj)], y_loc + (0.02*areas_range) * [0 1 1 0], 'LineWidth', 2 , 'Color', [0 0 0] );
            end
        end    
    end
    
    % overprint the p-value
    for jj = ii+1:areas2plot
        aux_idx = find(pvalues(ii,jj) < significance_levels, 1, 'last');
        if( ~isempty(aux_idx) ) 
            y_loc = max(areas(3,ii), areas(3,jj)) + (0.05*areas_range);

            if( bPrint )
                aux_hdl = text( box_loc(ii)+(box_loc(jj) - box_loc(ii))/2, y_loc + (0.07*areas_range), [ 'p < ' significance_levels_str{aux_idx} ], 'HorizontalAlignment', 'center', 'FontName', strFontName, 'FontSize', FontSizeTick, 'Visible', 'off');
                aux_extent = get(aux_hdl, 'Extent');
                delete(aux_hdl)
                aux_extent(1) = aux_extent(1) + 0.1*aux_extent(3);
                aux_extent(3) = 0.8*aux_extent(3);
                aux_extent(2) = aux_extent(2) + 0.25*aux_extent(4);
                aux_extent(4) = 0.7*aux_extent(4);
                fill( aux_extent(1)+ [ 0 0 aux_extent(3) aux_extent(3) ], aux_extent(2)+ [ 0 aux_extent(4) aux_extent(4) 0 ], [1 1 1], 'EdgeColor', [1 1 1]  );
                text( box_loc(ii)+(box_loc(jj) - box_loc(ii))/2, y_loc + (0.07*areas_range), [ 'p < ' significance_levels_str{aux_idx} ], 'HorizontalAlignment', 'center', 'FontName', strFontName, 'FontSize', FontSizeTick);
            else
                text( box_loc(ii)+(box_loc(jj) - box_loc(ii))/2, y_loc + (0.07*areas_range), [ 'p < ' significance_levels_str{aux_idx} ], 'HorizontalAlignment', 'center', 'FontName', strFontName, 'FontSize', FontSizeTick, 'BackgroundColor', [1 1 1]);
            end
        end    
    end
    
end

hold off

set(axes_hdl, 'Box', 'off' );
set(axes_hdl, 'Xtick', box_loc );
set(axes_hdl, 'XtickLabel', strLegend );
rotateticklabel(axes_hdl, 20);
set(axes_hdl, 'Ytick', linspace(areas_min, areas_max, 5) );
set(axes_hdl, 'YtickLabel', num2str(colvec(linspace(areas_min, areas_max, 5)), '%3.2f') );



set(axes_hdl, 'FontName', strFontName );
set(axes_hdl, 'FontSize', FontSizeTick );

title(strTitle, 'FontName', strFontName, 'FontSize', FontSizeTitle)

ylabel('AUC', 'FontName', strFontName, 'FontSize', FontSizeLabel)
