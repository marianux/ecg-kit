function tutorial_plot2svg
% First, let's look at the temperature in Germany depending on the month
% (averaged 1961-1990, source http://de.wikipedia.org/wiki/Klima)
name={'Januar','Februar','March','April','May','June','July','August','September','Oktober','November','December'};
avgT =[-0.5       0.5      3.7     7.6   12.2   15.5   17.1    16.9      13.8        9.4       4.2        0.9];
minT =[-3.0      -2.5      0.0     3.0    7.3   10.6   12.3    12.0       9.3        5.7       1.6       -1.5];
maxT =[ 2.0       3.4      7.5    12.1   17.2   20.4   22.0    21.9      18.4       13.1       6.9        3.2];
% Draw a new figure and set the default values for the text font
fig = figure;
set(fig,'DefaultAxesFontName','Arial')
set(fig,'DefaultAxesFontSize', 16)
% Let's plot the data
h = plot( 1:12, maxT, 1:12, avgT, 1:12, minT);
% Now we change the XTickLabels
set(gca, 'XTick', 1:12);
set(gca, 'XTickLabel', name);
% We change the LineWidth ...
set(h(2), 'LineWidth', 12);
set(h([1 3]), 'LineWidth', 8);
set(gca, 'LineWidth', 2);
% ... and add a legend
[h1,h2] = legend(h, {'Max T', 'Avg T', 'Min T'});
% Next we add a title, ylabel, grid, box on
title('Temperature in Germany (1961-1990)');
ylabel('T [°C]');
grid on
box on
axis([1 12 -5 25])

% Up to here we have created a standard Matlab figure.
% We create a svg plot as reference
plot2svg('temperature_standard.svg');
% Let's improve first the x-axis. The labels overlap. This doesn't look
% nice. Unfortunately, the tickmark labels do not have the same Matlab
% properties as standard text. Therfore, it is not possible to rotate them.
setting.svg.XTickLabelAngle = -40;
set(gca, 'UserData', setting);
% Just turning the labels won't help, as they would be cut at the bottom.
% Therfore, we shift the axes position up.
set(gca, 'Position', [0.13 0.21 0.720 0.715]);
% Now let's beautify the plot by adding some light and shadow
drop_shadow_lighting(h);
for i = 1:length(h2)
    if ~strcmp(get(h2(i),'Type'),'text')
        drop_shadow_lighting(h2(i));
    end
end
plot2svg('temperature_nicer.svg');
% But, this is not enough. To make it perfect we replace the line color by
% a background pixel image (1px x 50px). This pixel image is scaled to cover
% the whole axes region. 
colors = reshape(flipud(jet(50)),[50 1 3]);
imwrite(colors, 'gradient.png', 'png')
% To add the pixel image the used filter chain is modified.
drop_shadow_lighting_background(h, 'gradient.png');
for i = 1:length(h2)
    if ~strcmp(get(h2(i),'Type'),'text')
        drop_shadow_lighting_background(h2(i), 'gradient.png');
    end
end
plot2svg('temperature_perfect.svg');
% Filters are a very powerful feature of SVG. Unfortunately, not all SVG
% render tools give full support. The best browser for SVG filters at the
% moment seems to be Opera directly followed by Firefox. For Inkscape I had
% to add a workaround in the plot2svg as it does not correctly handle
% feImage filters. IE + RENESIS will not work as it does not support
% filters. IE + Adobe may work, but is not tested. The Safari browser may
% also fully support filters (not tested).



function drop_shadow_lighting(s)
% Draws a shadow and adds light effects.
% PRELIMINARY IMPLEMENTATION (Parameters may change)
% Parameters:
%   s : Array of plot object handles
svgBoundingBox(s, 'axes', 12, 'off');
svgGaussianBlur(s, 'SourceAlpha', 2, 'blur');
svgSpecularLightingDistant(s, 'blur', 1, 16, 2, 225, 45, 'lighting')
svgComposite(s, 'lighting', 'SourceGraphic', 'atop', 'obj');
svgGaussianBlur(s, 'SourceAlpha', 5, 'blur2');
svgOffset(s, 'blur2', [8 7], 'shade');
svgComposite(s, 'obj', 'shade', 'over', 'final');

function drop_shadow_lighting_background(s, background)
% Draws a shadow and adds light effects. The edge and face color of the
% elements is replaced by a background pixel graphic that is scaled to
% cover the axes
% PRELIMINARY IMPLEMENTATION (Parameters may change)
% Parameters:
%   s : Array of plot object handles
svgBoundingBox(s, 'axes', 12, 'off');
svgGaussianBlur(s, 'SourceAlpha', 2, 'blur');
svgSpecularLightingDistant(s, 'blur', 1, 16, 2, 225, 45, 'lighting')
svgImage(s, background, 'none', 'pic');
svgComposite(s, 'pic', 'SourceGraphic', 'atop', 'cut_pic');
svgComposite(s, 'lighting', 'cut_pic', 'atop', 'obj');
svgGaussianBlur(s, 'SourceAlpha', 5, 'blur2');
svgOffset(s, 'blur2', [8 7], 'shade');
svgComposite(s, 'obj', 'shade', 'over', 'final');


function svgBoundingBox(s, type, overlap, visible)
% Configures the bounding box of a SVG filter
% PRELIMINARY IMPLEMENTATION (Parameters may change)
% Parameters:
%   s : Array of plot object handles
%   type : [axes, element, relative]
%          Sets the filter bounding box to cover the axis reagion (axes), the
%          element extension (element or relative). Axes gives usually the
%          best results but may be slower.
%   overlap : Many filters need an overlap to work correctly.
%             Typical values for type 'axes' and 'element' -> 10
%             Typical values for type 'relative' -> 0.1
%   visible : Debugging functionality to see the bounding box used for an
%             object
for i = 1:length(s)
    userdata = get(s(i),'UserData');
    userdata.svg.BoundingBox.Visible = visible;    % Useful for debugging of bounding box for filters
    userdata.svg.BoundingBox.Type = type;          % [axes, element, relative]
    userdata.svg.BoundingBox.Overlap = overlap;
    set(s(i),'UserData', userdata);
end


function svgSpecularLightingDistant(s, source, specularConstant, specularExponent, surfaceScale, azimuth, elevation, result)
% Adds a feSpecularLighting SVG filter with distant light source
% PRELIMINARY IMPLEMENTATION (Parameters may change)
% Parameters:
%   s : Array of plot object handles
%   source : Any previous defined filter result string, 'SourceGraphic',
%            or 'SourceAlpha'.
%   specularConstant : Specular constant
%   specularExponent : Specular exponent
%   surfaceScale : Surface scaling factor
%   azimuth : Light azimuth angle [deg], typical 225.
%   elevation : Light elevation angle [deg], typical 45.
%   result : String that identifies the filter result for following filter
%            stages.
for i = 1:length(s)
    userdata = get(s(i),'UserData');
    if isfield(userdata, 'svg') && isfield(userdata.svg, 'Filter')
        next = length(userdata.svg.Filter) + 1;
    else
        next = 1;
    end
    userdata.svg.Filter(next).Subfilter.Type = 'feSpecularLighting';
    userdata.svg.Filter(next).Subfilter.Source = source;
    userdata.svg.Filter(next).Subfilter.Result = result;
    userdata.svg.Filter(next).Subfilter.SpecularConstant = specularConstant; 
    userdata.svg.Filter(next).Subfilter.SpecularExponent = specularExponent;
    userdata.svg.Filter(next).Subfilter.SurfaceScale = surfaceScale;
    userdata.svg.Filter(next).Subfilter.LightType = 'feDistantLight';
    userdata.svg.Filter(next).Subfilter.Azimuth = azimuth;
    userdata.svg.Filter(next).Subfilter.Elevation = elevation;
    set(s(i),'UserData', userdata);
end

function svgImage(s, file, aspectRatio, result)
% Adds a feImage SVG filter
% PRELIMINARY IMPLEMENTATION (Parameters may change)
% Parameters:
%   s : Array of plot object handles
%   file : Pixel graphics file name (png or jpeg) with extension.
%   aspectRatio: 'none' -> scale to bounding box limits
%                'xMinYMin meet', 'xMinYMin slice', 'xMidYMid meet', ...
%                -> see SVG 1.1 specification
%   result : String that identifies the filter result for following filter
%            stages.   
for i = 1:length(s)
    userdata = get(s(i),'UserData');
    if isfield(userdata, 'svg') && isfield(userdata.svg, 'Filter')
        next = length(userdata.svg.Filter) + 1;
    else
        next = 1;
    end
    userdata.svg.Filter(next).Subfilter.Type = 'feImage';
    userdata.svg.Filter(next).Subfilter.File = file;
    userdata.svg.Filter(next).Subfilter.AspectRatio = aspectRatio;
    userdata.svg.Filter(next).Subfilter.Result = result;
    set(s(i),'UserData', userdata);
end

function svgGaussianBlur(s, source, deviation, result)
% Adds a feGaussianBlur SVG filter
% PRELIMINARY IMPLEMENTATION (Parameters may change)
% Parameters:
%   s : Array of plot object handles
%   source : Any previous defined filter result string, 'SourceGraphic',
%            or 'SourceAlpha'.
%   deviation : Blur strength
%   result : String that identifies the filter result for following filter
%            stages.   
for i = 1:length(s)
    userdata = get(s(i),'UserData');
    if isfield(userdata, 'svg') && isfield(userdata.svg, 'Filter')
        next = length(userdata.svg.Filter) + 1;
    else
        next = 1;
    end
    userdata.svg.Filter(next).Subfilter.Type = 'feGaussianBlur';
    userdata.svg.Filter(next).Subfilter.Deviation = deviation;
    userdata.svg.Filter(next).Subfilter.Source = source;
    userdata.svg.Filter(next).Subfilter.Result = result;
    set(s(i),'UserData', userdata);
end


function svgOffset(s, source, offset, result)
% Adds a feOffset SVG filter
% PRELIMINARY IMPLEMENTATION (Parameters may change)
% Parameters:
%   s : Array of plot object handles
%   source : Any previous defined filter result string, 'SourceGraphic',
%            or 'SourceAlpha'.
%   offset : Offset value [x y]
%   result : String that identifies the filter result for following filter
%            stages.   
for i = 1:length(s)
    userdata = get(s(i),'UserData');
    if isfield(userdata, 'svg') && isfield(userdata.svg, 'Filter')
        next = length(userdata.svg.Filter) + 1;
    else
        next = 1;
    end
    userdata.svg.Filter(next).Subfilter.Type = 'feOffset';
    userdata.svg.Filter(next).Subfilter.Source = source;
    userdata.svg.Filter(next).Subfilter.Offset = offset;
    userdata.svg.Filter(next).Subfilter.Result = result;
    set(s(i),'UserData', userdata);
end

function svgComposite(s, source1, source2, operator, result)
% Adds a feComposite SVG filter
% PRELIMINARY IMPLEMENTATION (Parameters may change)
% Parameters:
%   s : Array of plot object handles
%   source1 : Any previous defined filter result string, 'SourceGraphic',
%             or 'SourceAlpha'.
%   source2 : Any previous defined filter result string, 'SourceGraphic',
%             or 'SourceAlpha'.
%   operator : Operator 'over','in','out','atop','xor','arithmetic'
%              -> see SVG 1.1 specification. 'arithmetic' is not yet
%              supported.
%   result : String that identifies the filter result for following filter
%            stages.   
for i = 1:length(s)
    userdata = get(s(i),'UserData');
    if isfield(userdata, 'svg') && isfield(userdata.svg, 'Filter')
        next = length(userdata.svg.Filter) + 1;
    else
        next = 1;
    end
    userdata.svg.Filter(next).Subfilter.Type = 'feComposite';
    userdata.svg.Filter(next).Subfilter.Source1 = source1;
    userdata.svg.Filter(next).Subfilter.Source2 = source2;
    userdata.svg.Filter(next).Subfilter.Operator = operator;
    userdata.svg.Filter(next).Subfilter.Result = result;
    set(s(i),'UserData', userdata);
end