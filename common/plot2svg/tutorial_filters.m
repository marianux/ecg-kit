function tutorial_filters
file = 'matterhorn_small.jpg';
fig = figure;
set(fig, 'DefaultAxesFontName', 'Arial')
set(fig, 'DefaultAxesFontSize', 6)
plotCounter = 1;
xPlots = 3;
yPlots = 4;
% a) Original plot with a line, patch, and text object
subplot(yPlots, xPlots, plotCounter)
plotCounter = plotCounter + 1;
h = create_plot;
title('a) Original')
% b) Adding a Gaussian blur filter for each object that uses the standard
%    color channels for the blur operation. The bounding box is covering
%    the axis region for each object. Note: there is no bounding box
%    overlap specified. This would lead to some distortion at the border if
%    a offset filter followed the blur filter.
subplot(yPlots, xPlots, plotCounter)
plotCounter = plotCounter + 1;
h = create_plot;
svgBoundingBox(h, 'axes', 0, 'off')
svgGaussianBlur(h, 'SourceGraphic', 2, 'blur');
title('b) Blur on (S = SourceGraphic)')
% c) Identical to b) but the blur filter uses the alpha channel. This mode
%    is very useful to create shadows.
subplot(yPlots, xPlots, plotCounter)
plotCounter = plotCounter + 1;
h = create_plot;
svgBoundingBox(h, 'axes', 0, 'off')
svgGaussianBlur(h, 'SourceAlpha', 2, 'blur');
title('c) Blur (S = SourceAlpha)')
% d) Adding a image filter that covers the whole axis region. The filter
%    replaces the object as no other combine filter is defined. The image
%    is scaled to the bounding box. This leads to a distortion of the
%    aspect ratio.
subplot(yPlots, xPlots, plotCounter)
plotCounter = plotCounter + 1;
h = create_plot;
svgBoundingBox(h, 'axes', 0, 'on')
svgImage(h, file, 'none', 'pic');
title('d) Pict (BB = axes, none)')
% e) Identical to d) but with correct aspect ratio of the image. The
%    setting 'xMidYMid slice' centers the image and scales x and y so that
%    both cover the bonding box region defined by the axis region. The
%    clipping of the axis object removes the overlap.
subplot(yPlots, xPlots, plotCounter)
plotCounter = plotCounter + 1;
h = create_plot;
svgBoundingBox(h, 'axes', 0, 'on')
svgImage(h, file, 'xMidYMid slice', 'pic');
title('e) Pict (BB = axes, xMidYMid slice)')
% f) Identical to e) but with image scaling 'xMidYMid meet'. This setting
%    also conserves the aspect ratio. However, the image may no more cover
%    the whole bounding box. 
subplot(yPlots, xPlots, plotCounter)
plotCounter = plotCounter + 1;
h = create_plot;
svgBoundingBox(h, 'axes', 0, 'on')
svgImage(h, file, 'xMidYMid meet', 'pic');
title('f) Pict (BB = axes, xMidYMid meet)')
% h) Identical to f) but with bounding box defined by the object. The
%    bounding box may be larger than the object extension due to line
%    width and marker size.
subplot(yPlots, xPlots, plotCounter)
plotCounter = plotCounter + 1;
h = create_plot;
svgBoundingBox(h, 'element', 0, 'on')
svgImage(h, file, 'xMidYMid meet', 'pic');
title('g) Pict (BB = element, xMidYMid meet)')
% g) Now we add a composite filter that combines the filter result from the
%    image filter and the source graphic. The keyword 'atop' restricts the
%    image filter to the object boundaries.
subplot(yPlots, xPlots, plotCounter)
plotCounter = plotCounter + 1;
h = create_plot;
svgBoundingBox(h, 'axes', 0, 'on')
svgImage(h, file, 'xMidYMid slice', 'pic');
svgComposite(h, 'pic', 'SourceGraphic', 'atop', 'obj');
title('h) Pict (BB = axes, xMidYMid slice)')
% i) Identical to g) but with a filter bounding box defined by the object
%    extension. The image scaling 'xMidYMid meet' is not well suited for
%    this application. The image is not covering the whole object.
subplot(yPlots, xPlots, plotCounter)
plotCounter = plotCounter + 1;
h = create_plot;
svgBoundingBox(h, 'element', 0, 'on')
svgImage(h, file, 'xMidYMid meet', 'pic');
svgComposite(h, 'pic', 'SourceGraphic', 'atop', 'obj');
title('i) Pict (BB = element, xMidYMid meet)')
% j) Identical to i) but with a better image scaling that conserves the
%    spect ratio and makes sure that the whole bounding box is covered.
subplot(yPlots, xPlots, plotCounter)
plotCounter = plotCounter + 1;
h = create_plot;
svgBoundingBox(h, 'element', 0, 'on')
svgImage(h, file, 'xMidYMid slice', 'pic');
svgComposite(h, 'pic', 'SourceGraphic', 'atop', 'obj');
title('j) Pict (BB = element, xMidYMid slice)')
% k) Identical to j) but with a bounding box overlap of 5 pixel. This is
%    important to avoid distortions at the border. 
subplot(yPlots, xPlots, plotCounter)
plotCounter = plotCounter + 1;
h = create_plot;
svgBoundingBox(h, 'element', 5, 'on')
svgImage(h, file, 'xMidYMid slice', 'pic');
svgComposite(h, 'pic', 'SourceGraphic', 'atop', 'obj');
title('k) Pict (BB = element, xMidYMid slice)')
% l) Identical to k) but the composition is based on the image filter and a
%    blur filter. In addition, the image filter covers the whole axis
%    region.
subplot(yPlots, xPlots, plotCounter)
plotCounter = plotCounter + 1;
h = create_plot;
svgBoundingBox(h, 'axes', 0, 'on')
svgImage(h, file, 'xMidYMid slice', 'pic');
svgGaussianBlur(h, 'SourceAlpha', 2, 'blur');
svgComposite(h, 'pic', 'blur', 'atop', 'obj');
title('l) Blur + Pict (BB = axes, xMidYMid slice)')
% Save the result
plot2svg('tutorial_filters.svg')

function handles = create_plot
% Helper function to plot the objects
hold on
s1 = text(3, 1, 'plot2svg', 'FontName', 'Arial', 'FontSize', 16, ...
    'Color', 'red', 'FontWeight', 'bold', 'Margin', 0.01, 'Clipping', 'on');
s2 = plot([6 8], [3 6], 'b', 'LineWidth', 20);
s3 = patch([1 2 3 4 3 2 1], [1 3 2 3 4 5 3],'black');
handles = [s1 s2 s3];
box on 
grid on
axis([0 10 0 6])
