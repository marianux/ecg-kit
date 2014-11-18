function demo_svg_water
% Create a new figure and remember the figure handle
fig = figure;
% Set default font and font size
set(fig, 'DefaultAxesFontName', 'Arial')
set(fig, 'DefaultAxesFontSize', 16)
% This is the data [year water_consumption]
data = [
1980	230
1982	265
1984	265
1986	266
1988	262
1990	259
1992	277
1994	247
1996	238
1998	250
2000	250
2002	235
2004	233
2006	228
];
% Let's plot the data
hold on
s = bar(data(:, 1), data(:, 2));
% We do not want to apply the filter to the edge of the bars. Therefore, we
% plot them in addition on top
h = bar(data(:, 1), data(:, 2));
set(h, 'FaceColor', 'none');
set(h, 'EdgeColor', 'black');
% Add the text inside each bar with correct rotation and alignment
t = text(data(:, 1), data(:, 2) - 10, num2str(data(:, 2)));
set(t, 'Rotation', -90);
set(t, 'FontSize', 16);
set(t, 'FontName', 'Arial');
set(t, 'VerticalAlignment', 'middle');
set(t, 'FontWeight', 'bold');
% Adding labels and tick marks
axis([1979 2008 0 400])
set(gca, 'XTick', data(:, 1));
set(gca, 'XTickLabel', num2str(data(:, 1)));
title('Water Consumption per Person and Year')
ylabel('Water consumption [l]')
grid on
box on
set(gca, 'Position', [0.13 0.17 0.80 0.72]);
% Now we add the filters
% The bounding box with extension axes makes sure that we cover the full
% axis region with the background images. Due to the shadow we have to
% define an overlap region of 12px. Otherwise, distortions at the border of
% the axis reagion may be visible.
svgBoundingBox(s, 'axes', 12, 'off');
% Now we add the image. To keep the aspect ratio we use 'xMidYMid slice'.
% This setting centers the picture for both x and y. The kewyord 'slice'
% makes sure that the picture is scaled so that the full axis region is
% covered.
svgImage(s, 'water_stones.jpg', 'xMidYMid slice', 'pic');
% Create a composite between bars and picture by union.
svgComposite(s, 'pic', 'SourceGraphic', 'atop', 'obj');
% To create the shadow we use the bars and use a blur filter
svgGaussianBlur(s, 'SourceAlpha', 5, 'blur');
% The blurred bars are now shifted by [10 10].
svgOffset(s, 'blur', [10 10], 'shade');
% Combine the shadow and 'picture bars' by applying the later on top of the
% shadow.
svgComposite(s, 'obj', 'shade', 'over', 'front');
% For the background we use the same picture but make it brighter
svgComposite(s, 'pic', 'pic', 'arithmetic', 'background', [0 0.2 0 0]);
% As last step, we combine the foreground and background
svgComposite(s, 'front', 'background', 'over', 'final');
% Beautify the tick mark labels by rotating them by 90 deg
setting.svg.XTickLabelAngle = -90;
set(gca, 'UserData', setting);
% Finally, we save the result into a SVG file
plot2svg('demo_svg_water.svg')