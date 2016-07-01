%% (Internal) Plot the ROC curve
%
%       fig_hdl = plot_roc( roc_array, strTitle, strLegend, bPrint, fig_hdl )
% 
% 
% Arguments:
% 
%      + roc_array: string to parse
%             
%      + strTitle: string to parse
%             
%      + strLegend: string to parse
%             
%      + bPrint: string to parse
%             
%      + fig_hdl: string to parse
% 
% Output:
% 
%      + fig_hdl : 
% 
% Example:
% 
% 
% See also ECGwrapper
% 
% Author: Mariano Llamedo Soria (llamedom at frba.utn.edu.ar)
% Version: 0.1 beta
% Birthdate  : 30/7/2014
% Last update: 30/7/2014
% Copyright 2008-2015
% 
function fig_hdl = plot_roc( roc_array, strTitle, strLegend, bPrint, fig_hdl )

if( nargin < 3  ) 
    strLegend = [];
end

if( nargin < 2  ) 
    strTitle = [];
end

if( nargin < 4 || isempty(bPrint) ) 
    bPrint = false;
end

if( nargin < 5 || isempty(fig_hdl) ) 
    fig_hdl = figure();
else
    fig_hdl = figure(fig_hdl);
%     clf(fig_hdl );
end

if( isstruct(roc_array) )
    roc_array = {roc_array};
elseif(~iscell(roc_array))
    error('roc_array should be a struct, output of prtools roc() function; or a cell containing several structs.')
end

rocs2plot = length(roc_array);
colors = my_colormap(max(10,rocs2plot));
Xvalues = [ ];
Yvalues = [ ];

for ii = 1:rocs2plot
    this_roc = roc_array{ii};
    Xvalues = [ Xvalues 1-colvec(this_roc.specificity) ];
    Yvalues = [ Yvalues colvec(this_roc.sensitivity) ];
end

axes_hdl = gca;
set(axes_hdl, 'ColorOrder', colors);
hold(axes_hdl, 'on')

plot(axes_hdl, Xvalues, Yvalues, 'Linewidth', 1.5)

set(axes_hdl, 'Box', 'off' );
set(axes_hdl, 'Xtick', 0:0.1:1 );
set(axes_hdl, 'XtickLabel', num2str(colvec(100:-10:0)) );
set(axes_hdl, 'Ytick', 0.1:0.1:1 );
set(axes_hdl, 'YtickLabel', num2str(colvec(10:10:100)) );

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

set(axes_hdl, 'FontName', strFontName );
set(axes_hdl, 'FontSize', FontSizeTick );


title(strTitle, 'FontName', strFontName, 'FontSize', FontSizeTitle)
legend(strLegend, 'Location', 'East', 'FontName', strFontName, 'FontSize', FontSizeTick)

ylabel('Sensitivity', 'FontName', strFontName, 'FontSize', FontSizeLabel)
xlabel('Specificity', 'FontName', strFontName, 'FontSize', FontSizeLabel)
