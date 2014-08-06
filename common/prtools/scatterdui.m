% SCATTERDUI Scatter plot with user interactivity
%
%   SCATTERDUI(A)
%   SCATTERDUI(A,DIM,S,CMAP,FONTSIZE,'label','both','legend','gridded')
%
% INPUT
%    DATA  Dataset
%   ...    See SCATTERD
%
% OUTPUT
%
% DESCRIPTION
% SCATTERDUI is a wrapper around SCATTERD (see SCATTERD for the options). If
% the user clicks on a sample in the plot, the corresponding index in a
% dataset is written nearby. A right-button click clears all printed indices. 
% Buttons along axes allow for browsing through the dataset dimensions. 
% Selected points are remembered when the plotted dimension changes.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% SCATTERD 

% This script is based on the segmentgui of Cris Luengo 
% <cris@ph.tn.tudelft.nl>.
% Copyright: Pavel Paclik, pavel@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: scatterdui.m,v 1.3 2006/10/23 10:21:52 davidt Exp $

function fig_hnd=scatterdui (varargin)

	% First argument is dataset.
	a = varargin{1};

	% Place all remaining arguments in a string, for the call to SCATTERD.
	%args = '';
  %for i = 2:size(varargin,2)
    %eval(sprintf('p%d = varargin{%d};',i,i)); args = [args sprintf(',p%d',i)];
	%end
  args = varargin(2:end);

	% Are we called through one of the callback handles?
	if (ischar(a))
 		switch (a)
		 case 'denclick_scatter'
		  scatterdui_inspect;
		 case 'scatterdui_change_dim'
		  scatterdui_change_dim(varargin{2});
		end
		return
	end

	% Get figure handle, clear window and start new axis.
	%	fig_hnd = gcf; clf; ui_data.axis = subplot(1,1,1); 

        % open a new figure
        figure; fig_hnd = gcf; ui_data.axis = subplot(1,1,1); 

	% tag this figure
	set(fig_hnd,'tag','scatterdui');

	set(fig_hnd,'busyaction','cancel','DoubleBuffer','on');
	
	% Call SCATTERD.
	ui_data.args = args; %eval(['scatterd(a' args,');']);
  scatterd(a, args{:});

	% Store the data for later reference.
	ui_data.a = a;
	
	% Define viable neighborhood for sample ID back-reading.
	range_x = get(ui_data.axis,'xlim'); range_y = get(ui_data.axis,'ylim');
	ui_data.neighborhood = [0.01*abs(range_x(1)-range_x(end)) ...
		    0.01*abs(range_y(1)-range_y(end))];
	ui_data.point_id = []; ui_data.text_hnd = [];

	% Add callback function to axis.
	set(ui_data.axis,'ButtonDownFcn','scatterdui(''denclick_scatter'')');
	set(get(ui_data.axis,'Children'),'ButtonDownFcn','scatterdui(''denclick_scatter'')');

	% Initialise current X- and Y-dimensions.
	ui_data.xdim = 1; ui_data.ydim = 2;

	% Place buttons for increasing/decreasing the dimensions (features) of 
	% the dataset shown on the X- and Y-axes of the plot.
	pos = get(fig_hnd,'position');
	ui_data.x_dim_inc = uicontrol('style','pushbutton', ...
				      'string','->', ...
				      'units','normalized', ...
				      'position',[0.9 0.02 0.05 0.05], ...
				      'callback','scatterdui(''scatterdui_change_dim'',''xinc'')');
	ui_data.x_dim_dec = uicontrol('style','pushbutton', ...
				      'string','<-', ...
				      'units','normalized', ...
				      'position',[0.85 0.02 0.05 0.05], ...
				      'callback','scatterdui(''scatterdui_change_dim'',''xdec'')');

	ui_data.y_dim_inc = uicontrol('style','pushbutton', ...
				      'string','^', ...
				      'units','normalized', ...
				      'position',[0.02 0.9 0.05 0.05], ...
				      'callback','scatterdui(''scatterdui_change_dim'',''yinc'')');
	ui_data.y_dim_dec = uicontrol('style','pushbutton', ...
				      'string','v', ...
				      'units','normalized', ...
				      'position',[0.02 0.85 0.05 0.05], ...
				      'callback','scatterdui(''scatterdui_change_dim'',''ydec'')');

	% Draw feature labels.
	featlabs = getfeat(ui_data.a);
	%RD Take care that we have proper char labels
	if isempty(featlabs)
		featlabs = num2str([1:size(ui_data.a,2)]');
	elseif (isnumeric(featlabs))
		featlabs = num2str(featlabs); 
	elseif iscellstr(featlabs)
		featlabs = char(featlabs);
	end;
	%RD and store them in the dataset
  ui_data.a = setfeatlab(ui_data.a,featlabs);
	xlabel(sprintf('%d : %s',ui_data.xdim,featlabs(ui_data.xdim,:)));
	ylabel(sprintf('%d : %s',ui_data.ydim,featlabs(ui_data.ydim,:)));
  
  % Save user interface data into figure window. First clear it, to
	% avoid a (sometimes) very slow update.
	set(fig_hnd,'userdata',[]); set(fig_hnd,'userdata',ui_data); 
	hold on;

	if nargout == 0
		clear('fig_hnd');
	end
	
	return

% SCATTERDUI_CHANGE_DIM (CODE)
%
% Call-back for the buttons in the window: increases or decreases the
% feature number (dimension) the plot represents on the x- or y-axis. CODE
% can be 'xinc', 'xdec', 'yinc', 'ydec'.

function scatterdui_change_dim (code)

	fig_hnd = gcbf; ui_data = get(fig_hnd,'userdata');

	[m,k,c] = getsize(ui_data.a);

	% Increase/decrease feature shown in X- or Y-axis. Loop around:
	% feature k+1 -> 1, features 1-1 -> k.

	switch (code)
	 case 'xinc'
	  ui_data.xdim = ui_data.xdim + 1;
	  if (ui_data.xdim > k)   ui_data.xdim = 1; end		
	 case 'xdec'
	  ui_data.xdim = ui_data.xdim - 1;
	  if (ui_data.xdim == 0), ui_data.xdim = k; end
	 case 'yinc'
	  ui_data.ydim = ui_data.ydim + 1;
	  if (ui_data.ydim > k),  ui_data.ydim = 1; end
	 case 'ydec'
	  ui_data.ydim = ui_data.ydim - 1;
	  if (ui_data.ydim == 0), ui_data.ydim = k; end
	end

	% Redraw figure.
	cla; ui_data.axis = subplot(1,1,1);
	% I had to add this to make it work!!!:
	scatterd(ui_data.a(:,[1 1]));
	
	%eval(['scatterd(ui_data.a(:,[ui_data.xdim,ui_data.ydim])' ui_data.args,');']);
  scatterd(ui_data.a(:,[ui_data.xdim,ui_data.ydim]), ui_data.args{:});

	% Draw feature labels.
	featlabs = getfeat(ui_data.a);
	if (isnumeric(featlabs)), featlabs = num2str(featlabs); end;

%PP!! deal with cell labels here: I don't know how, yet
	if (iscellstr(featlabs))
	   featlabs = (1:size(featlabs,1))'; 
	   featlabs = num2str(featlabs);
	end
	xlabel(sprintf('%d : %s',ui_data.xdim,featlabs(ui_data.xdim,:)));
	ylabel(sprintf('%d : %s',ui_data.ydim,featlabs(ui_data.ydim,:)));

	% Define viable neighborhood for sample ID back-reading.
	range_x = get(ui_data.axis,'xlim'); range_y = get(ui_data.axis,'ylim');
	ui_data.neighborhood = [0.01*abs(range_x(1)-range_x(end)) ...
		    0.01*abs(range_y(1)-range_y(end))];
	
 	% Redraw selected points.
	ui_data.text_hnd = ...
	    scatterdui_add_labels(ui_data.point_id, ...
				  +ui_data.a(:,[ui_data.xdim,ui_data.ydim]), ...
				  ui_data.neighborhood );
	
	% Add callback function to axis.
	set(ui_data.axis,'ButtonDownFcn','scatterdui(''denclick_scatter'')');
	set(get(ui_data.axis,'Children'),'ButtonDownFcn',...
			  'scatterdui(''denclick_scatter'')');
	
	% Save user interface data into figure window.   
	set(fig_hnd,'userdata',[]); set(fig_hnd,'userdata',ui_data); 

	return

% SCATTERDUI_INSPECT
%
% Callback for mouse-click in axis. Finds all samples in the neighbourhood
% of the clicked point and plots them, with text labels containing the index.
% Right-click clears all selected points.

function scatterdui_inspect

	fig_hnd = gcbf; ui_data = get(fig_hnd,'userdata');

	if (strcmp(get(fig_hnd,'SelectionType'),'alt'))

		% Clear all selected points.
		delete(ui_data.text_hnd); ui_data.text_hnd = []; ui_data.point_id = [];

	else

		point = get(ui_data.axis,'CurrentPoint');

		% Get all points close to the selected point from the dataset.
		% 'Close' means inside a box around POINT defined by UI_DATA.NEIGHBORHOOD.
		a = +ui_data.a; a = a(:,[ui_data.xdim,ui_data.ydim]);
		ind = find((a(:,1) >= (point(1) - ui_data.neighborhood(1))) & ...
			   (a(:,1) <= (point(1) + ui_data.neighborhood(1))) & ...
			   (a(:,2) >= (point(3) - ui_data.neighborhood(2))) & ...
			   (a(:,2) <= (point(3) + ui_data.neighborhood(2))));

		% If any points fall inside the box, plot them and add them (and their
		% text handles) to the user interface data.
		if (length(ind) > 0)
			text_hnd = scatterdui_add_labels(ind,a,ui_data.neighborhood);
			ui_data.point_id = [ui_data.point_id ind'];
			ui_data.text_hnd = [ ui_data.text_hnd text_hnd ];
		end
	end

	% Save user interface data into figure window.   
	set(fig_hnd,'userdata',[]); set(fig_hnd,'userdata',ui_data); 

	return

% HND = SCATTERDUI_ADD_LABELS (IND,DATA,NEIGHBORHOOD)
%
% Plots samples with indices IND in DATA, and places the indices as text
% labels next to the points (at a (x,y)-distance defined by NEIGHBORHOOD).
% Returns handles to all text labels in HND.

function hnd = scatterdui_add_labels (ind,data,neighborhood)

	hold on;

	% Plot the data points, plus their index in the dataset as text.
	hnd = plot(data(ind,1),data(ind,2),'gh');
	for i = 1:length(ind)
		hnd=[hnd text(data(ind(i),1) + neighborhood(1), ...
			      data(ind(i),2) + neighborhood(2), ...
			      sprintf('%d',ind(i)))];
	end

	% Add callback function to each of the texts.
	set(hnd,'ButtonDownFcn','scatterdui(''denclick_scatter'')');

	return
