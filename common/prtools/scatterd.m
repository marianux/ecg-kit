%SCATTERD Display scatterplot
% 
%   H = SCATTERD(A)
%   H = A*SCATTERD
%
%   H = SCATTERD(A,DIM,S,CMAP,FONTSIZE,'label','both','legend','gridded')
%   H = A*SCATTERD([],DIM,S,CMAP,FONTSIZE,'label','both','legend','gridded')
%   H = A*SCATTERD(DIM,S,CMAP,FONTSIZE,'label','both','legend','gridded')
%
% INPUT
%   A     Dataset or matrix
%   DIM   Number of dimensions: 1,2 or 3 (optional; default: 2)
%   S     String specifying the colors and markers (optional)
%   CMAP  Matrix with a color map (optional)
% 
% OUTPUT
%   H     Vector of handles
%
% DESCRIPTION
% SCATTERD(A) displays a 2D scatterplot of the first two features of the 
% dataset A. If the number of dimensions DIM is specified (1..3), it plots 
% the first D features in a D-dimensional plot (D<4). If the plot string S 
% is provided, e.g. S = 'w+', all points are plotted accordingly. If given,
% different plot strings are used for different classes. See PLOT for the
% specification of plot strings.
%   
% If CMAP is specified, the color of the object symbols is determined by 
% CMAP indexed by the object labels. A colormap has the size [C x 3], where 
% C is the number of classes. The three components of CMAP(I,:) determine 
% the red, green and blue components of the color. For instance: 
%
%   MAP = HSV; [M,K] = SIZE(A); LABELS = CEIL(64*[1:M]'/M); 
%   A = PRDATASET(A,LABELS); SCATTERD(A,'.','COLORMAP',MAP); 
%
% This may be used for tracking ordered objects.
%  
% FONTSIZE may be a vector with three elements: fontsize, markersize and
% size of the label font in case of a label plot.
%  
% Various other options are:
%   'LABEL'  : plot labels instead of symbols
%   'BOTH'   : plot labels next to each sample
%   'LEGEND' : place a legend in the figure
%   'GRIDDED': make a grid of 2D scatterplots of each pair of features
%  
% All the parameters, except for the dataset A can be specified in any 
% order or can be left out.
%  
% Classifiers can be plot in the scatterplot by PLOTC.
% Note that PLOTC does not work for 1D and 3D scatterplots.
%
% EXAMPLES
% See PREX_CONFMAT, PREX_DENSITY, PREX_PLOTC, PREX_MCPLOT.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, PLOTC, SCATTERN

% Copyright: D. de Ridder, R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: scatterd.m,v 1.8 2010/06/25 09:52:40 duin Exp $

% REVISIONS
% DR1 - Dick, 05-10-2004
% Added plotting of unlabeled data as 'k.'.

function handle = scatterd(varargin)

  argin = shiftargin(varargin,'scalar');
  argin = shiftargin(argin,'char');
  argin = setdefaults(argin,[],[],[],[],[],[],[],[],[]);
  if mapping_task(argin,'definition')
    % standard return, name = filename
    handle = define_mapping(argin,'fixed');
    return
  end
  
  par = cell(1,8);
  [a,par{:}] = deal(argin{:});
	a = prdataset(a); 		% Allow for a non-dataset data
	a = remclass(a);
  
		% Defaults
	d = min(size(a,2),2);		% Dimensionality of plot
	s = []; 					% Plot symbol(s)
	map = []; 					% Color map
	plotlab = 0;				% Put text labels instead of or next to samples
	plotsym = 1;				% Plot symbols for samples
	plotlegend = 0; 		    % Plot legend
	gridded = 0;				% Make a gridded plot
	gridrun = 0;				% Inner loop in a gridded plot?
	font_size = [];
	mark_size = [];
	lab_size  = [];
	hold_axis = ishold; % A flag to check if 'hold on' is set for the current axis
% 
% 	if (nargin < 9), par{8} = []; else, par{8} = p8; end
% 	if (nargin < 8), par{7} = []; else, par{7} = p7; end
% 	if (nargin < 7), par{6} = []; else, par{6} = p6; end
% 	if (nargin < 6), par{5} = []; else, par{5} = p5; end
% 	if (nargin < 5), par{4} = []; else, par{4} = p4; end
% 	if (nargin < 4), par{3} = []; else, par{3} = p3; end
% 	if (nargin < 3), par{2} = []; else, par{2} = p2; end
% 	if (nargin < 2), par{1} = []; else, par{1} = p1; end

	% Set up default values.
	for i = 1:5
		if (~isempty(par{i})) 
			if (length(par{i}) == 1) & par{i} < 5 & (~ischar(par{i})) 	% Dimensionality
				d   = par{i}; par{i} = 2; 		% Necessary for gridded: D needs to be 2.
			elseif ((size(par{i},1) > 1) & (size(par{i},2)==3) & (~ischar(par{i}))) 	% Color map
				map = par{i};
			elseif (strcmp(par{i},'label'))
				plotlab = 1; plotsym = 0;
			elseif (strcmp(par{i},'both'))
				plotlab = 1; plotsym = 1;
			elseif (strcmp(par{i},'legend'))
				plotlegend = 1;
			elseif (strcmp(par{i},'gridded'))
				gridded = 1; 
				par{i} = 'gridrun'; 					% Necessary for gridded: otherwise an infinite recursion :)
			elseif (strcmp(par{i},'gridrun'))
				gridrun = 1;
			elseif ~isstr(par{i}) & length(par{i}) <= 3 & par{i}(1) >= 5
				font_size = par{i}(1); 
				if length(par{i}) >= 2
					mark_size = par{i}(2);
				end
				if length(par{i}) == 3
					lab_size = par{i}(3);
				end
			else
				s	= par{i};
			end
		end
	end

	if (gridrun)
		if isempty(font_size), font_size = 10; end
		if isempty(mark_size), mark_size =  5; end
		if isempty(lab_size),   lab_size =  8; end
	elseif (~hold_axis)
		%clf; 
		cla;
		if isempty(font_size), font_size = 16; end
		if isempty(mark_size), mark_size =  7; end
		if isempty(lab_size),   lab_size = 14; end
	else
		if isempty(font_size), font_size = get(gca,'fontsize'); end
		if isempty(mark_size), mark_size = font_size/2; end
		if isempty(lab_size),   lab_size = font_size-2; end
	end

	feats = getfeatlab(a,'string');
	
	if isempty(feats)
		for i=1:size(a,2)
			feats = strvcat(feats,sprintf('Feature %d',i));
		end
	else
		if size(feats,2) == 1
			feats = [repmat('Feature ',size(feats,1),1) feats];
		end
	end

	if (gridded)
		clf;
		gs = size(a,2);
		for i = 1:gs
			for j = 1:gs
				subplot(gs,gs,(i-1)*gs+j);
				gridrun = 1;
				h = feval(mfilename,a(:,[j i]),par{1},par{2},par{3},par{4},par{5},par{6},par{7});
				gridrun = 0;
				if (i~=gs), xlabel(''); end
				if (j~=1),  ylabel(''); end
			end
			subplot(111); 		% Takes care of clf for the next scatterplot.
		end
		return;
	end

	if (isa(a,'prdataset')) & (~isempty(getlablist(a)))
		[m,k,c] = getsize(a);
		lab = getnlab(a); 
		lablist = getlablist(a,'string');
		ident = getident(a); %DXD: needed for prcursor
		dataset_name = getname(a);
		a = double(a);
	else
		[m,k] = size(a);
		lab = ones(m,1);
		ident = (1:m)'; %DXD: needed for prcursor
		dataset_name = [];
		c = 1;
		a = double(a);
	end

	% Character string defining the plotting setup in terms of symbols and colors.
	if (isempty(s))
		vers = version;
		if (str2num(vers(1)) < 5)
			col = 'brmw';
			sym = ['+*xo.']';
			i   = [1:20];
			ss  = [col(i-floor((i-1)/4)*4)' sym(i-floor((i-1)/5)*5)];
		else
			col = 'brmk';
			sym = ['+*oxsdv^<>p']';
			i   = [1:44];
			ss  = [col(i-floor((i-1)/4)*4)' sym(i-floor((i-1)/11)*11)];
		end
    ss  = ['k.'; ss];  								% DR1 - Add symbol for "unlabeled" data.
		[ms,ns] = size(ss);
		if ms == 1, ss = setstr(ones(m,1)*ss); end
	else
		if size(s,1) == 1
			ss = repmat(s,c,1); s = [];
		else
			ss = s; s = [];
		end
    ss  = char('k.', ss);  								% DR1 - Add symbol for "unlabeled" data.
	                                                %DXD - changed [.] into char(.)
	end

	% Define some 'space' OY to be added around the data plotted in symbols.
	oy = zeros(1,3);		
	if (plotsym)
		oy = 0.02*(max(a)-min(a));
	else
		s = 'w.'; 				% Plot white spot instead of symbols.
	end
	oy(2) = 0;

	% Make a plot.
	lhandle = []; thandle = [];

  % Also plot label "0" (unlabeled).
  ss = [ss; ss; ss; ss];
	for i = 0:c
		J = find(lab==i);
		if (isempty(s)), symbol = ss(i+1,:); else, symbol = s; end
		if ((d == 1) & ~isempty(J))
			h = plot(a(J,1),zeros(length(J),1),symbol);
			hold on;
			set(h,'markersize',mark_size);
			lhandle = [lhandle h];
			if (plotlab)
				for j = 1:length(J)
					h = text(a(J(j),1)+oy(1),oy(2),lablist(lab(J(j)),:));
					set(h,'fontsize',lab_size);
					thandle = [thandle h]; if (~isempty(map)), set (h, 'color', prmap(i+1,:)); end
				end
			end
		elseif ((d == 2) & ~isempty(J))
			h = plot(a(J,1),a(J,2),symbol);
			hold on;
			set(h,'markersize',mark_size);
			lhandle = [lhandle h];
			if (plotlab)
				for j = 1:length(J)
					h = text(a(J(j),1)+oy(1),a(J(j),2)+oy(2),lablist(lab(J(j)),:));
					set(h,'fontsize',lab_size);
					thandle = [thandle h]; if (~isempty(map)), set (h, 'color', prmap(i+1,:)); end
				end
			end
		elseif (~isempty(J))
			h = plot3(a(J,1),a(J,2),a(J,3),symbol);
			hold on;
			set(h,'markersize',mark_size);
			lhandle = [lhandle h];
			if (plotlab)
				for j = 1:length(J)
					h = text(a(J(j),1)+oy(1),a(J(j),2)+oy(2),a(J(j),3)+oy(3),lablist(lab(J(j)),:));
					set(h,'fontsize',lab_size);
					thandle = [thandle h]; if (~isempty(map)), set (h, 'color', prmap(i+1,:)); end			
				end
			end
		end
		%DXD: store the object identifiers in the userdata such that you
		%can retrieve them by prcursor:
		if ~isempty(J)
			ud = get(h,'UserData');
			ud.ident = ident(J);
			set(h,'UserData',ud);
		end
	end

	if (plotsym)
		if (~isempty(map))
			for i = 0:c, set (lhandle(i+1), 'color', prmap(i+1,:)); end
		end
		if (plotlegend), 
			[ht, hl] = legend(lhandle(:),lablist);
			hl = hl(:)';
			%set(hl(1:2*c),'markersize',mark_size);
			lhandle = [lhandle hl];
			thandle = [thandle ht(:)'];
		end
	end

	% !%_%*!_% Matlab

	set(gca,'fontsize',font_size);

	if (~hold_axis)
		dd = (max(a,1) - min(a,1))*0.05;		% offset, avoiding points on plot box.
		J = find(dd==0);
		dd(J) = dd(J)+1;

%feats
		if (d == 1),     
			axis ([min(a(:,1))-dd(1) max(a(:,1))+dd(1) -0.5 0.5]); 
			hx = xlabel(feats(1,:));
			thandle = [thandle hx];
		elseif (d == 2),
			axis ([min(a(:,1))-dd(1) max(a(:,1))+dd(1) min(a(:,2))-dd(2) max(a(:,2))+dd(2)]); 
			hx = xlabel(feats(1,:)); hy = ylabel(feats(2,:));
			thandle = [thandle(:)' hx hy];
			view(2); 
		else		% D = 3            
			axis ([min(a(:,1))-dd(1) max(a(:,1))+dd(1) min(a(:,2))-dd(2) max(a(:,2))+dd(2) min(a(:,3))-dd(3) max(a(:,3))+dd(3)]); 
			hx = xlabel(feats(1,:)); hy = ylabel(feats(2,:)); hz = zlabel(feats(3,:));
			thandle = [thandle hx hy hz];
			view(3); 
		end
	end

	%if (~gridrun) & (length(get(gcf,'children')) < 3 & any(get(gca,'position')>0.80))
	if (~gridrun) & (length(get(gcf,'children')) < 3) 
		P = get(gca,'position');
		if any(P) > 0.80 & all(P(3:4)>0.50) % don't do this in case of subplots
			set(gca,'position',[0.13 0.13 0.79 0.78]); % make axis labels visible
		end
	end
		
	if (~gridrun) & (~isempty(dataset_name))
		title(dataset_name);
	end

	hold off;

	if (nargout > 0)
		handle = [lhandle thandle];
	end

	%DXD this should be standard:	
	vers = version;
	if (str2num(vers(1)) > 6)
    prcursor(gcf);
    datacursormode('off');
  end

return;

