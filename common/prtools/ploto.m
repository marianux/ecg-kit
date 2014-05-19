%PLOTO Plot objects as 1-D functions of the feature number
% 
%   [HH HO HC] = PLOTO(A,N)
%
% INPUT
%   A   Dataset
%   N   Integer
%
% OUTPUT
%   HH  Lines handles
%   HO  Object identifier handles
%   HC  Class number handles
% 
% DESCRIPTION
% Produces 1-D function plots for all the objects in dataset A. The plots
% are organised as subplots, N on a row. Default is the squareroot of the
% number of objects. Object identifiers and class numbers are written in
% the correspopnding plots.
%
% See also DATASETS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [h_out1,h_out2,h_out3] = ploto(a,p)

		
	if nargin < 2, p = []; end
	[m,k,c] = getsize(a);
	nlab = getnlab(a);

	% Define the color for each of the classes:
	if c == 2
		map = [0 0 1; 1 0 0];
	else
		map = hsv(c);
	end

	% Make subplots for each object, so a grid of p x q subplots is
	% defined
	h = [];
	if ~isempty(p)
		q = ceil(m/p);
	elseif m > 3
		p = ceil(sqrt(m)); q = ceil(m/p);
	else
		p = m; q = 1;
	end
	% Get the object labels
	labs = getlabels(a);
	ymin = min(a.data(:));
	ymax = max(a.data(:));
	V = [1 k ymin ymax];
	% Make the plot for each of the objects:
	h = [];
	ho = [];
	hc = [];
	s = sprintf('Plot %i objects: ',m);
	prwaitbar(m,s);
	for j = 1:m
		if isdatafile(a) | 1
			prwaitbar(m,j,[s int2str(j)]);
			b = +prdataset(a(j,:));
			ymin = min(b);
			ymax = max(b);
			k = length(b);
			V = [1 k ymin ymax];
		else
			b = +a(j,:);
		end
		% Create the subplots with the correct sizes:
		subplot(q,p,j)
		hh = plot(b);
		set(gca,'xtick',[]);
		set(gca,'ytick',[]);
		axis(gca,V);
		ho = [ho text(2,ymax-0.15*(ymax-ymin),getident(a(j,:),'string'))];
		hc = [hc text(3*k/4,ymax-0.15*(ymax-ymin),num2str(nlab(j)))];
		h = [h hh];
		hold on
	end
	prwaitbar(0);
	
	% The last details to take care of:
	if nargout > 0
		h_out1 = h;
		h_out2 = ho;
		h_out3 = hc;
	end

	return
