%NEWFIG Create new figure on given position
%
%	NEWFIG(FIGURE_NUMBER,FIGURE_PER_ROW)
%
% INPUT
%   FIGURE_NUMBER   Number of the figure
%   FIGURE_PER_ROW  Figures per row (default: 4)
%
% OUTPUT
%
% DESCRIPTION
% Creates figure number FIGURE_NUMBER and places it on the screen,
% such that (when sufficient figures are created), the screen is
% covered by an array of figures, with FIGURE_PER_ROW figures per row.
%
% SEE ALSO
% SCATTERD, PLOTC

% $Id: newfig.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function newfig(fign,n);
	
	if nargin < 2, n=4; end
	if nargin < 1, fign=gcf+1; end
	
	% if the figure already exist, kill it!
	if any(get(0,'children') == fign)
		delete(fign);
	end
	
	% Create the figure,
	figure(fign);
	set(gcf,'menubar','none');
	% and place it in the correct position:
	ny = 1-ceil(fign/n)/n;
	nx = (fign -1 - n*floor((fign-1)/n))/n;
	d = 1/n;
	set(gcf,'units','normalized','position',[nx ny d*0.95 d*0.85]);

return
	
