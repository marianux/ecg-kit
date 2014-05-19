%MARKSIZE Change markersize in lineplot or scatterplot
%
%	marksize(size,marker)
%
% The marker fontsize (default set by matlab to 6) is reset
% to size for the given marker, default: all markers.
% 
function marksize(siz,marker)
if nargin == 1, marker = ' '; end
h = get(gca,'Children')';
for i = h
	if strcmp(get(i,'Type'),'line')
		if nargin == 1 | get(i,'Marker') == marker
			set(i,'MarkerSize',siz);
		end
	end
end
