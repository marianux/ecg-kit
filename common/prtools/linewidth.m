%LINEWIDTH Set linewidth in plot
%
%	linewidth(width)
%Set linewidth for current figure

function linewidth(width)
if strcmp(get(gca,'type'),'line')
	set(gca,'linewidth',width);
end
children = get(gca,'children');
set_linewidth_children(children,width)
return

function set_linewidth_children(children,width)
if isempty(children), return, end
for i = children(:)'
	if length(i) > 1
		set_linewidth_children(i,width)
		return
	end
	if strcmp(get(i,'type'),'line') | strcmp(get(i,'type'),'patch')
		set(i,'linewidth',width);
	end
	children2 = get(i,'children');
	set_linewidth_children(children2,width)
end
