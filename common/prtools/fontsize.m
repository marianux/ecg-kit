%FONT_SIZE Set large graphic font
%
%       font_size(size)
%
% Set font size for current figure

function font_size(size)

  V = axis;
  H = get(gcf,'Children');
	c1 = [];
	for h = H'
	if strcmp(get(h,'type'),'axes')
  	set(get(h,'XLabel'), 'FontSize', size);
  	set(get(h,'YLabel'), 'FontSize', size);
  	set(get(h,'Title'),  'FontSize', size);
  	set(h, 'FontSize', size);
  	c1 = [c1; get(gca, 'Children')];
	end
	end
 	axis(V);
  for h1 = c1'
    v1 = get (h1);
    if (isfield (v1, 'FontSize'))
      set (h1, 'FontSize', size);
    end;
    c2 = get (h1, 'Children');
    for h2 = c2'
      v2 = get (h2);
      if (isfield (v2, 'FontSize'))
        set (h2, 'FontSize', size);
      end;
    end;
  end;

return

