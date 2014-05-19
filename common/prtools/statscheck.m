%STATSCHECK  Check whether the STATS toolbox is in the path

function n = statscheck
	
  n = ~isempty(ver('stats'));
  if (nargout == 0) & (n == 0)
		error([newline 'The Matlab STATS toolbox is missing.'])
  end
		
return