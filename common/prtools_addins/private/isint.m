%ISINT Test on number(s) on integer >= 0
% $Id: isint.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function n = isint(m)

		
	if all(all(m == round(m) & m >= 0.5))
		n = 1;
	else
		n = 0;
	end

	return
